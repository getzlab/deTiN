import numpy as np
from scipy.stats import beta

np.seterr(all='ignore')


class model:
    """Model of tumor in normal (TiN) based on only candidate SSNVs. This estimate is most
    reliable when there are greater then 6 mutations and TiN is less then ~30%. Previously
     this estimate has been used on its own to assess TiN in targeted panel data where
     copy number data is usually limited but SSNVs are well measured.
     TiN estimate : model.TiN
     Somatic classification of SSNVs : model.E_z (E_z > 0.5 -> somatic)"""

    def __init__(self, candidate_sites, p_somatic):
        # variables follow notation:
        # ac = allele count n = normal t = tumor

        # Variables for SSNV fit
        self.TiN_range = np.linspace(0, 1, num=101)
        self.af = np.linspace(0.005, 1, num=200)

        # observed data
        self.contig = candidate_sites['contig']
        self.position = candidate_sites['position']
        self.genomic_coord_x = candidate_sites['genomic_coord_x']
        self.n_alt_count = np.array(candidate_sites['n_alt_count'])
        self.n_ref_count = np.array(candidate_sites['n_ref_count'])
        self.n_depth = self.n_alt_count + self.n_ref_count
        self.normal_f = np.true_divide(self.n_alt_count, self.n_depth)
        self.t_alt_count = np.array(candidate_sites['t_alt_count'])
        self.t_ref_count = np.array(candidate_sites['t_ref_count'])
        self.t_depth = self.t_alt_count + self.t_ref_count
        self.tumor_f = np.true_divide(self.t_alt_count, self.t_depth)
        self.number_of_sites = len(self.n_alt_count)

        # hyperparameter
        self.p_somatic = p_somatic

        # parameter
        self.TiN = 0
        self.CI_tin_high = []
        self.CI_tin_low = []
        self.E_z = np.zeros([self.number_of_sites, 1])

        # expected allele fraction of minor allele given allelic copy data
        self.psi = .5 - np.array(candidate_sites['f_acs'])
        self.t_het_direction = self.tumor_f < 0.5
        self.t_het_direction = self.t_het_direction * -1
        self.t_het_direction[self.t_het_direction == 0] = 1

        # determine ratio of tumor to normal copies given tau and TiN at each locus
        self.tau = candidate_sites['tau']
        self.tin_correct_tau = np.multiply(self.TiN_range, candidate_sites['tau'][:, np.newaxis])
        self.tin_correct_normal_tau = np.multiply((1 - self.TiN_range), 2)
        self.CN_ratio = np.divide(self.tin_correct_tau, np.array(self.tin_correct_tau + self.tin_correct_normal_tau))

        # random variables
        self.rv_normal_af = beta(self.n_alt_count + 1, self.n_ref_count + 1)
        self.rv_tumor_af = beta(self.t_alt_count + 1, self.t_ref_count + 1)

        # conditionals
        self.p_TiN_given_S = np.zeros([self.number_of_sites, 101])
        self.p_TiN_given_G = np.zeros([self.number_of_sites, 101])
        self.p_TiN_given_het = np.zeros([self.number_of_sites, 101])
        self.p_artifact = np.zeros([self.number_of_sites, 1])

        # likelihood
        self.TiN_likelihood = np.zeros([101, 1])

    def generate_conditional_ps(self):
        # p(TiN|Somatic) and p(TiN|Germline)
        n_af_w = np.zeros([self.number_of_sites, len(self.af)])
        t_af_w = np.zeros([self.number_of_sites, len(self.af)])
        t_het_direction = np.ones([self.number_of_sites, len(self.af)])
        t_het_direction[:, 0:np.int(np.round(np.true_divide(len(self.af),2)))] = -1
        for i, f in enumerate(self.af):
            n_af_w[:, i] = self.rv_normal_af.cdf(f ) - self.rv_normal_af.cdf(f - 0.005)
            t_af_w[:, i] = self.rv_tumor_af.cdf(f ) - self.rv_tumor_af.cdf(f - 0.005)
            # ac given somatic
            t_af = np.multiply(f, self.n_depth)
            n_ac_given_tin = np.multiply(t_af[:, np.newaxis], self.CN_ratio)

            # ac given heterozygous
            f_t_af = 0.5 - np.abs((.5 - f))
            psi_t_af = 0.5 - f_t_af
            psi_t_af = np.multiply(psi_t_af, t_het_direction[:, i])
            exp_f = 0.5 + np.multiply(psi_t_af[:, np.newaxis], self.CN_ratio)
            n_het_ac_given_tin = np.multiply(exp_f, self.n_depth[:, np.newaxis])

            for TiN_idx, TiN in enumerate(self.TiN_range):
                self.p_TiN_given_S[:, TiN_idx] += np.multiply(beta.pdf(self.normal_f[:], n_ac_given_tin[:, TiN_idx] + 1,
                                                                       self.n_depth[:] - n_ac_given_tin[:,
                                                                                         TiN_idx] + 1) * 0.01,
                                                              t_af_w[:, i])
                self.p_TiN_given_het[:, TiN_idx] += np.multiply(
                    beta.pdf(self.normal_f[:], n_het_ac_given_tin[:, TiN_idx] + 1,
                             self.n_depth[:] - n_het_ac_given_tin[:, TiN_idx] + 1) * 0.01,
                    t_af_w[:, i])

        self.p_artifact = self.rv_tumor_af.pdf(self.normal_f) * 0.01
        self.p_TiN_given_G = np.multiply(1 - self.p_artifact[:, np.newaxis], self.p_TiN_given_het) + np.multiply(
            self.p_artifact[:, np.newaxis], 1 - self.p_TiN_given_het)
        self.p_TiN_given_G = np.true_divide(self.p_TiN_given_G,np.nansum(self.p_TiN_given_G,axis=1)[:,np.newaxis])
        self.p_TiN_given_S = np.true_divide(self.p_TiN_given_S,np.nansum(self.p_TiN_given_S,axis=1)[:,np.newaxis])
    def expectation_of_z_given_TiN(self):
        # E step
        numerator = self.p_somatic * self.p_TiN_given_S[:, self.TiN]
        denominator = numerator + np.array([1 - self.p_somatic] * self.p_TiN_given_G[:, self.TiN])
        self.E_z = np.nan_to_num(np.true_divide(numerator, denominator))

    def maximize_TiN_likelihood(self):
        # M step
        self.TiN_likelihood = np.nansum(np.multiply(self.E_z[:, np.newaxis], np.ma.log(self.p_TiN_given_S)), axis=0) + \
                              np.nansum(np.multiply(1 - self.E_z[:, np.newaxis], np.ma.log(self.p_TiN_given_G)))
        self.TiN = np.argmax(self.TiN_likelihood)

    def perform_inference(self):
        print 'pre-processing SSNV data'
        self.generate_conditional_ps()
        TiN_last = []
        iteration = 0
        print 'initialized TiN to 0'
        while self.TiN != TiN_last and iteration <= 100:
            iteration += 1
            TiN_last = self.TiN
            self.expectation_of_z_given_TiN()
            self.maximize_TiN_likelihood()
            print 'TiN inference after ' + str(iteration) + ' iterations = ' + str(self.TiN_range[self.TiN])
        print 'SSNV based TiN estimate converged: TiN = ' + str(self.TiN_range[self.TiN])
        self.TiN = self.TiN_range[self.TiN]

        posterior = np.exp(self.TiN_likelihood- np.nanmax(self.TiN_likelihood))
        self.CI_tin_low = self.TiN_range[
            next(x[0] for x in enumerate(np.cumsum(np.ma.masked_array((np.true_divide(posterior, np.nansum(posterior))))))
             if x[1] > 0.025)]
        self.CI_tin_high = self.TiN_range[
            next(x[0] for x in enumerate(np.cumsum(np.ma.masked_array((np.true_divide(posterior, np.nansum(posterior))))))
            if x[1] > 0.975)]
