import numpy as np
import pandas as pd
from scipy.stats import beta
import deTiN_utilities as du

np.seterr(all='ignore')


class model:
    """Model of tumor in normal (TiN) based on only candidate SSNVs. This estimate is most
    reliable when there are greater then 6 mutations and TiN is less then ~30%. Previously
     this estimate has been used on its own to assess TiN in targeted panel data where
     copy number data is usually limited but SSNVs are well measured.
     TiN estimate : model.TiN
     Somatic classification of SSNVs : model.E_z (E_z > 0.5 -> somatic)"""

    def __init__(self, candidate_sites, p_somatic, resolution=101, f_thresh=0.15, depth=15, hot_spots_file = 'NA'):
        # variables follow notation:
        # ac = allele count n = normal t = tumor

        # Variables for SSNV fit
        self.TiN_range = np.linspace(0, 1, num=resolution)
        self.af = np.linspace(0.005, 1, num=200)

        # observed data
        self.contig = candidate_sites['contig']
        self.position = candidate_sites['position']
        self.genomic_coord_x = candidate_sites['genomic_coord_x']
        self.n_alt_count = np.array(candidate_sites['n_alt_count'])
        self.n_ref_count = np.array(candidate_sites['n_ref_count'])
        self.n_depth = self.n_alt_count + self.n_ref_count
        self.normal_f = np.nan_to_num(np.true_divide(self.n_alt_count, self.n_depth))
        self.t_alt_count = np.array(candidate_sites['t_alt_count'])
        self.t_ref_count = np.array(candidate_sites['t_ref_count'])
        self.t_depth = self.t_alt_count + self.t_ref_count
        self.tumor_f = np.true_divide(self.t_alt_count, self.t_depth)
        self.number_of_sites = len(self.n_alt_count)
        self.candidate_sites = np.logical_and(np.logical_and(self.tumor_f > f_thresh, self.t_depth > depth),
                                              self.n_depth > depth)
        # hyperparameter
        self.p_somatic = np.zeros([self.number_of_sites,1]) + p_somatic
        if hot_spots_file != 'NA':
            hot_spots = pd.read_csv(hot_spots_file,sep='\t',low_memory=False,index_col=False)
            if type(hot_spots['Chromosome'][0]) == str:
                hot_spots['contig'] = du.chr2num(np.array(hot_spots['Chromosome']))
            else:
                hot_spots['contig'] = np.array(hot_spots['Chromosome']) - 1
            hot_spots = hot_spots[np.isfinite(hot_spots['contig'])]
            hot_spots['genomic_coord_x'] = du.hg19_to_linear_positions(
                np.array(hot_spots['contig']), np.array(hot_spots['Position']))
            for index,hot_spot in hot_spots.iterrows():
                if np.size(np.where(self.genomic_coord_x==hot_spot['genomic_coord_x'])) > 0:
                    print 'Using user provided probabilities for cancer hot spots:'
                    print hot_spot['Chromosome'] + ' ' + hot_spot['Position']
                    self.p_somatic[np.where(self.genomic_coord_x==hot_spot['genomic_coord_x'])] = hot_spot['Probability']

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
        self.p_TiN_given_S = np.zeros([self.number_of_sites, resolution])
        self.p_TiN_given_G = np.zeros([self.number_of_sites, resolution])
        self.p_TiN_given_het = np.zeros([self.number_of_sites, resolution])
        self.p_artifact = np.zeros([self.number_of_sites, 1])

        # likelihood
        self.TiN_likelihood = np.zeros([resolution, 1])

    def generate_conditional_ps(self):
        # p(TiN|Somatic) and p(TiN|Germline)
        t_het_direction = np.ones([self.number_of_sites, len(self.af)])
        t_het_direction[:, 0:np.int(np.round(np.true_divide(len(self.af), 2)))] = -1
        self.afexp = np.repeat(np.expand_dims(self.af, 1), self.number_of_sites, axis=1).T
        t_af_w = beta._cdf(self.afexp, np.expand_dims(self.t_alt_count + 1, 1),
                           np.expand_dims(self.t_ref_count + 1, 1)) - beta._cdf(self.afexp - 0.005,
                                                                                np.expand_dims(self.t_alt_count + 1, 1),
                                                                                np.expand_dims(self.t_ref_count + 1, 1))
        f_t_af = 0.5 - np.abs(0.5 - self.afexp)
        t_af = np.multiply(self.afexp, np.expand_dims(self.n_depth, 1))

        psi_t_af = 0.5 - f_t_af
        psi_t_af = np.multiply(psi_t_af, t_het_direction)
        for TiN_idx, TiN in enumerate(self.TiN_range):
            n_ac_given_tin = np.multiply(t_af, np.expand_dims(self.CN_ratio[:, TiN_idx], 1))
            exp_f = 0.5 + np.multiply(psi_t_af, np.expand_dims(self.CN_ratio[:, TiN_idx], 1))
            n_het_ac_given_tin = np.multiply(exp_f, self.n_depth[:, np.newaxis])

            self.p_TiN_given_S[:, TiN_idx] += np.sum(
                np.multiply(beta._cdf(np.expand_dims(self.normal_f[:] + .01, 1), n_ac_given_tin + 1,
                                      self.n_depth[:, np.newaxis] - n_ac_given_tin + 1) -
                            beta._cdf(np.expand_dims(self.normal_f[:], 1), n_ac_given_tin + 1,
                                      self.n_depth[:, np.newaxis] - n_ac_given_tin + 1), t_af_w), axis=1)
            self.p_TiN_given_het[:, TiN_idx] += np.sum(
                np.multiply(beta._cdf(np.expand_dims(self.normal_f[:] + .01, 1), n_het_ac_given_tin + 1,
                                      self.n_depth[:, np.newaxis] - n_het_ac_given_tin + 1) -
                            beta._cdf(np.expand_dims(self.normal_f[:], 1), n_het_ac_given_tin + 1,
                                      self.n_depth[:, np.newaxis] - n_het_ac_given_tin + 1), t_af_w), axis=1)
        self.p_artifact = beta._cdf(self.normal_f + .01, self.t_alt_count + 1, self.t_ref_count + 1) - beta._cdf(
            self.normal_f, self.t_alt_count + 1, self.t_ref_count + 1)
        self.p_TiN_given_G = np.multiply(1 - self.p_artifact[:, np.newaxis], self.p_TiN_given_het) + np.multiply(
            self.p_artifact[:, np.newaxis], 1 - self.p_TiN_given_het)

    def expectation_of_z_given_TiN(self):
        # E step
        numerator = self.p_somatic * np.expand_dims(self.p_TiN_given_S[:, self.TiN],1)
        denominator = numerator + (1 - self.p_somatic) * np.expand_dims(np.nan_to_num(self.p_TiN_given_G[:, self.TiN]),1)
        self.E_z = np.nan_to_num(np.true_divide(numerator, denominator))
    def maximize_TiN_likelihood(self):
        # M step
        self.TiN_likelihood = np.nansum(np.multiply(self.E_z[self.candidate_sites],
                                                    np.ma.log(self.p_TiN_given_S[self.candidate_sites, :])), axis=0) + \
                              np.nansum(np.multiply(1 - self.E_z[self.candidate_sites],
                                                    np.ma.log(self.p_TiN_given_G[self.candidate_sites, :])), axis=0)
        self.TiN = np.argmax(self.TiN_likelihood)

    def perform_inference(self):
        # perform EM procedure over
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
        print 'SSNV based TiN estimate converged: TiN = ' + str(self.TiN_range[self.TiN]) + ' based on ' + str(np.sum(self.candidate_sites)) + ' sites'
        self.TiN = self.TiN_range[self.TiN]

        posterior = np.exp(self.TiN_likelihood - np.nanmax(self.TiN_likelihood))
        self.CI_tin_low = self.TiN_range[
            next(x[0] for x in
                 enumerate(np.cumsum(np.ma.masked_array((np.true_divide(posterior, np.nansum(posterior))))))
                 if x[1] > 0.025)]
        self.CI_tin_high = self.TiN_range[
            next(x[0] for x in
                 enumerate(np.cumsum(np.ma.masked_array((np.true_divide(posterior, np.nansum(posterior))))))
                 if x[1] > 0.975)]
