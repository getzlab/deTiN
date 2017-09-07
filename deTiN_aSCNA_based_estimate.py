import numpy as np
from scipy.stats import beta
from scipy.stats import norm
from scipy.cluster.vq import kmeans2

np.seterr(all='ignore')

class model:
    """Model of tumor in normal (TiN) based on only aSCNAs. This estimate is most
       only usable when tumors have sufficient allele imbalance (>200 probes and >10 SNPs under loh).
        TiN estimate : model.TiN"""
    def __init__(self,aSCNA_segs,aSCNA_hets):

        self.segs = aSCNA_segs
        self.hets = aSCNA_hets
        self.TiN_range = np.linspace(0, 1, num=101)
        self.af = np.linspace(0, 1, num=101)
        self.mu_af_n = np.mean(self.hets['AF_N'])
        self.tau = self.hets['tau']
        self.tin_correct_tau = np.multiply(self.TiN_range, self.hets['tau'][:, np.newaxis])
        self.tin_correct_normal_tau = np.multiply((1 - self.TiN_range), 2)
        self.CN_ratio = np.divide(self.tin_correct_tau, np.array(self.tin_correct_tau + self.tin_correct_normal_tau))
        self.p_TiN = np.zeros([len(self.hets),len(self.TiN_range)])
        self.seg_likelihood = dict()
        self.TiN_likelihood = np.zeros([len(self.segs), 101])

    def calculate_TiN_likelihood(self):
        t_af_w = np.zeros([len(self.hets), 101])
        rv_tumor_af = beta(self.hets['ALT_COUNT_T']+1,self.hets['REF_COUNT_T']+1)
        rv_normal_af = beta(self.hets['ALT_COUNT_N'] + 1, self.hets['REF_COUNT_N'] + 1)

        for i, f in enumerate(self.af):
            t_af_w[:, i] = rv_tumor_af.cdf(f + .01) - rv_tumor_af.cdf(f)
            f_t_af = self.mu_af_n- np.abs((self.mu_af_n - f))
            psi_t_af = self.mu_af_n - f_t_af
            psi_t_af = np.multiply(psi_t_af, self.hets['d'])
            exp_f = self.mu_af_n + np.multiply(psi_t_af[:, np.newaxis], self.CN_ratio)
            exp_f[exp_f < 0] = 0
            for TiN_idx, TiN in enumerate(self.TiN_range):
                self.p_TiN[:, TiN_idx] += np.multiply(rv_normal_af.logpdf(exp_f[:, TiN_idx]) * 0.01, t_af_w[:, i])

        seg_var = np.zeros([len(self.segs),0])
        TiN_MAP = np.zeros([len(self.segs),0])
        TiN_likelihood = np.zeros([len(self.segs),101])
        counter = 0
        for seg_id,seg in self.segs.iterrows():
            self.seg_likelihood[seg_id] = np.sum(self.p_TiN[self.hets['seg_id'] == seg_id],axis=0)
            seg_var[counter] = np.var(np.argmax(self.p_TiN[self.hets['seg_id'] == seg_id],axis=0))
            TiN_MAP[counter] = np.argmax(self.seg_likelihood[seg_id])
            TiN_likelihood[counter,:] = np.sum(self.p_TiN[self.hets['seg_id']==seg_id],axis=0)
            counter += 1
        self.segs['TiN_var'] = seg_var
        self.segs['TiN_MAP'] = TiN_MAP
        self.TiN_likelihood = TiN_likelihood

    def cluster_segments(self):
        MAP_estimate_one_cluster = np.argmax(np.sum(self.TiN_likelihood,axis=0))
        mu_c1 = np.mean(self.segs['TiN_MAP'])
        var_c1 = np.true_divide(np.sum(self.segs['TiN_var']),len(self.segs)-1)
        ll_one_cluster = np.sum(norm.logpdf(np.argmax(self.TiN_likelihood,axis=0),mu_c1,var_c1))

        kmeans2(self.segs['TiN_MAP'],iter=1000,k=2)
        mu_c21 = np.mean(self.segs['TiN_MAP'])
        mu_c22 = np.mean(self.segs['TiN_MAP'])
        var_c1 = np.true_divide(np.sum(self.segs['TiN_var']), len(self.segs) - 1)
        ll_one_cluster = np.sum(norm.logpdf(np.argmax(self.TiN_likelihood, axis=0), mu_c1, var_c1))

        return

    def perfom_inference(self):
       # self.calculate_TiN_likelihood()
       # self.cluster_segments()
        return
