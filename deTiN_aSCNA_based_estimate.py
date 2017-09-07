import numpy as np
from scipy.stats import beta

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

    def calculate_map_estimate(self):

        for seg_id,seg in self.segs.iterrows():
            seg_hets = self.hets[self.hets['seg_id'] == seg_id]
            t_af_w = np.zeros([len(seg_hets), 101])
            rv_tumor_af = beta(seg_hets['ALT_COUNT_T']+1,seg_hets['REF_COUNT_T']+1)
            for i, f in enumerate(self.af):
                t_af_w[:, i] = rv_tumor_af.cdf(f + .01) - rv_tumor_af.cdf(f)
                f_t_af = 0.5 - np.abs((.5 - f))
                psi_t_af = 0.5 - f_t_af
                psi_t_af = np.multiply(psi_t_af, t_het_direction[:, i])
                exp_f = 0.5 + np.multiply(psi_t_af[:, np.newaxis], self.CN_ratio)