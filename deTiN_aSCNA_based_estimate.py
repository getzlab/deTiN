import numpy as np
from scipy.stats import beta
from scipy.stats import norm
from scipy.cluster.vq import kmeans
from scipy.stats import mode
from itertools import combinations
import pandas as pd

np.seterr(all='ignore')
pd.options.mode.chained_assignment = None  # default='warn'

class model:
    """Model of tumor in normal (TiN) based on only aSCNAs. This estimate is most
       only usable when tumors have sufficient allele imbalance (>200 probes and >10 SNPs under loh).
        TiN estimate : model.TiN"""

    def __init__(self, aSCNA_segs, aSCNA_hets):

        # input data
        self.segs = aSCNA_segs
        self.hets = aSCNA_hets

        # Variables for fit
        self.TiN_range = np.linspace(0, 1, num=101)
        self.af = np.linspace(0, 1, num=101)

        # model parameters
        self.mu_af_n = np.mean(self.hets['AF_N'])
        self.tau = self.hets['tau']
        self.tin_correct_tau = np.multiply(self.TiN_range, self.hets['tau'][:, np.newaxis])
        self.tin_correct_normal_tau = np.multiply((1 - self.TiN_range), 2)
        self.CN_ratio = np.divide(self.tin_correct_tau, np.array(self.tin_correct_tau + self.tin_correct_normal_tau))
        self.p_TiN = np.zeros([len(self.hets), len(self.TiN_range)])
        self.seg_likelihood = dict()
        self.TiN_likelihood_matrix = np.zeros([len(self.segs), 101])
        self.reporting_cluster = 'mode'

        # model outputs
        self.TiN_likelihood = np.zeros([101, 1])
        self.TiN = 0
        self.mean_sum_squared_distance = np.zeros([3,1])
        self.cluster_assignment = np.zeros([len(self.segs['Chromosome']),1])
        self.centroids = np.zeros([3, 1])

    def calculate_TiN_likelihood(self):
        t_af_w = np.zeros([len(self.hets), 101])
        rv_tumor_af = beta(self.hets['ALT_COUNT_T'] + 1, self.hets['REF_COUNT_T'] + 1)
        rv_normal_af = beta(self.hets['ALT_COUNT_N'] + 1, self.hets['REF_COUNT_N'] + 1)

        for i, f in enumerate(self.af):
            t_af_w[:, i] = rv_tumor_af.cdf(f + .01) - rv_tumor_af.cdf(f)
            f_t_af = self.mu_af_n - np.abs((self.mu_af_n - f))
            psi_t_af = self.mu_af_n - f_t_af
            psi_t_af = np.multiply(psi_t_af, self.hets['d'])
            exp_f = self.mu_af_n + np.multiply(psi_t_af[:, np.newaxis], self.CN_ratio)
            exp_f[exp_f < 0] = 0
            for TiN_idx, TiN in enumerate(self.TiN_range):
                self.p_TiN[:, TiN_idx] += np.multiply(rv_normal_af.logpdf(exp_f[:, TiN_idx]) * 0.01, t_af_w[:, i])

        seg_var = np.zeros([len(self.segs), 1])
        TiN_MAP = np.zeros([len(self.segs), 1])
        TiN_likelihood = np.zeros([len(self.segs), 101])
        counter = 0
        for seg_id, seg in self.segs.iterrows():
            self.seg_likelihood[seg_id] = np.sum(self.p_TiN[self.hets['seg_id'] == seg_id], axis=0)
            seg_var[counter] = np.nanvar(np.argmax(self.p_TiN[self.hets['seg_id'] == seg_id], axis=0))
            TiN_MAP[counter] = np.nanargmax(self.seg_likelihood[seg_id])
            TiN_likelihood[counter, :] = np.sum(self.p_TiN[self.hets['seg_id'] == seg_id], axis=0)
            counter += 1
        self.segs.loc[:, ('TiN_var')] = seg_var
        self.segs.loc[:, ('TiN_MAP')] = TiN_MAP
        self.TiN_likelihood_matrix = TiN_likelihood

    def cluster_segments(self):
        K = range(1, 4)
        km = [kmeans(self.segs['TiN_MAP'], k, iter=1000) for k in K]
        N = len(self.segs['TiN_MAP'])
        centroids = [cent for (cent, var) in km]
        squared_distance_to_centroids = [np.power(np.subtract(self.segs['TiN_MAP'][:, np.newaxis], cent), 2) for cent in
                                         centroids]
        self.mean_sum_squared_distance = [sum(np.min(d,axis=1))/N for d in squared_distance_to_centroids]
        cluster_assignment = [np.argmin(d, axis=1) for d in squared_distance_to_centroids]

        cl_var = np.zeros([3, 1])
        ll_cluster = np.zeros([3, 1])
        for m in range(3):
            cl_var[m] = (1.0 / (N - m)) * sum([sum(
                np.power(np.subtract(self.segs['TiN_MAP'][cluster_assignment[m] == i], centroids[m][i]), 2))
                for i in range(m + 1)])
            ll_cluster[m] = sum(
                [np.sum(norm.logpdf(self.segs['TiN_MAP'][cluster_assignment[m] == i], centroids[m][i], cl_var[m]))
                 for i in range(m + 1)])
            p = [2, 6, 12]
            bic = np.multiply(2, ll_cluster).T + np.multiply(p, np.log(N))
            dist_btwn_c3 = np.min([abs(i - j) for i, j in combinations(centroids[2], 2)])
            dist_btwn_c2 = np.abs(np.diff(centroids[1]))
            if dist_btwn_c3 < 2 * np.sqrt(cl_var[2]) and dist_btwn_c2 > 2 * np.sqrt(cl_var[1]):
                solution_idx = np.argmax(bic[0:1])
                self.cluster_assignment = cluster_assignment[solution_idx]
                self.centroids = centroids[solution_idx]
            if dist_btwn_c3 < 2 * np.sqrt(cl_var[2]) and dist_btwn_c2 < 2 * np.sqrt(cl_var[1]):
                self.cluster_assignment = cluster_assignment[0]
                self.centroids = centroids[0]
            else:
                solution_idx = np.argmax(bic)
                self.cluster_assignment = cluster_assignment[solution_idx]
                self.centroids = centroids[solution_idx]

    def perform_inference(self):
        print 'calculating aSCNA based TiN estimate using data from chromosomes: ' + str(np.unique(self.segs['Chromosome']))
        self.calculate_TiN_likelihood()
        self.cluster_segments()
        if np.max(self.cluster_assignment) > 0:
            print 'detected ' + str(np.max(self.cluster_assignment)+1) +' clusters'
            if self.reporting_cluster == 'mode':
                mode_cluster = mode(self.cluster_assignment)[0][0]
                self.TiN = self.TiN_range[np.nanargmax(np.sum(self.TiN_likelihood_matrix[self.cluster_assignment==mode_cluster,:],axis=0))]
                self.TiN_likelihood = np.sum(self.TiN_likelihood_matrix[self.cluster_assignment==mode_cluster,:],axis=0)
                print 'aSCNA based TiN estimate from modal TiN cluster :  ' +\
                      str(self.TiN)
        else:
            self.TiN = self.TiN_range[np.nanargmax(np.sum(self.TiN_likelihood_matrix,axis=0))]
            self.TiN_likelihood = np.sum(self.TiN_likelihood_matrix,axis=0)
            print 'aSCNA based TiN estimate: TiN =  ' + str(self.TiN)
