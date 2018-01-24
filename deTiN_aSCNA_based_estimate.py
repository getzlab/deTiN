import numpy as np
from scipy.stats import beta
from scipy.stats import norm
from sklearn import preprocessing
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

    def __init__(self, aSCNA_segs, aSCNA_hets, resolution = 101):

        # input data
        self.segs = aSCNA_segs
        self.hets = aSCNA_hets
        self.resolution = resolution
        # Variables for fit
        self.TiN_range = np.linspace(0, 1, num=resolution)
        self.af = np.linspace(0.005, 1, num=200)

        # model parameters
        self.mu_af_n = np.mean(self.hets['AF_N'])
        self.tau = self.hets['tau']
        self.tin_correct_tau = np.multiply(self.TiN_range, self.hets['tau'][:, np.newaxis])
        self.tin_correct_normal_tau = np.multiply((1 - self.TiN_range), 2)
        self.CN_ratio = np.divide(self.tin_correct_tau, np.array(self.tin_correct_tau + self.tin_correct_normal_tau))
        self.p_TiN = np.zeros([len(self.hets), len(self.TiN_range)])
        self.seg_likelihood = dict()
        self.TiN_likelihood_matrix = np.zeros([len(self.segs), resolution])
        self.reporting_cluster = 'mode'

        # model outputs
        self.TiN_likelihood = np.zeros([resolution, 1])
        self.TiN = 0
        self.sum_squared_distance = np.zeros([3, 1])
        self.cluster_assignment = np.zeros([len(self.segs['Chromosome']), 1])
        self.centroids = np.zeros([3, 1])
        self.bic = np.zeros([3, 1])

    def calculate_TiN_likelihood(self):
        t_af_w = np.zeros([len(self.hets), len(self.af)])
        rv_tumor_af = beta(self.hets['ALT_COUNT_T'] + 1, self.hets['REF_COUNT_T'] + 1)
        rv_normal_af = beta(self.hets['ALT_COUNT_N'] + 1, self.hets['REF_COUNT_N'] + 1)

        for i, f in enumerate(self.af):
            t_af_w[:, i] = rv_tumor_af.cdf(f) - rv_tumor_af.cdf(f - 0.005)
            f_t_af = self.mu_af_n - np.abs((self.mu_af_n - f))
            psi_t_af = self.mu_af_n - f_t_af
            psi_t_af = np.multiply(psi_t_af, self.hets['d'])
            exp_f = self.mu_af_n + np.multiply(psi_t_af[:, np.newaxis], self.CN_ratio)
            exp_f[exp_f < 0] = 0
            for TiN_idx, TiN in enumerate(self.TiN_range):
                self.p_TiN[:, TiN_idx] += np.multiply(rv_normal_af.pdf(exp_f[:, TiN_idx]) * 0.01, t_af_w[:, i])

        seg_var = np.zeros([len(self.segs), 1])
        TiN_MAP = np.zeros([len(self.segs), 1],dtype=int)
        TiN_likelihood = np.zeros([len(self.segs), self.resolution])
        TiN_post = np.zeros([len(self.segs), self.resolution])
        counter = 0
        for seg_id, seg in self.segs.iterrows():
            self.seg_likelihood[seg_id] = np.sum(
                np.log(self.p_TiN[np.array(self.hets['seg_id'] == seg_id, dtype=bool)]), axis=0)
            seg_var[counter] = np.nanvar(
                np.argmax(self.p_TiN[np.array(self.hets['seg_id'] == seg_id, dtype=bool)], axis=0))
            TiN_MAP[counter] = np.nanargmax(self.seg_likelihood[seg_id])
            TiN_likelihood[counter, :] = np.sum(np.log(self.p_TiN[np.array(self.hets['seg_id'] == seg_id, dtype=bool)]),
                                                axis=0)
            prior = np.true_divide(np.ones([1, self.resolution]), self.resolution)
            TiN_post[counter, :] = TiN_likelihood[counter, :] + np.log(prior)
            TiN_post[counter, :] = TiN_post[counter, :] + (1 - np.max(TiN_post[counter, :]))
            TiN_post[counter, :] = np.exp(TiN_post[counter, :])
            TiN_post[counter, :] = np.true_divide(TiN_post[counter, :], np.sum(TiN_post[counter, :]))
            counter += 1
        self.TiN_post_seg = TiN_post
        self.segs.loc[:, ('TiN_var')] = seg_var
        self.segs.loc[:, ('TiN_MAP')] = self.TiN_range[TiN_MAP]*100
        self.TiN_likelihood_matrix = TiN_likelihood

    def cluster_segments(self):
        if len(self.segs) > 1:
            K = range(1, 4)
            N = len(self.hets['seg_id'])
            self.segs.reset_index(inplace=True, drop=False)
            tin_data = np.array(self.segs['TiN_MAP'])
            km = [kmeans(tin_data, k, iter=1000) for k in K]
            centroids = [cent for (cent, var) in km]
            squared_distance_to_centroids = [np.power(np.subtract(tin_data[:, np.newaxis], cent), 2) for cent in
                                             centroids]
            self.sum_squared_distance = [sum(np.min(d, axis=1)) / N for d in squared_distance_to_centroids]
            cluster_assignment = [np.argmin(d, axis=1) for d in squared_distance_to_centroids]
            het_tin_map = np.argmax(self.p_TiN, axis=1)
            self.cl_distance_points = np.zeros([3, 3])
            for k, clust in enumerate(cluster_assignment):
                for idx, row in self.segs.iterrows():
                    self.cl_distance_points[k, clust[idx]] += np.sum(
                        np.power(het_tin_map[self.hets['seg_id'] == row['index']] - centroids[k][clust[idx]], 2))

            self.cl_var = np.sqrt(
                np.true_divide(self.cl_distance_points, len(self.hets['seg_id'])))
            p = [1, 2, 3]
            delta_bic = [0, 10, 20]
            self.bic = (
                       np.multiply(N, np.log(np.true_divide(np.sum(self.cl_distance_points, axis=1), N))) + np.multiply(
                           p,
                           np.log(
                               N))) + delta_bic
            dist_btwn_c3 = np.min([abs(i - j) for i, j in combinations(centroids[2], 2)])
            dist_btwn_c2 = np.abs(np.diff(centroids[1]))
            if dist_btwn_c3 < 2 * np.nanmax(self.cl_var[2, :]) and dist_btwn_c2 > 2 * np.nanmax(self.cl_var[1, :]):
                solution_idx = np.nanargmin(self.bic[0:1])
                self.cluster_assignment = cluster_assignment[solution_idx]
                self.centroids = centroids[solution_idx]
            if dist_btwn_c3 < 2 * np.nanmax(self.cl_var[2, :]) and dist_btwn_c2 < 2 * np.nanmax(self.cl_var[1, :]):
                self.cluster_assignment = cluster_assignment[0]
                self.centroids = centroids[0]
            else:
                solution_idx = np.nanargmin(self.bic)
                self.cluster_assignment = cluster_assignment[solution_idx]
                self.centroids = centroids[solution_idx]
        else:
            self.cluster_assignment = 0
            self.centroids = self.segs['TiN_MAP']

    def perform_inference(self):
        # MAP estimation of TiN using copy number data
        print 'calculating aSCNA based TiN estimate using data from chromosomes: ' + str(
            np.unique(self.segs['Chromosome']) + 1)
        # calculate likelihood function for TiN in each segment
        self.calculate_TiN_likelihood()
        # perform k-means clustering on TiN segment data
        self.cluster_segments()
        if np.max(self.cluster_assignment) > 0:
            print 'detected ' + str(np.max(self.cluster_assignment) + 1) + ' clusters'
            if self.reporting_cluster == 'mode':
                mode_cluster = mode(self.cluster_assignment)[0][0]
            elif self.reporting_cluster == 'min':
                mode_cluster = self.cluster_assignment[np.argmin(self.centroids)]

            self.TiN = self.TiN_range[
                np.nanargmax(np.sum(self.TiN_likelihood_matrix[self.cluster_assignment == mode_cluster, :], axis=0))]
            self.TiN_likelihood = np.sum(self.TiN_likelihood_matrix[self.cluster_assignment == mode_cluster, :], axis=0)
            posterior = np.exp(self.TiN_likelihood - np.nanmax(self.TiN_likelihood))
            self.CI_tin_low = self.TiN_range[
                next(x[0] for x in
                     enumerate(np.cumsum(np.ma.masked_array(np.true_divide(posterior, np.nansum(posterior)))))
                     if x[1] > 0.025)]
            self.CI_tin_high = self.TiN_range[
                next(x[0] for x in
                     enumerate(np.cumsum(np.ma.masked_array(np.true_divide(posterior, np.nansum(posterior)))))
                     if x[1] > 0.975)]
            print 'aSCNA based TiN estimate from selected TiN cluster :  ' + str(self.TiN)
        else:
            self.TiN = self.TiN_range[np.nanargmax(np.sum(self.TiN_likelihood_matrix, axis=0))]
            self.TiN_likelihood = np.sum(self.TiN_likelihood_matrix, axis=0)
            posterior = np.exp(self.TiN_likelihood - np.nanmax(self.TiN_likelihood))
            self.CI_tin_low = self.TiN_range[
                next(x[0] for x in
                     enumerate(np.cumsum(np.ma.masked_array(np.true_divide(posterior, np.nansum(posterior)))))
                     if x[1] > 0.025)]
            self.CI_tin_high = self.TiN_range[
                next(x[0] for x in
                     enumerate(np.cumsum(np.ma.masked_array(np.true_divide(posterior, np.nansum(posterior)))))
                     if x[1] > 0.975)]
            print 'aSCNA based TiN estimate: TiN =  ' + str(self.TiN)
