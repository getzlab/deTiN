from __future__ import division
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
        self.segs = aSCNA_segs.copy()
        self.hets = aSCNA_hets.copy()
        self.n_segs = self.segs.shape[0]
        self.n_hets = self.hets.shape[0]
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
        self.p_TiN = np.zeros([self.n_hets, len(self.TiN_range)])
        self.seg_likelihood = dict()
        self.TiN_likelihood_matrix = np.zeros([self.n_segs, resolution])
        self.reporting_cluster = 'mode'

        # model outputs
        self.TiN_likelihood = np.zeros([len(self.TiN_range), 1])
        self.TiN = 0
        self.sum_squared_distance = np.zeros([3, 1])
        self.cluster_assignment = np.zeros([self.n_segs, 1])
        self.centroids = np.zeros([3, 1])
        self.bic = np.zeros([3, 1])

    def calculate_TiN_likelihood(self):
        self.t_alt_count = self.hets.as_matrix(['ALT_COUNT_T'])
        self.t_ref_count = self.hets.as_matrix(['REF_COUNT_T'])
        self.afexp = np.repeat(np.expand_dims(self.af, 1), self.n_hets, axis=1).T
        t_af_w = beta._cdf(self.afexp, self.t_alt_count + 1, self.t_ref_count + 1) - beta._cdf(self.afexp-0.005, self.t_alt_count + 1, self.t_ref_count + 1)
        f_t_af = self.mu_af_n - np.abs(self.mu_af_n - self.afexp)
        psi_t_af = self.mu_af_n - f_t_af
        psi_t_af = np.multiply(psi_t_af, np.expand_dims(self.hets['d'], 1))
        self.n_alt_count = np.squeeze(self.hets.as_matrix(['ALT_COUNT_N']))
        self.n_ref_count = np.squeeze(self.hets.as_matrix(['REF_COUNT_N']))
        self.p_TiN = np.zeros([self.n_hets, self.resolution])
        for i, f in enumerate(self.af):
            exp_f = self.mu_af_n + np.multiply(np.expand_dims(psi_t_af[:, i], 1), self.CN_ratio)
            exp_f[exp_f < 0] = 0
            self.p_TiN += np.multiply(beta._pdf(exp_f, np.expand_dims(self.n_alt_count + 1, 1),
                                                np.expand_dims(self.n_ref_count + 1, 1)) * 0.01,
                                      np.expand_dims(t_af_w[:, i], 1))
        seg_var = np.zeros([self.n_segs, 1])
        TiN_MAP = np.zeros([self.n_segs, 1], dtype=int)
        TiN_ci_h = np.zeros([self.n_segs, 1], dtype=int)
        TiN_ci_l = np.zeros([self.n_segs, 1], dtype=int)
        TiN_likelihood = np.zeros([self.n_segs, self.resolution])
        TiN_post = np.zeros([self.n_segs, self.resolution])
        counter = 0
        for seg_id, seg in self.segs.iterrows():
            self.seg_likelihood[seg_id] = np.sum(
                np.log(self.p_TiN[np.array(self.hets['seg_id'] == seg_id, dtype=bool)]), axis=0)
            seg_var[counter] = np.nanvar(
                np.argmax(self.p_TiN[np.array(self.hets['seg_id'] == seg_id, dtype=bool)], axis=0))
            TiN_MAP[counter] = np.nanargmax(self.seg_likelihood[seg_id])
            TiN_likelihood[counter, :] = np.nansum(np.log(self.p_TiN[np.array(self.hets['seg_id'] == seg_id, dtype=bool)]),
                                                axis=0)
            prior = np.true_divide(np.ones([1, self.resolution]), self.resolution)
            TiN_post[counter, :] = TiN_likelihood[counter, :] + np.log(prior)
            TiN_post[counter, :] = TiN_post[counter, :] + (1 - np.max(TiN_post[counter, :]))
            TiN_post[counter, :] = np.exp(TiN_post[counter, :])
            TiN_post[counter, :] = np.true_divide(TiN_post[counter, :], np.nansum(TiN_post[counter, :]))
            TiN_ci_l[counter, :] = self.TiN_range[map(lambda x: x>0.025, np.cumsum(TiN_post[counter, :])).index(True)]*100
            TiN_ci_h[counter, :] = self.TiN_range[map(lambda x: x>0.975, np.cumsum(TiN_post[counter, :])).index(True)]*100
            counter += 1
        self.TiN_post_seg = TiN_post
        self.segs.loc[:, ('TiN_var')] = seg_var
        self.segs.loc[:, ('TiN_MAP')] = self.TiN_range[TiN_MAP] * 100
        self.TiN_likelihood_matrix = TiN_likelihood
        self.segs.loc[:,('TiN_ci_h')] = TiN_ci_h
        self.segs.loc[:, ('TiN_ci_l')] = TiN_ci_l
    def cluster_segments(self):
        if self.n_segs >= 3:
            K = range(1, 4)
            N = len(self.hets['seg_id'])
            self.segs.reset_index(inplace=True, drop=False)
            tin_data = np.nanargmax(self.TiN_likelihood_matrix,axis=1).astype(float)
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
            if len(centroids[2])>2:
                dist_btwn_c3 = np.min([abs(i - j) for i, j in combinations(centroids[2], 2)])
            else:
                dist_btwn_c3 = 0
            if len(centroids[1]) > 1:
                dist_btwn_c2 = np.abs(np.diff(centroids[1]))
            else:
                dist_btwn_c2 = 0
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
            np.unique(self.segs['Chromosome']) +1)
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
