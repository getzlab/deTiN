import argparse
import os
import numpy as np
import pandas as pd
from itertools import compress

import deTiN_utilities as du
import deTiN_SSNV_based_estimate as dssnv
import deTiN_aSCNA_based_estimate as dascna


class input:
    """class which holds the required detin somatic data prior to model"""

    def __init__(self, args):

        # related to inputs from command line
        self.call_stats_file = args.mutation_data_path
        self.seg_file = args.cn_data_path
        self.tumor_het_file = args.tumor_het_data_path
        self.normal_het_file = args.normal_het_data_path
        self.exac_db_file = args.exac_data_path
        self.mutation_prior = args.mutation_prior
        self.TiN_prior = args.TiN_prior
        self.output_path = args.output_dir
        self.output_name = args.output_name

        # related to inputs from class functions
        self.call_stats_table = []
        self.seg_table = []
        self.het_table = []
        self.candidates = []

    def read_call_stats_file(self):
        try:
            self.call_stats_table = pd.read_csv(self.call_stats_file, '\t', index_col=False, low_memory=False,comment='#')
        except (ValueError, LookupError):
            print 'Error reading call stats skipping first two rows and trying again'
            self.call_stats_table = pd.read_csv(self.call_stats_file, '\t', index_col=False, low_memory=False,
                                                comment='#',skiprows=2)
        if type(self.call_stats_table['contig'][0]) == str:
            self.call_stats_table['Chromosome'] = du.chr2num(np.array(self.call_stats_table['contig']))
        else:
            self.call_stats_table['Chromosome'] = np.array(self.call_stats_table['contig'])-1
        self.call_stats_table = self.call_stats_table[np.isfinite(self.call_stats_table['Chromosome'])]
        self.call_stats_table['genomic_coord_x'] = du.hg19_to_linear_positions(
            np.array(self.call_stats_table['Chromosome']), np.array(self.call_stats_table['position']))
        self.n_calls_in = len(self.call_stats_table)

    def read_het_file(self):
        tumor_het_table = pd.read_csv(self.tumor_het_file, '\t', index_col=False, low_memory=False,comment='#')
        normal_het_table = pd.read_csv(self.normal_het_file, '\t', index_col=False, low_memory=False,comment='#')
        tumor_het_table = du.fix_het_file_header(tumor_het_table)
        normal_het_table = du.fix_het_file_header(normal_het_table)
        if type(tumor_het_table['CONTIG'][0]) == str:
            tumor_het_table['Chromosome'] = du.chr2num(np.array(tumor_het_table['CONTIG']))
        else:
            tumor_het_table['Chromosome'] = np.array(tumor_het_table['CONTIG'])

        if type(normal_het_table['CONTIG'][0]) == str:
            normal_het_table['Chromosome'] = du.chr2num(np.array(normal_het_table['CONTIG']))
        else:
            normal_het_table['Chromosome'] = np.array(normal_het_table['CONTIG'])
        tumor_het_table = tumor_het_table[np.isfinite(tumor_het_table['Chromosome'])]
        tumor_het_table['genomic_coord_x'] = du.hg19_to_linear_positions(np.array(tumor_het_table['Chromosome']),
                                                                         np.array(tumor_het_table['POSITION']))
        normal_het_table = normal_het_table[np.isfinite(normal_het_table['Chromosome'])]
        normal_het_table['genomic_coord_x'] = du.hg19_to_linear_positions(np.array(normal_het_table['Chromosome']),
                                                                          np.array(normal_het_table['POSITION']))
        tumor_het_table['AF'] = np.true_divide(tumor_het_table['ALT_COUNT'],
                                               tumor_het_table['ALT_COUNT']+tumor_het_table['REF_COUNT'])
        normal_het_table['AF'] = np.true_divide(normal_het_table['ALT_COUNT'],
                                                normal_het_table['ALT_COUNT']+normal_het_table['REF_COUNT'])
        self.het_table = pd.merge(normal_het_table, tumor_het_table, on='genomic_coord_x', suffixes=('_N', '_T'))

    def read_seg_file(self):
        self.seg_table = pd.read_csv(self.seg_file, '\t', index_col=False, low_memory=False,comment='#')
        self.seg_table = du.fix_seg_file_header(self.seg_table)
        if not du.is_number(self.seg_table['Chromosome'][0]):
            self.seg_table['Chromosome'] = du.chr2num(np.array(self.seg_table['Chromosome']))
        else:
            self.seg_table['Chromosome'] = self.seg_table['Chromosome'] - 1
        self.seg_table['genomic_coord_start'] = du.hg19_to_linear_positions(np.array(self.seg_table['Chromosome']),
                                                                            np.array(self.seg_table['Start.bp']))
        self.seg_table['genomic_coord_end'] = du.hg19_to_linear_positions(np.array(self.seg_table['Chromosome']),
                                                                          np.array(self.seg_table['End.bp']))

    def annotate_call_stats_with_allelic_cn_data(self):
        f_acs = np.zeros([self.n_calls_in, 1]) + 0.5
        tau = np.zeros([self.n_calls_in, 1]) + 2
        for i, r in self.seg_table.iterrows():
            f_acs[np.logical_and(np.array(self.call_stats_table['genomic_coord_x']) >= r['genomic_coord_start'],
                                 np.array(self.call_stats_table['genomic_coord_x']) <= r['genomic_coord_end'])] = r.f
            tau[np.logical_and(np.array(self.call_stats_table['genomic_coord_x']) >= r['genomic_coord_start'],
                               np.array(self.call_stats_table['genomic_coord_x']) <= r[
                                   'genomic_coord_end'])] = r.tau + 0.001
        self.call_stats_table['tau'] = tau
        self.call_stats_table['f_acs'] = f_acs

    def annotate_het_table(self):
        seg_id = np.zeros([len(self.het_table), 1]) - 1
        tau = np.zeros([len(self.het_table),1]) + 2
        f = np.zeros([len(self.het_table),1]) + 0.5
        for seg_index, seg in self.seg_table.iterrows():
            het_index = np.logical_and(self.het_table['genomic_coord_x'] >= seg['genomic_coord_start'],
                                        self.het_table['genomic_coord_x'] <= seg['genomic_coord_end'])
            ix = list(compress(xrange(len(het_index)), het_index))
            seg_id[ix] = seg_index
            tau[ix] = seg['tau']
            f[ix] = seg['f']
        self.het_table['seg_id'] = seg_id
        self.het_table['tau'] = tau
        self.het_table['f'] = f
        d = np.ones([len(self.het_table), 1])
        d[np.array(self.het_table['AF_T'] <= 0.5,dtype=bool)] = -1
        self.het_table['d'] = d

    def read_and_preprocess_data(self):
        self.read_call_stats_file()
        self.read_seg_file()
        self.annotate_call_stats_with_allelic_cn_data()
        self.read_het_file()
        self.seg_table = du.filter_segments_based_on_size_f_and_tau(self.seg_table)
        self.annotate_het_table()
        self.het_table = du.remove_sites_near_centromere_and_telomeres(self.het_table)


class output:
    """ combined from deTiN's models
    reclassified SSNVs based on TiN estimate are labeled KEEP in judgement column
    self.SSNVs['judgement'] ==  KEEP

    confidence intervals (CI_tin_high/low) represent 95% interval
    """

    def __init__(self,input,ssnv_based_model,ascna_based_model):

        # previous results
        self.input = input
        self.ssnv_based_model = ssnv_based_model
        self.ascna_based_model = ascna_based_model

        # useful outputs
        self.SSNVs = input.candidates
        self.joint_log_likelihood = np.zeros([101, 1])
        self.joint_posterior = np.zeros([101, 1])
        self.CI_tin_high = []
        self.CI_tin_low = []
        self.TiN = []

        # variables
        self.TiN_range = np.linspace(0, 1, num=101)

    def calculate_joint_estimate(self):
        if self.ssnv_based_model.TiN <= 0.3:
            self.joint_log_likelihood = self.ascna_based_model.TiN_likelihood + self.ssnv_based_model.TiN_likelihood

            self.joint_posterior = np.exp(self.ascna_based_model.TiN_likelihood + self.ssnv_based_model.TiN_likelihood
                                      - np.nanmax(self.ascna_based_model.TiN_likelihood + self.ssnv_based_model.TiN_likelihood))
            self.joint_posterior = np.true_divide(self.joint_posterior,np.nansum(self.joint_posterior))
            self.CI_tin_low = self.TiN_range[next(x[0] for x in enumerate(np.cumsum(np.ma.masked_array(np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))))) if
                x[1] > 0.025)]
            self.CI_tin_high = self.TiN_range[
                next(x[0] for x in enumerate(np.cumsum(np.ma.masked_array(np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))))) if
                x[1] > 0.975)]

            self.TiN = self.TiN_range[np.nanargmax(self.joint_posterior)]
        else:
            print 'SSNV based TiN estimate exceed 0.3 using only aSCNA based estimate'
            self.joint_log_likelihood = self.ascna_based_model.TiN_likelihood
            self.joint_posterior = np.exp(self.ascna_based_model.TiN_likelihood - np.nanmax(self.ascna_based_model.TiN_likelihood))
            self.CI_tin_low = self.TiN_range[next(x[0] for x in enumerate(
                np.cumsum(np.ma.masked_array(np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))))) if
                                                  x[1] > 0.025)]
            self.CI_tin_high = self.TiN_range[
                next(x[0] for x in enumerate(np.cumsum(
                    np.ma.masked_array(np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))))) if
                     x[1] > 0.975)]
            self.TiN = self.TiN_range[np.nanargmax(self.joint_posterior)]

    def reclassify_mutations(self):
        # calculate E_z given joint TiN point estimate
        numerator = self.ssnv_based_model.p_somatic * self.ssnv_based_model.p_TiN_given_S[:, np.nanargmax(self.joint_posterior)]
        denominator = numerator + np.array([1 - self.ssnv_based_model.p_somatic] * self.ssnv_based_model.p_TiN_given_G[:, np.nanargmax(self.joint_posterior)])
        self.SSNVs.loc[:, ('p_somatic_given_TiN')] = np.nan_to_num(np.true_divide(numerator, denominator))

        # remove outliers mutations p(af_n >= E[af_n|TiN]) < 0.05
        af_n_given_TiN = np.multiply(self.ssnv_based_model.tumor_f,self.ssnv_based_model.CN_ratio[:,np.nanargmax(self.joint_posterior)])
        self.SSNVs.loc[:,'p_outlier'] = self.ssnv_based_model.rv_normal_af.cdf(af_n_given_TiN)
        self.SSNVs['judgement'][np.logical_and(self.SSNVs['p_somatic_given_TiN']>0.5,self.SSNVs['p_outlier']>=0.01)] = 'KEEP'


__version__ = '1.0'


def main():
    """ deTiN pipeline. Method operates in two stages (1) estimating tumor in normal via candidate SSNVs and SCNAS.
        (2) Performing variant re-classification using bayes rule.
    """

    parser = argparse.ArgumentParser(description='Estimate tumor in normal (TiN) using putative somatic'
                                                 ' events see Taylor-Weiner & Stewart et al. 2017')
    parser.add_argument('--mutation_data_path',
                        help='Path to mutation data.'
                             'Supported formats: MuTect call-stats, Strelka VCF', required=True)
    parser.add_argument('--cn_data_path',
                        help='Path to copy number data.'
                             'Supported format: GATK4CNV .seg file or AllelicCapseg .seg file', required=True)
    parser.add_argument('--tumor_het_data_path',
                        help='Path to heterozygous site allele count data in tumor.'
                             'Required columns: CONTIG,POS,REF_COUNT and ALT_COUNT', required=True)
    parser.add_argument('--normal_het_data_path',
                        help='Path to heterozygous site allele count data in normal.'
                             'Required columns: CONTIG,POS,REF_COUNT and ALT_COUNT', required=True)
    parser.add_argument('--exac_data_path',
                        help='Path to exac vcf file or pickle.', required=False)
    parser.add_argument('--output_name', required=False,
                        help='sample name')
    parser.add_argument('--mutation_prior', help='prior expected ratio of somatic mutations to rare germline events'
                        , required=False, default=0.08)
    parser.add_argument('--mutation_data_source',
                        help='Variant caller used.'
                             'Supported values: mutect,strelka,varscan,and somaticsniper', required=False,
                        default='MuTect')

    parser.add_argument('--TiN_prior', help='expected frequency of TiN contamination in sequencing setting',
                        required=False, default=0.1)
    parser.add_argument('--output_dir', help='directory to put plots and TiN solution', required=False, default='.')

    args = parser.parse_args()
    di = input(args)
    di.read_and_preprocess_data()

    # identify candidate mutations based on MuTect flags.
    # kept sites are flagged as KEEP or rejected for normal lod and/or alt_allele_in_normal
    di.candidates = du.select_candidate_mutations(di.call_stats_table)

    # generate SSNV based model using candidate sites
    ssnv_based_model = dssnv.model(di.candidates, di.mutation_prior)
    ssnv_based_model.perform_inference()

    # identify aSCNAs and filter hets
    di.aSCNA_hets = du.ensure_balanced_hets(di.seg_table,di.het_table)
    di.aSCNA_segs = du.identify_aSCNAs(di.seg_table,di.aSCNA_hets)

    # generate aSCNA based model
    ascna_based_model = dascna.model(di.aSCNA_segs, di.aSCNA_hets)
    ascna_based_model.perform_inference()

    # combine models and reclassify mutations
    do = output(di,ssnv_based_model,ascna_based_model)
    do.calculate_joint_estimate()
    do.reclassify_mutations()
    do.SSNVs.drop('Chromosome', axis=1, inplace=True)

    # make output directory if needed
    if args.output_dir != '.':
        os.makedirs(args.output_dir, exist_ok=True)
    # write deTiN reclassified SSNVs
    do.SSNVs.to_csv(path_or_buf=do.input.output_path + '/' + do.input.output_name + '.deTiN_SSNVs.txt', sep='\t',index=None)
    # write plots
    du.plot_kmeans_info(ascna_based_model, do.input.output_path, do.input.output_name)
    du.plot_TiN_models(do)
    du.plot_SSNVs(do)
    # write TiN and CIs
    file = open(do.input.output_path +'/' + do.input.output_name + '.TiN_estimate.txt', 'w')
    file.write('%s' % (do.TiN))
    file.close()

    file = open(do.input.output_path + '/' + do.input.output_name + '.TiN_estimate_CI.txt', 'w')
    file.write('%s - %s' % (str(do.CI_tin_low), str(do.CI_tin_high)))
    file.close()

if __name__ == "__main__":
    main()
