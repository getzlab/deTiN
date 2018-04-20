import numpy as np
import sys
from scipy.stats import beta
from scipy.stats import fisher_exact
from itertools import compress
import gzip
import random
import pandas as pd
import matplotlib
import pickle
matplotlib.use('agg')
import matplotlib.pyplot as plt

random.seed(1)


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def is_member(a, b):
    # based on the matlab is_member function
    # code from stackoverflow user Joe Kington
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return [bind.get(itm, np.nan) for itm in a]


def chr2num(chr):
    # convert chromosome from strings to ints
    chr[chr == 'X'] = '23'
    chr[chr == 'Y'] = '24'
    chr[np.array(chr == 'MT') | np.array(chr == 'M')] = '25'
    chromosomes = np.array(range(1, 26))
    return np.array(is_member(chr, chromosomes.astype(np.str)))


def remove_sites_near_centromere_and_telomeres(het_table):
    positions = het_table['genomic_coord_x']
    centromere_positions = [125000001, 1718373143, 1720573143, 1867507890,
                            1869607890, 1984214406, 1986714406, 2101066301,
                            2102666301, 2216036179, 2217536179, 2323085719, 2326285719, 2444417111,
                            2446417111, 2522371864, 2524171864, 2596767074, 2598567074, 2683844322,
                            2685944322, 339750621, 342550621, 2744173305, 2746073305, 2792498825, 2794798825,
                            2841928720,
                            2844428720, 580349994, 583449994, 738672424, 740872424, 927726700, 930026700, 1121241960,
                            1123541960, 1291657027, 1293557027, 1435895690, 1438395690, 1586459712, 1588159712]

    lengths = np.array([249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                        146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                        102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                        155270560, 59373566, 16569])  # chromosome lengths from genome-mysql.cse.ucsc.edu
    telomere_positions = np.append(1, np.cumsum(lengths))
    distance_centromere = np.zeros([len(positions), len(centromere_positions)])
    distance_telomere = np.zeros([len(positions), len(telomere_positions)])

    for i, centromere in enumerate(centromere_positions):
        distance_centromere[:, i] = np.abs(positions - centromere)
    distance_centromere = np.min(distance_centromere, axis=1)

    for i, telomere in enumerate(telomere_positions):
        distance_telomere[:, i] = np.abs(positions - telomere)
    distance_telomere = np.min(distance_telomere, axis=1)
    het_table = het_table[np.logical_and(distance_centromere > 5000000, distance_telomere > 5000000)]
    het_table.reset_index(inplace=True, drop=True)
    return het_table


def filter_hets_based_on_coverage(het_table, depth = 10):
    het_table = het_table[np.logical_and(het_table['READ_DEPTH_N'] > depth, het_table['READ_DEPTH_T'] > depth)]
    het_table.reset_index(inplace=True, drop=True)
    return het_table


def filter_segments_based_on_size_f_and_tau(seg_table, aSCNA_thresh, n_probes = 200):
    seg_table = seg_table[np.logical_and.reduce(np.array([np.array(seg_table['f']) < 0.5 - aSCNA_thresh,
                                                          seg_table['n_probes'] > n_probes, seg_table['tau'] > 0]))]
    seg_table.reset_index(inplace=True, drop=True)
    return seg_table

def alternate_file_headers():

    headers = {'alternate_headers_position' : ['Start', 'Start_bp', 'start','position','pos','POS','Start_position'],
               'alternate_headers_start_position' : ['Start', 'Start_bp', 'start','position','pos','POS','Start_position'],
               'alternate_headers_end_position' : ['End', 'End_bp', 'end'],
               'alternate_headers_chromosome' : ['Contig', 'chrom', 'CONTIG', 'chr', 'Chrom', 'CHROMOSOME','Chromosome','contig'],
               'alternate_headers_f' : ['f_acs', 'MAF_Post_Mode'],
               'alternate_headers_tau' : ['CN', 'Segment_Mean_Post_Mode'],
               'alternate_headers_alt_count' : ['t_alt_count', 'n_alt_count', 'alt_count', 'i_t_alt_count', 'i_n_alt_count'],
               'alternate_headers_ref_count' : ['t_ref_count', 'n_ref_count', 'ref_count', 'i_t_ref_count', 'i_n_ref_count'] }
    return headers

def read_file_header(text_file):

    headers = alternate_file_headers()
    with open(text_file, 'rt') as f:
        for header_lines, line in enumerate(f):
            line = line.strip()
            if not line[0] == '#':
                break
    file_head = line.split('\t')
    try :
        headers['alternate_headers_chromosome'].index(file_head[0])
    except ValueError:
        sys.exit('The first column of all input files should be chromosome: could not find any of the chromosome headers in the first column of '+
                 text_file)
    return file_head

def identify_aSCNAs(seg_table, het_table, aSCNA_thresh = 0.1, n_snps = 20, var_thresh = 0.025):
    # identify aSCNAs based on minor allele fraction of segments
    mu_af_n = np.mean(het_table['AF_N'])
    f_detin = np.zeros([len(seg_table), 1])
    f_variance = np.zeros([len(seg_table), 1])
    n_snps_above_mu = np.zeros([len(seg_table), 1])
    n_snps_below_mu = np.zeros([len(seg_table), 1])
    fishers_p_convergent_seg = np.ones([len(seg_table), 1])
    thresh_snps = np.round(np.true_divide(n_snps,2))
    for seg_id, seg in seg_table.iterrows():
        seg_hets = het_table[het_table['seg_id'] == seg_id]
        f_detin[seg_id] = mu_af_n - np.mean(np.abs(seg_hets['AF_T'] - mu_af_n))
        f_variance[seg_id] = np.var(np.abs(seg_hets['AF_T'] - mu_af_n))
        n_snps_above_mu[seg_id] = np.sum(seg_hets['AF_T'] > mu_af_n)
        n_snps_below_mu[seg_id] = np.sum(seg_hets['AF_T'] <= mu_af_n)
        try:
            fe_tuple = fisher_exact([[np.sum(np.logical_and(seg_hets['AF_T'] > mu_af_n,
                                                            seg_hets['AF_N'] > mu_af_n)),
                                      np.sum(np.logical_and(seg_hets['AF_T'] > mu_af_n,
                                                            seg_hets['AF_N'] <= mu_af_n))],
                                     [np.sum(np.logical_and(seg_hets['AF_T'] <= mu_af_n,
                                                            seg_hets['AF_N'] > mu_af_n)),
                                      np.sum(np.logical_and(seg_hets['AF_T'] <= mu_af_n,
                                                            seg_hets['AF_N'] <= mu_af_n))]], 'less')
            fishers_p_convergent_seg[seg_id] = fe_tuple[1]
        except ValueError:
            fishers_p_convergent_seg[seg_id] = 1
    seg_table['f_detin'] = f_detin
    seg_table['f_variance'] = f_variance
    seg_table['n_snps_above_mu'] = n_snps_above_mu
    seg_table['n_snps_below_mu'] = n_snps_below_mu
    seg_table['fishers_p_convergent_seg'] = fishers_p_convergent_seg
    if any((seg_table['fishers_p_convergent_seg'] * len(seg_table)) < 0.05):
        segs = (seg_table['fishers_p_convergent_seg'] * len(seg_table)) < 0.05
        ix = list(compress(xrange(len(segs)), segs))
        print 'identified convergent aSCNA in normal on chromosomes:' + str(np.unique(seg_table['Chromosome'][ix] + 1))
        convergent_segs = seg_table[seg_table['fishers_p_convergent_seg'] * len(seg_table) <= 0.05]
    else:
        convergent_segs = None
    aSCNAs = seg_table[
        np.logical_and.reduce(np.array([np.array(seg_table['fishers_p_convergent_seg'] * len(seg_table)) > 0.05,
                                        seg_table['n_snps_above_mu'] > thresh_snps,
                                        seg_table['n_snps_below_mu'] > thresh_snps,
                                        seg_table['f_detin'] <= 0.5 - aSCNA_thresh,
                                        seg_table['f_variance'] < var_thresh]))]
    return aSCNAs,convergent_segs


def ensure_balanced_hets(seg_table, het_table):
    seg_table['aSCNA'] = np.zeros([len(seg_table), 1])
    aSCNA_hets = []
    for seg_id, seg in seg_table.iterrows():
        seg_hets = het_table[het_table['seg_id'] == seg_id]
        if np.sum(seg_hets['d'] == -1) > 10 and np.sum(seg_hets['d'] == 1) > 10:
            if sum(seg_hets['AF_T'] > 0.5) < sum(seg_hets['AF_T'] <= 0.5):
                sites = seg_hets['AF_T'] <= 0.5
                index = list(compress(xrange(len(sites)), sites))
                ixs = random.sample(index, (sum(seg_hets['AF_T'] <= 0.5) - sum(seg_hets['AF_T'] > 0.5)))
                seg_hets = seg_hets.drop(seg_hets.index[[ixs]])
                seg_hets.reset_index(inplace=True, drop=True)

            if sum(seg_hets['AF_T'] > 0.5) > sum(seg_hets['AF_T'] <= 0.5):
                sites = seg_hets['AF_T'] > 0.5
                index = list(compress(xrange(len(sites)), sites))
                ixs = random.sample(index, (sum(seg_hets['AF_T'] > 0.5) - sum(seg_hets['AF_T'] <= 0.5)))
                seg_hets = seg_hets.drop(seg_hets.index[[ixs]])
                seg_hets.reset_index(inplace=True, drop=True)
            if len(aSCNA_hets) == 0:
                aSCNA_hets = seg_hets
            else:
                aSCNA_hets = pd.concat([aSCNA_hets, seg_hets])
                aSCNA_hets.reset_index(inplace=True, drop=True)
    return aSCNA_hets


def plot_kmeans_info(ascna_based_model, output_path, sample_name):
    # method for plotting clustering results of aSCNA TiN estimates
    X = np.array(ascna_based_model.segs['TiN_MAP'])
    X_low = np.array(ascna_based_model.segs['TiN_ci_l'])
    X_high = np.array(ascna_based_model.segs['TiN_ci_h'])
    Y = np.array(ascna_based_model.segs['Chromosome'])
    kIdx = np.max(ascna_based_model.cluster_assignment)
    K = range(1, 4)

    # variance explained by incorporating additional clusters
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(K, ascna_based_model.sum_squared_distance, 'b.-')
    ax.plot(K[kIdx], ascna_based_model.sum_squared_distance[kIdx], marker='o', markersize=12,
            markeredgewidth=2, markeredgecolor='r', markerfacecolor='None')
    plt.grid(True)
    plt.xlabel('Number of clusters')
    plt.ylabel('Average within-cluster sum of squares')
    plt.title('KMeans residual')
    plt.xticks([1,2,3])
    fig.set_dpi(150)
    fig.savefig(output_path + '/' + sample_name + '_KmeansEval_plot.png', bbox_inches='tight')

    # scatter plot of TiN estimates per segment by chromosome location and cluster
    fig = plt.figure()
    ax = fig.add_subplot(111)
    clr = ['b', 'g', 'r']
    if len(X) > 1:
        for i in range(K[kIdx]):
            ind = (ascna_based_model.cluster_assignment == i)
            ax.errorbar(X[ind], Y[ind], xerr=[X[ind]-X_low[ind],X_high[ind]-X[ind]] , c=clr[i], label='Cluster %d' % i,ls='None',marker='.')
    else:
        ax.errorbar(X,Y+1,xerr=[X-X_low,X_high-X],c='b',label='Cluster 1',ls='None',marker='.')

    plt.xlabel('MAP tumor in normal estimate (%)')
    plt.ylabel('Chromosome')
    plt.title('Cluster by chromosome and TiN')
    plt.yticks(np.arange(min(Y) , max(Y) + 2, 2.0))
    plt.xticks(np.arange(0, max(X) + 1, np.max([np.round(np.true_divide(np.max(X), 10)), 1])))
    ax.set_xlim([-2, np.max(X) + 2])

    fig.set_dpi(150)
    fig.savefig(output_path + '/' + sample_name + '_KmeansEval_scatter_plot.png', bbox_inches='tight')

def plot_aSCNA_het_data(do):
    fig, ax = plt.subplots(1, 1)
    ax.plot(do.input.het_table['genomic_coord_x'], do.input.het_table['AF_T'], c=[0.5, 0.5, 0.5], marker='.', ls='None',
            ms=1, alpha=0.5)
    tumor_af = ax.plot(do.ascna_based_model.hets['genomic_coord_x'], do.ascna_based_model.hets['AF_T'], c=[0, 0, 1], marker='.',
            ls='None', ms=5)
    normal_af = ax.plot(do.ascna_based_model.hets['genomic_coord_x'], do.ascna_based_model.hets['AF_N'], c=[1, 0, 0], marker='.',
            ls='None', ms=5)
    fig.set_dpi(300)
    chrs = hg19_to_linear_positions(np.linspace(0, 23, 24), np.ones([23]))
    for c in chrs:
        ax.plot([c, c], [0, 1], 'k--')
    plt.legend(handles=[tumor_af[0], normal_af[0]],labels=['Tumor', 'Normal'])
    ax.set_xticks((chrs[1:] + chrs[:-1]) / 2)
    ax.set_xticklabels((np.linspace(1, 24, 24, dtype=int)), size=5, rotation=90)
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.set_yticklabels(np.linspace(0, 1, 5), size=5)
    ax.set_xlabel('Chromosomes')
    ax.set_ylabel('Allele fraction')
    fig.set_dpi(150)
    fig.savefig(do.input.output_path + '/' + do.input.output_name + '_TiN_hets_aSCNA_model.png', bbox_inches='tight')
    fig.savefig(do.input.output_path + '/' + do.input.output_name + '_TiN_hets_aSCNA_model.eps', bbox_inches='tight')

def plot_TiN_models(do):
    fig, ax = plt.subplots(1, 1)
    TiN_range = np.linspace(0, 1, num=do.input.resolution)
    if ~np.isnan(do.ascna_based_model.TiN):
        ascna = ax.plot(TiN_range,
                        np.true_divide(np.exp(
                            do.ascna_based_model.TiN_likelihood - np.nanmax(do.ascna_based_model.TiN_likelihood)),
                                       np.nansum(np.exp(do.ascna_based_model.TiN_likelihood - np.nanmax(
                                           do.ascna_based_model.TiN_likelihood))))
                        , 'r--', lw=1)
    ssnv = ax.plot(TiN_range,
                   np.true_divide(
                       np.exp(do.ssnv_based_model.TiN_likelihood - np.nanmax(do.ssnv_based_model.TiN_likelihood)),
                       np.nansum(
                           np.exp(do.ssnv_based_model.TiN_likelihood - np.nanmax(do.ssnv_based_model.TiN_likelihood))))
                   , 'b--', lw=1)

    joint = ax.plot(TiN_range, do.joint_posterior
                    , 'k-', lw=2)
    plt.xlabel('Tumor in normal estimate')
    plt.ylabel('p(TiN=x)')
    plt.title('TiN estimate posterior')
    if ~np.isnan(do.ascna_based_model.TiN):
        plt.legend(handles=[ascna[0], ssnv[0], joint[0]], labels=['aSCNA', 'SSNV', 'Joint Est.'])
    else:
        plt.legend(handles=[ssnv[0], joint[0]], labels=['SSNV', 'Joint Est.'])
    fig.set_dpi(150)
    fig.savefig(do.input.output_path + '/' + do.input.output_name + '_TiN_models_plot.png', bbox_inches='tight')


def plot_SSNVs(do):
    fig, ax = plt.subplots(1, 1)
    TiN_fit = ax.plot(np.linspace(0, 1, do.input.resolution), np.multiply(do.TiN, np.linspace(0, 1, do.input.resolution)), '--', lw=1, alpha=1,
                      color='#1D1D1D')
    background = ax.plot(do.ssnv_based_model.tumor_f, do.ssnv_based_model.normal_f
                         , '.', lw=0.1, alpha=0.75, color=[0.75, 0.75, 0.75])
    nod_kept = np.logical_and(do.SSNVs['judgement'] == 'KEEP', do.SSNVs.isnull()['failure_reasons']).values
    cis = do.ssnv_based_model.rv_normal_af.interval(0.6825)

    kept_def = ax.plot(do.ssnv_based_model.tumor_f[nod_kept], do.ssnv_based_model.normal_f[nod_kept],
                       'b.', lw=0.1)
    d_kept = np.logical_and(do.SSNVs['judgement'] == 'KEEP', ~do.SSNVs.isnull()['failure_reasons']).values
    yerr_low = do.ssnv_based_model.normal_f[d_kept] - cis[0][d_kept]
    yerr_low[yerr_low<0] = 0
    detin_kept = ax.errorbar(do.ssnv_based_model.tumor_f[d_kept], do.ssnv_based_model.normal_f[d_kept],
                             yerr=[yerr_low,
                                   cis[1][d_kept] - do.ssnv_based_model.normal_f[d_kept]], fmt='r.', capsize=2)

    plt.xlabel('Tumor AF')
    plt.ylabel('Normal AF')
    plt.title('SSNVs considered and recovered')
    plt.legend(handles=[background[0], kept_def[0], detin_kept[0], TiN_fit[0]],
               labels=['Candidate Sites', 'Called w/o deTiN ', 'deTiN recovered', 'TiN_fit'])
    fig.set_dpi(300)
    fig.savefig(do.input.output_path + '/' + do.input.output_name + '_SSNVs_plot.png', bbox_inches='tight')
    fig.savefig(do.input.output_path + '/' + do.input.output_name + '_SSNVs_plot.eps', format='eps',
                bbox_inches='tight')


def select_candidate_mutations(call_stats_table, exac_db_file):
    # filter sites in call stats table to those only rejected for presence in the normal
    failure_reasons = np.array(call_stats_table['failure_reasons'])

    candidate_sites = call_stats_table[np.logical_or.reduce(np.array([np.array(call_stats_table['judgement']) == 'KEEP',
                                                                      failure_reasons == 'normal_lod,alt_allele_in_normal',
                                                                      failure_reasons == 'alt_allele_in_normal']))]
    candidate_sites['t_depth'] = candidate_sites['t_alt_count'] + candidate_sites['t_ref_count']
    candidate_sites['n_depth'] = candidate_sites['n_alt_count'] + candidate_sites['n_ref_count']
    candidate_sites.reset_index(inplace=True, drop=True)
    candidate_sites = remove_exac_sites_from_call_stats(candidate_sites, exac_db_file)

    candidate_sites.reset_index(inplace=True, drop=True)
    return candidate_sites


def hg19_to_linear_positions(chromosome, position, **keyword_parameters):
    # type: (nparray, nparray,string) -> nparray
    """
    Change chromosome-position to continuous linear coordinates
    """
    if ('build' in keyword_parameters):
        build = keyword_parameters['build']
    else:
        build = 'hg19'
    if build == 'hg19':
        L = np.array([249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                      146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                      102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                      155270560, 59373566, 16569])  # chromosome lengths from genome-mysql.cse.ucsc.edu
    if build == 'hg38':
        # add support for hg38
        sys.exit('support still missing for hg38')

    C = np.append(1, np.cumsum(L))
    x = np.array([chromosome[int(i)] for i in np.arange(0, len(position))])
    return C[[x.astype(int)]] + position


def fix_het_file_header(het_file):
    # allowing flexibility in het file headers to accommodate changing versions of GATK4 and other CN tools
    # in order to add support for your het file headers please modify the alternate header lists above
    headers = alternate_file_headers()

    required_headers = ['CONTIG', 'POSITION', 'ALT_COUNT', 'REF_COUNT']

    if np.sum(np.isfinite((is_member(required_headers, het_file.columns)))) == 4:
        return het_file
    else:
        missing_idx = np.where(~np.isfinite((is_member(required_headers, het_file.columns))))
        for i in missing_idx[0]:
            if required_headers[i] == 'POSITION':
                if np.sum(np.isfinite(is_member(headers['alternate_headers_position'], het_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(headers['alternate_headers_position'], het_file.columns))) > 1:
                    sys.exit('missing required header POSITION and could not replace with POS,position, or pos!')
                else:
                    idx_replace = np.where(np.isfinite(is_member(headers['alternate_headers_position'], het_file.columns)))
                    het_file.rename(columns={headers['alternate_headers_position'][idx_replace[0][0]]: 'POSITION'}, inplace=True)
                    print 'changing header of het file from ' + headers['alternate_headers_position'][
                        idx_replace[0][0]] + ' to POSITION'

            if required_headers[i] == 'CONTIG':
                if np.sum(np.isfinite(is_member(headers['alternate_headers_chromosome'], het_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(headers['alternate_headers_chromosome'], het_file.columns))) > 1:
                    sys.exit(
                        'missing required header CONTIG and could not replace with any one of CHR, chrom, Chromosome, chr, Chrom!')
                else:
                    idx_replace = np.where(np.isfinite(is_member(headers['alternate_headers_chromosome'], het_file.columns)))
                    het_file.rename(columns={headers['alternate_headers_chromosome'][idx_replace[0][0]]: 'CONTIG'}, inplace=True)
                    print 'changing header of het file from ' + headers['alternate_headers_chromosome'][
                        idx_replace[0][0]] + ' to CONTIG'

            if required_headers[i] == 'ALT_COUNT':
                if np.sum(np.isfinite(is_member(headers['alternate_headers_alt_count'], het_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(headers['alternate_headers_alt_count'], het_file.columns))) > 1:
                    sys.exit(
                        'missing required header ALT_COUNT and could not replace with any one of t_alt_count, n_alt_count, alt_count')
                else:
                    idx_replace = np.where(np.isfinite(is_member(headers['alternate_headers_alt_count'], het_file.columns)))
                    het_file.rename(columns={headers['alternate_headers_alt_count'][idx_replace[0][0]]: 'ALT_COUNT'}, inplace=True)
                    print 'changing header of het file from ' + headers['alternate_headers_alt_count'][
                        idx_replace[0][0]] + ' to ALT_COUNT'

            if required_headers[i] == 'REF_COUNT':
                if np.sum(np.isfinite(is_member(headers['alternate_headers_ref_count'], het_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(headers['alternate_headers_ref_count'], het_file.columns))) > 1:
                    sys.exit(
                        'missing required header ALT_COUNT and could not replace with any one of t_ref_count, n_ref_count, ref_count')
                else:
                    idx_replace = np.where(np.isfinite(is_member(headers['alternate_headers_ref_count'], het_file.columns)))
                    het_file.rename(columns={headers['alternate_headers_ref_count'][idx_replace[0][0]]: 'REF_COUNT'}, inplace=True)
                    print 'changing header of het file from ' + headers['alternate_headers_ref_count'][
                        idx_replace[0][0]] + ' to REF_COUNT'

        return het_file


def fix_seg_file_header(seg_file):
    # allowing flexibility in seg file headers to accommodate changing versions of GATK4 and other CN tools
    # in order to add support for your seg file headers please modify the alternate header lists above

    headers = alternate_file_headers()

    required_headers = ['Chromosome', 'Start.bp', 'End.bp', 'f', 'tau']

    if np.sum(np.isfinite((is_member(required_headers, seg_file.columns)))) == 5:
        return seg_file
    else:
        missing_idx = np.where(~np.isfinite((is_member(required_headers, seg_file.columns))))
        for i in missing_idx[0]:
            if required_headers[i] == 'Start.bp':
                if np.sum(np.isfinite(is_member(headers['alternate_headers_start_position'], seg_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(headers['alternate_headers_start_position'], seg_file.columns))) > 1:
                    sys.exit('missing required header Start.bp and could not replace with Start or Start_bp')
                else:
                    idx_replace = np.where(np.isfinite(is_member(headers['alternate_headers_start_position'], seg_file.columns)))
                    seg_file.rename(columns={headers['alternate_headers_start_position'][idx_replace[0][0]]: 'Start.bp'},
                                    inplace=True)
                    print 'changing header of seg file from ' + headers['alternate_headers_start_position'][
                        idx_replace[0][0]] + ' to Start.bp'

            if required_headers[i] == 'End.bp':
                if np.sum(np.isfinite(is_member(headers['alternate_headers_end_position'], seg_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(headers['alternate_headers_end_position'], seg_file.columns))) > 1:
                    sys.exit('missing required header End.bp and could not replace with End or End_bp')
                else:
                    idx_replace = np.where(np.isfinite(is_member(headers['alternate_headers_end_position'], seg_file.columns)))
                    seg_file.rename(columns={headers['alternate_headers_end_position'][idx_replace[0][0]]: 'End.bp'}, inplace=True)
                    print 'changing header of seg file from ' + headers['alternate_headers_end_position'][
                        idx_replace[0][0]] + ' to End.bp'

            if required_headers[i] == 'Chromosome':
                if np.sum(np.isfinite(is_member(headers['alternate_headers_chromosome'], seg_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(headers['alternate_headers_chromosome'], seg_file.columns))) > 1:
                    sys.exit(
                        'missing required header Chromosome and could not replace with any other header')
                else:
                    idx_replace = np.where(np.isfinite(is_member(headers['alternate_headers_chromosome'], seg_file.columns)))
                    seg_file.rename(columns={headers['alternate_headers_chromosome'][idx_replace[0][0]]: 'Chromosome'},
                                    inplace=True)
                    print 'changing header of seg file from ' + headers['alternate_headers_chromosome'][
                        idx_replace[0][0]] + ' to Chromosome'

            if required_headers[i] == 'f':
                if np.sum(np.isfinite(is_member(headers['alternate_headers_f'], seg_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(headers['alternate_headers_f'], seg_file.columns))) > 1:
                    sys.exit(
                        'missing required header f and could not replace with any one of f_acs')
                else:
                    idx_replace = np.where(np.isfinite(is_member(headers['alternate_headers_f'], seg_file.columns)))
                    seg_file.rename(columns={headers['alternate_headers_f'][idx_replace[0][0]]: 'f'}, inplace=True)
                    print 'changing header of het file from ' + headers['alternate_headers_f'][idx_replace[0][0]] + ' to f'

            if required_headers[i] == 'tau':
                if np.sum(np.isfinite(is_member(headers['alternate_headers_tau'], seg_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(headers['alternate_headers_tau'], seg_file.columns))) > 1:
                    sys.exit(
                        'missing required header tau and could not replace with any one of CN')
                else:
                    idx_replace = np.where(np.isfinite(is_member(headers['alternate_headers_tau'], seg_file.columns)))
                    seg_file.rename(columns={headers['alternate_headers_tau'][idx_replace[0][0]]: 'tau'}, inplace=True)
                    print 'changing header of het file from ' + headers['alternate_headers_tau'][idx_replace[0][0]] + ' to tau'

        return seg_file


def read_indel_vcf(vcf,seg_table,indel_type):
    # read strelka vcf
    # headerline should be in this format: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
    cols_type = {0: str}
    with open(vcf) as f:
        content = f.readlines()
    for line in content:
        if line[0] == '#' and line[1] != '#':
            headerline = line.split('\t')
            break

    if indel_type.lower() == 'strelka':
        indel_table = pd.read_csv(vcf, sep='\t', comment='#', header=None, low_memory=False, dtype=cols_type)
        indel_table.rename(columns={0: 'contig', 1: 'position',2:'ID',3:'REF',4:'ALT',5:'QUAL',7:'INFO', 8: 'format', 6: 'filter', 9: headerline[9].lower(), 10: headerline[10][0:-1].lower()},
                       inplace=True)
        counts_format = indel_table['format'][0].split(':')
        depth_ix = counts_format.index('DP')
        alt_indel_ix = counts_format.index('TIR')
        ref_indel_ix = counts_format.index('TAR')
        indel_table = indel_table[np.isfinite(is_member(indel_table['filter'], ['PASS', 'QSI_ref']))]
        indel_table.reset_index(inplace=True, drop=True)

    elif indel_type.lower() == 'mutect2':
        indel_table = pd.read_csv(vcf, sep='\t', comment='#', header=None, low_memory=False, dtype=cols_type)
        # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR	NORMAL
        normal_sample = 'normal'
        tumor_sample = 'tumor'
        for line in content:
            if line[0:15] == '##normal_sample':
                normal_sample = line.split('=')[1][0:-1]
            if line[0:14] == '##tumor_sample':
                tumor_sample = line.split('=')[1][0:-1]
        if tumor_sample == 'tumor' and normal_sample == 'normal':
            indel_table.rename(
                columns={0: 'contig', 1: 'position', 2: 'ID', 3: 'REF', 4: 'ALT', 5: 'QUAL', 7: 'INFO', 8: 'format',
                         6: 'filter', 9: 'tumor', 10: 'normal'},
                inplace=True)
        else:
            if tumor_sample == headerline[9]:
                indel_table.rename(
                        columns={0: 'contig', 1: 'position', 2: 'ID', 3: 'REF', 4: 'ALT', 5: 'QUAL', 7: 'INFO', 8: 'format',
                         6: 'filter', 9: 'tumor', 10: 'normal'},
                        inplace=True)
            elif tumor_sample == headerline[10][0:-1]:
                indel_table.rename(
                    columns={0: 'contig', 1: 'position', 2: 'ID', 3: 'REF', 4: 'ALT', 5: 'QUAL', 7: 'INFO', 8: 'format',
                             6: 'filter', 9: 'normal', 10: 'tumor'},
                    inplace=True)
            else:
                print 'failed to read MuTect 2 indels VCF'
                sys.exit()
        counts_format = indel_table['format'][0].split(':')
        depth_ix = counts_format.index('AD')
        indel_table = indel_table[np.isfinite(is_member(indel_table['filter'], ['PASS', 'alt_allele_in_normal','artifact_in_normal']))]
        indel_table.reset_index(inplace=True, drop=True)

    elif indel_type.lower() == 'sanger':
        indel_table = pd.read_csv(vcf, sep='\t', comment='#', header=None, low_memory=False, dtype=cols_type)
        # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOUR
        indel_table.rename(columns={0: 'contig', 1: 'position',2:'ID',3:'REF',4:'ALT',5:'QUAL',7:'INFO',8: 'format', 6: 'filter', 9: headerline[9].lower(), 10: headerline[10][0:-1].lower()},
                           inplace=True)
        b1 = np.logical_or.reduce([indel_table['filter'] == 'F012', indel_table['filter'] == 'F012;F015'])
        b2 = np.logical_or.reduce([indel_table['filter'] == 'PASS', indel_table['filter'] == 'F015'])
        indel_table = indel_table[np.logical_or.reduce([b1, b2])]
        indel_table.reset_index(inplace=True,drop=True)
        format_string = indel_table['format'][0].split(':')
        total_depth_idx = [format_string.index('PR'), format_string.index('NR')]
        alt_count_idx = [format_string.index('PU'), format_string.index('NU')]

    # parsing format line and file to determine required alt and ref columns
    # we use "tier 1" read counts for varaints
    n_depth = np.zeros([len(indel_table), 1])
    n_alt_count = np.zeros([len(indel_table), 1])
    n_ref_count = np.zeros([len(indel_table), 1])

    t_depth = np.zeros([len(indel_table), 1])
    t_alt_count = np.zeros([len(indel_table), 1])
    t_ref_count = np.zeros([len(indel_table), 1])

    for index, row in indel_table.iterrows():
        spl_n = row['normal'].split(':')
        spl_t = row['tumor'].split(':')
        if indel_type.lower() == 'strelka':
            n_depth[index] = int(spl_n[depth_ix])
            n_alt_count[index] = int(spl_n[alt_indel_ix].split(',')[0])
            n_ref_count[index] = int(spl_n[ref_indel_ix].split(',')[0])
            t_depth[index] = int(spl_t[depth_ix])
            t_alt_count[index] = int(spl_t[alt_indel_ix].split(',')[0])
            t_ref_count[index] = int(spl_t[ref_indel_ix].split(',')[0])
        if indel_type.lower() == 'mutect2':
            n_alt_count[index] = int(spl_n[depth_ix].split(',')[1])
            n_ref_count[index] = int(spl_n[depth_ix].split(',')[0])
            n_depth[index] = n_alt_count[index]+n_ref_count[index]
            t_alt_count[index] = int(spl_t[depth_ix].split(',')[1])
            t_ref_count[index] = int(spl_t[depth_ix].split(',')[0])
            t_depth[index] = t_alt_count[index] + t_ref_count[index]
        if indel_type.lower() == 'sanger':
            n_depth[index] = np.sum([int(spl_n[i]) for i in total_depth_idx])
            n_alt_count[index] = np.sum([int(spl_n[i]) for i in alt_count_idx])
            n_ref_count[index] = n_depth[index] - n_alt_count[index]
            t_depth[index] = np.sum([int(spl_t[i]) for i in total_depth_idx])
            t_alt_count[index] = np.sum([int(spl_t[i]) for i in alt_count_idx])
            t_ref_count[index] = t_depth[index] - t_alt_count[index]
    if len(indel_table) == 0:
        indel_table = pd.DataFrame(index=[0],columns=['contig', 'position','ID','REF','ALT','QUAL','INFO','format', 'filter',headerline[9].lower(), headerline[10][0:-1].lower(),
                                                      't_depth','t_alt_count','t_ref_count','n_alt_count','n_depth','n_ref_count','tau','f_acs','Chromosome','genomic_coord_x'])
    else:
        indel_table['t_depth'] = t_alt_count + t_ref_count
        indel_table['t_alt_count'] = t_alt_count
        indel_table['t_ref_count'] = t_ref_count

        indel_table['n_depth'] = n_alt_count + n_ref_count
        indel_table['n_alt_count'] = n_alt_count
        indel_table['n_ref_count'] = n_ref_count
    # only consider sites which were rejected as germline or were passed

        if type(indel_table['contig'][0]) == str :
            indel_table['Chromosome'] = chr2num(indel_table['contig'])
        else:
            indel_table['Chromosome'] = indel_table['contig']-1
    # add linear position field and consider only sites which are rejected as germline i.e. PASS or QSI_ref
        indel_table = indel_table[np.isfinite(indel_table['Chromosome'])]
        indel_table.reset_index(inplace=True, drop=True)
        indel_table['genomic_coord_x'] = hg19_to_linear_positions(indel_table['Chromosome'], indel_table['position'])
    # annotate with acs data
        f_acs = np.zeros([len(indel_table), 1]) + 0.5
        tau = np.zeros([len(indel_table), 1]) + 2
        for i, r in seg_table.iterrows():
            f_acs[np.logical_and(np.array(indel_table['genomic_coord_x']) >= r['genomic_coord_start'],
                             np.array(indel_table['genomic_coord_x']) <= r['genomic_coord_end'])] = r.f
            tau[np.logical_and(np.array(indel_table['genomic_coord_x']) >= r['genomic_coord_start'],
                           np.array(indel_table['genomic_coord_x']) <= r[
                               'genomic_coord_end'])] = r.tau + 0.001
        indel_table['tau'] = tau
        indel_table['f_acs'] = f_acs

    return indel_table

def build_exac_pickle(exac_file):
    # create ExAC site dictionary from VCF file
    exac_site_info = {}
    print 'Filtering ExAC sites from candidate mutations'
    with gzip.open(exac_file, "rb") as vcf_file:
        for line_index, line in enumerate(vcf_file):
            if line_index % 10000 == 0:
                print 'processed ' + str(line_index) + ' ExAC sites'
            spl = line.strip("\n").split("\t")

            # line is a comment
            if line[0] == '#':
                continue
            site = spl[0] + '_' + spl[1]
            info = spl[7]
            info_dict = {}
            info_dict['ref_allele'] = spl[3]
            info_dict['alt_allele'] = spl[4]
            for field in info.strip("\n").split(";"):
                if field.split("=")[0] not in ['AC', 'AF']: continue
                try:
                    info_dict[field.split("=")[0]] = field.split("=")[1]
                except:
                    pass
                    # print 'bad field:', field
            # select only sites where population allele fractions exceeds 0.01
            if np.sum(np.array(info_dict['AF'].split(','), dtype=float)) >= 0.01:
                exac_site_info[site] = info_dict
    with open('exac.pickle', 'wb') as handle:
        pickle.dump(exac_site_info, handle, protocol=pickle.HIGHEST_PROTOCOL)


def remove_exac_sites_from_call_stats(call_stats_table, exac_file):
    # use ExAC vcf to filter likely germline variants out of candidate sites
    with open(exac_file, 'rb') as handle:
        exac_dict = pickle.load(handle)
        keep = np.ones_like(call_stats_table['position'], dtype=bool)
    for index, row in call_stats_table.iterrows():
        key = str(row['contig']) + '_' + str(row['position'])
        try:
            exac_dict[key]
            # print 'removing site '+ key+ ' minor allele fraction = ' + str(exac_dict[key]['AF'])
            keep[index] = False
        except KeyError:
            pass
    return call_stats_table[keep]
