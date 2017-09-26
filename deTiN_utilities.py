import numpy as np
import sys
from scipy.stats import beta
from scipy.stats import fisher_exact
from itertools import compress
import random
import pandas as pd
import matplotlib

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


def filter_hets_based_on_coverage(het_table):
    het_table = het_table[np.logical_and(het_table['READ_DEPTH_N'] > 10, het_table['READ_DEPTH_T'] > 10)]
    het_table.reset_index(inplace=True, drop=True)
    return het_table


def filter_segments_based_on_size_f_and_tau(seg_table):
    seg_table = seg_table[np.logical_and.reduce(np.array([np.array(seg_table['f']) < 0.40,
                                                          seg_table['n_probes'] > 200, seg_table['tau'] > 0]))]
    seg_table.reset_index(inplace=True, drop=True)
    return seg_table


def load_exac(exac_vcf):
    return


def identify_aSCNAs(seg_table, het_table):
    # identify aSCNAs based on minor allele fraction of segments
    mu_af_n = np.mean(het_table['AF_N'])
    f_detin = np.zeros([len(seg_table), 1])
    f_variance = np.zeros([len(seg_table), 1])
    n_snps_above_mu = np.zeros([len(seg_table), 1])
    n_snps_below_mu = np.zeros([len(seg_table), 1])
    fishers_p_convergent_seg = np.ones([len(seg_table), 1])

    for seg_id, seg in seg_table.iterrows():
        seg_hets = het_table[het_table['seg_id'] == seg_id]
        f_detin[seg_id] = mu_af_n - np.mean(np.abs(seg_hets['AF_T'] - mu_af_n))
        f_variance[seg_id] = np.var(np.abs(seg_hets['AF_T'] - mu_af_n))
        n_snps_above_mu[seg_id] = np.mean([np.sum(seg_hets['AF_T'] > mu_af_n),
                                           np.sum(seg_hets['AF_N']) > mu_af_n])
        n_snps_below_mu[seg_id] = np.mean([np.sum(seg_hets['AF_T'] <= mu_af_n),
                                           np.sum(seg_hets['AF_N']) <= mu_af_n])
        fe_tuple = fisher_exact([[np.sum(np.logical_and(seg_hets['AF_T'] > mu_af_n,
                                                        seg_hets['AF_N'] > mu_af_n)),
                                  np.sum(np.logical_and(seg_hets['AF_T'] > mu_af_n,
                                                        seg_hets['AF_N'] <= mu_af_n))],
                                 [np.sum(np.logical_and(seg_hets['AF_T'] <= mu_af_n,
                                                        seg_hets['AF_N'] > mu_af_n)),
                                  np.sum(np.logical_and(seg_hets['AF_T'] <= mu_af_n,
                                                        seg_hets['AF_N'] <= mu_af_n))]], 'less')
        fishers_p_convergent_seg[seg_id] = fe_tuple[1]
    seg_table['f_detin'] = f_detin
    seg_table['f_variance'] = f_variance
    seg_table['n_snps_above_mu'] = n_snps_above_mu
    seg_table['n_snps_below_mu'] = n_snps_below_mu
    seg_table['fishers_p_convergent_seg'] = fishers_p_convergent_seg
    if any((seg_table['fishers_p_convergent_seg'] * len(seg_table)) < 0.05):
        segs = (seg_table['fishers_p_convergent_seg'] * len(seg_table)) < 0.05
        ix = list(compress(xrange(len(segs)), segs))
        print 'identified convergent aSCNA in normal on chromosomes:' + str(np.unique(seg_table['Chromosome'][ix] + 1))
    aSCNAs = seg_table[
        np.logical_and.reduce(np.array([np.array(seg_table['fishers_p_convergent_seg'] * len(seg_table)) > 0.05,
                                        seg_table['n_snps_above_mu'] > 10,
                                        seg_table['n_snps_below_mu'] > 10,
                                        seg_table['f_detin'] <= 0.4,
                                        seg_table['f_variance'] < 0.025]))]
    return aSCNAs


def ensure_balanced_hets(seg_table, het_table):
    seg_table['aSCNA'] = np.zeros([len(seg_table), 1])
    aSCNA_hets = []
    for seg_id, seg in seg_table.iterrows():
        seg_hets = het_table[het_table['seg_id'] == seg_id]
        if np.sum(seg_hets['d'] == -1) > 10 and np.sum(seg_hets['d'] == 1) > 10:
            alts = np.concatenate([np.array(seg_hets['ALT_COUNT_T'][np.array(seg_hets['d'] == -1)]),
                                   np.array(seg_hets['REF_COUNT_T'][np.array(seg_hets['d'] == 1)])])
            refs = np.concatenate([np.array(seg_hets['ALT_COUNT_T'][np.array(seg_hets['d'] == 1)]),
                                   np.array(seg_hets['REF_COUNT_T'][np.array(seg_hets['d'] == -1)])])

            f = np.mean(np.true_divide(alts, alts + refs))
            seg_hets = seg_hets[
                np.logical_and(beta.sf(f, alts + 1, refs + 1) < 0.995, beta.sf(f, alts + 1, refs + 1) > 0.005)]
            if sum(seg_hets['AF_N'] > 0.5) < sum(seg_hets['AF_N'] <= 0.5):
                sites = seg_hets['AF_N'] <= 0.5
                index = list(compress(xrange(len(sites)), sites))
                ixs = random.sample(index, (sum(seg_hets['AF_N'] <= 0.5) - sum(seg_hets['AF_N'] > 0.5)))
                seg_hets = seg_hets.drop(seg_hets.index[[ixs]])
                seg_hets.reset_index(inplace=True, drop=True)

            if sum(seg_hets['AF_N'] > 0.5) > sum(seg_hets['AF_N'] <= 0.5):
                sites = seg_hets['AF_N'] > 0.5
                index = list(compress(xrange(len(sites)), sites))
                ixs = random.sample(index, (sum(seg_hets['AF_N'] > 0.5) - sum(seg_hets['AF_N'] <= 0.5)))
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
    fig.set_dpi(150)
    fig.savefig(output_path + '/' + sample_name + '_KmeansEval_plot.png', bbox_inches='tight')

    # scatter plot of TiN estimates per segment by chromosome location and cluster
    fig = plt.figure()
    ax = fig.add_subplot(111)
    clr = ['b', 'g', 'r']
    for i in range(K[kIdx]):
        ind = (ascna_based_model.cluster_assignment == i)
        ax.scatter(X[ind], Y[ind] + 1, s=30, c=clr[i], label='Cluster %d' % i)
    plt.xlabel('MAP tumor in normal estimate (%)')
    plt.ylabel('Chromosome')
    plt.title('Cluster by chromosome and TiN')
    plt.yticks(np.arange(min(Y) + 1, max(Y) + 2, 2.0))
    plt.xticks(np.arange(0, max(X) + 1, np.max([np.round(np.true_divide(np.max(X), 10)), 1])))
    ax.set_xlim([-2, np.max(X) + 2])

    fig.set_dpi(150)
    fig.savefig(output_path + '/' + sample_name + '_KmeansEval_scatter_plot.png', bbox_inches='tight')


def plot_TiN_models(do):
    fig, ax = plt.subplots(1, 1)
    ascna = ax.plot(do.ascna_based_model.TiN_range,
                    np.exp(do.ascna_based_model.TiN_likelihood - np.nanmax(do.ascna_based_model.TiN_likelihood))
                    , 'r--', lw=1)
    ssnv = ax.plot(do.ascna_based_model.TiN_range,
                   np.exp(do.ssnv_based_model.TiN_likelihood - np.nanmax(do.ssnv_based_model.TiN_likelihood))
                   , 'b--', lw=1)

    joint = ax.plot(do.ascna_based_model.TiN_range, do.joint_posterior
                    , 'k-', lw=2)
    plt.xlabel('Tumor in normal estimate')
    plt.ylabel('p(TiN=x)')
    plt.title('TiN estimate posterior')
    plt.legend(handles=[ascna[0], ssnv[0], joint[0]], labels=['aSCNA', 'SSNV', 'Joint Est.'])
    fig.set_dpi(150)
    fig.savefig(do.input.output_path + '/' + do.input.output_name + '_TiN_models_plot.png', bbox_inches='tight')


def plot_SSNVs(do):
    fig, ax = plt.subplots(1, 1)
    TiN_fit = ax.plot(np.linspace(0, 1, 101), np.multiply(do.TiN, np.linspace(0, 1, 101)), '--', lw=1, alpha=1,
                      color='#1D1D1D')
    background = ax.plot(do.ssnv_based_model.tumor_f, do.ssnv_based_model.normal_f
                         , '.', lw=0.1, alpha=0.75, color=[0.75, 0.75, 0.75])
    nod_kept = np.logical_and(do.SSNVs['judgement'] == 'KEEP', do.SSNVs.isnull()['failure_reasons'])
    cis = do.ssnv_based_model.rv_normal_af.interval(0.682)

    kept_def = ax.plot(do.ssnv_based_model.tumor_f[nod_kept], do.ssnv_based_model.normal_f[nod_kept],
                       'b.', lw=0.1)
    d_kept = np.logical_and(do.SSNVs['judgement'] == 'KEEP', ~do.SSNVs.isnull()['failure_reasons'])
    detin_kept = ax.errorbar(do.ssnv_based_model.tumor_f[d_kept], do.ssnv_based_model.normal_f[d_kept],
                             yerr=[do.ssnv_based_model.normal_f[d_kept] - cis[0][d_kept],
                                   cis[1][d_kept] - do.ssnv_based_model.normal_f[d_kept]], fmt='r.', capsize=2)

    plt.xlabel('Tumor AF')
    plt.ylabel('Normal AF')
    plt.title('SSNVs considered and recovered')
    plt.legend(handles=[background[0], kept_def[0], detin_kept[0], TiN_fit[0]],
               labels=['Candidate Sites', 'Called w/o deTiN ', 'deTiN recovered', 'TiN_fit'])
    fig.set_dpi(150)
    fig.savefig(do.input.output_path + '/' + do.input.output_name + '_SSNVs_plot.png', bbox_inches='tight')


def select_candidate_mutations(call_stats_table):
    # filter sites in call stats table to those only rejected for presence in the normal
    failure_reasons = np.array(call_stats_table['failure_reasons'])

    candidate_sites = call_stats_table[np.logical_or.reduce(np.array([np.array(call_stats_table['judgement']) == 'KEEP',
                                                                      failure_reasons == 'normal_lod,alt_allele_in_normal',
                                                                      failure_reasons == 'normal_lod',
                                                                      failure_reasons == 'alt_allele_in_normal']))]
    candidate_sites['t_depth'] = candidate_sites['t_alt_count'] + candidate_sites['t_ref_count']
    candidate_sites['n_depth'] = candidate_sites['n_alt_count'] + candidate_sites['n_ref_count']


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
    # in order to add support for your het file headers please modify the alternate header lists below

    alternate_headers_position = ['POS', 'position', 'pos']
    alternate_headers_chromosome = ['CHR', 'chrom', 'Chromosome', 'chr', 'Chrom']
    alternate_headers_alt_count = ['t_alt_count', 'n_alt_count', 'alt_count']
    alternate_headers_ref_count = ['t_ref_count', 'n_ref_count', 'ref_count']

    required_headers = ['CONTIG', 'POSITION', 'ALT_COUNT', 'REF_COUNT']

    if np.sum(np.isfinite((is_member(required_headers, het_file.columns)))) == 4:
        return het_file
    else:
        missing_idx = np.where(~np.isfinite((is_member(required_headers, het_file.columns))))
        for i in missing_idx[0]:
            if required_headers[i] == 'POSITION':
                if np.sum(np.isfinite(is_member(alternate_headers_position, het_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(alternate_headers_position, het_file.columns))) > 1:
                    sys.exit('missing required header POSITION and could not replace with POS,position, or pos!')
                else:
                    idx_replace = np.where(np.isfinite(is_member(alternate_headers_position, het_file.columns)))
                    het_file.rename(columns={alternate_headers_position[idx_replace[0][0]]: 'POSITION'}, inplace=True)
                    print 'changing header of het file from ' + alternate_headers_position[
                        idx_replace[0][0]] + ' to POSITION'

            if required_headers[i] == 'CONTIG':
                if np.sum(np.isfinite(is_member(alternate_headers_chromosome, het_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(alternate_headers_chromosome, het_file.columns))) > 1:
                    sys.exit(
                        'missing required header CONTIG and could not replace with any one of CHR, chrom, Chromosome, chr, Chrom!')
                else:
                    idx_replace = np.where(np.isfinite(is_member(alternate_headers_chromosome, het_file.columns)))
                    het_file.rename(columns={alternate_headers_chromosome[idx_replace[0][0]]: 'CONTIG'}, inplace=True)
                    print 'changing header of het file from ' + alternate_headers_chromosome[
                        idx_replace[0][0]] + ' to CONTIG'

            if required_headers[i] == 'ALT_COUNT':
                if np.sum(np.isfinite(is_member(alternate_headers_alt_count, het_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(alternate_headers_alt_count, het_file.columns))) > 1:
                    sys.exit(
                        'missing required header ALT_COUNT and could not replace with any one of t_alt_count, n_alt_count, alt_count')
                else:
                    idx_replace = np.where(np.isfinite(is_member(alternate_headers_alt_count, het_file.columns)))
                    het_file.rename(columns={alternate_headers_alt_count[idx_replace[0][0]]: 'ALT_COUNT'}, inplace=True)
                    print 'changing header of het file from ' + alternate_headers_alt_count[
                        idx_replace[0][0]] + ' to ALT_COUNT'

            if required_headers[i] == 'REF_COUNT':
                if np.sum(np.isfinite(is_member(alternate_headers_ref_count, het_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(alternate_headers_ref_count, het_file.columns))) > 1:
                    sys.exit(
                        'missing required header ALT_COUNT and could not replace with any one of t_ref_count, n_ref_count, ref_count')
                else:
                    idx_replace = np.where(np.isfinite(is_member(alternate_headers_ref_count, het_file.columns)))
                    het_file.rename(columns={alternate_headers_ref_count[idx_replace[0][0]]: 'REF_COUNT'}, inplace=True)
                    print 'changing header of het file from ' + alternate_headers_ref_count[
                        idx_replace[0][0]] + ' to REF_COUNT'

        return het_file


def fix_seg_file_header(seg_file):
    # allowing flexibility in seg file headers to accommodate changing versions of GATK4 and other CN tools
    # in order to add support for your seg file headers please modify the alternate header lists below

    alternate_headers_start_position = ['Start', 'Start_bp', 'start']
    alternate_headers_end_position = ['End', 'End_bp', 'end']
    alternate_headers_chromosome = ['Contig', 'chrom', 'CONTIG', 'chr', 'Chrom', 'CHROMOSOME']
    alternate_headers_f = ['f_acs','MAF_Post_Mode']
    alternate_headers_tau = ['CN','Segment_Mean_Post_Mode']

    required_headers = ['Chromosome', 'Start.bp', 'End.bp', 'f', 'tau']

    if np.sum(np.isfinite((is_member(required_headers, seg_file.columns)))) == 5:
        return seg_file
    else:
        missing_idx = np.where(~np.isfinite((is_member(required_headers, seg_file.columns))))
        for i in missing_idx[0]:
            if required_headers[i] == 'Start.bp':
                if np.sum(np.isfinite(is_member(alternate_headers_start_position, seg_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(alternate_headers_start_position, seg_file.columns))) > 1:
                    sys.exit('missing required header Start.bp and could not replace with Start or Start_bp')
                else:
                    idx_replace = np.where(np.isfinite(is_member(alternate_headers_start_position, seg_file.columns)))
                    seg_file.rename(columns={alternate_headers_start_position[idx_replace[0][0]]: 'Start.bp'},
                                    inplace=True)
                    print 'changing header of seg file from ' + alternate_headers_start_position[
                        idx_replace[0][0]] + ' to Start.bp'

            if required_headers[i] == 'End.bp':
                if np.sum(np.isfinite(is_member(alternate_headers_end_position, seg_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(alternate_headers_end_position, seg_file.columns))) > 1:
                    sys.exit('missing required header End.bp and could not replace with End or End_bp')
                else:
                    idx_replace = np.where(np.isfinite(is_member(alternate_headers_end_position, seg_file.columns)))
                    seg_file.rename(columns={alternate_headers_end_position[idx_replace[0][0]]: 'End.bp'}, inplace=True)
                    print 'changing header of seg file from ' + alternate_headers_end_position[
                        idx_replace[0][0]] + ' to End.bp'

            if required_headers[i] == 'Chromosome':
                if np.sum(np.isfinite(is_member(alternate_headers_chromosome, seg_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(alternate_headers_chromosome, seg_file.columns))) > 1:
                    sys.exit(
                        'missing required header Chromosome and could not replace with any other header')
                else:
                    idx_replace = np.where(np.isfinite(is_member(alternate_headers_chromosome, seg_file.columns)))
                    seg_file.rename(columns={alternate_headers_chromosome[idx_replace[0][0]]: 'Chromosome'},
                                    inplace=True)
                    print 'changing header of seg file from ' + alternate_headers_chromosome[
                        idx_replace[0][0]] + ' to Chromosome'

            if required_headers[i] == 'f':
                if np.sum(np.isfinite(is_member(alternate_headers_f, seg_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(alternate_headers_f, seg_file.columns))) > 1:
                    sys.exit(
                        'missing required header f and could not replace with any one of f_acs')
                else:
                    idx_replace = np.where(np.isfinite(is_member(alternate_headers_f, seg_file.columns)))
                    seg_file.rename(columns={alternate_headers_f[idx_replace[0][0]]: 'f'}, inplace=True)
                    print 'changing header of het file from ' + alternate_headers_f[idx_replace[0][0]] + ' to f'

            if required_headers[i] == 'tau':
                if np.sum(np.isfinite(is_member(alternate_headers_tau, seg_file.columns))) == 0 or np.sum(
                        np.isfinite(is_member(alternate_headers_tau, seg_file.columns))) > 1:
                    sys.exit(
                        'missing required header tau and could not replace with any one of CN')
                else:
                    idx_replace = np.where(np.isfinite(is_member(alternate_headers_tau, seg_file.columns)))
                    seg_file.rename(columns={alternate_headers_tau[idx_replace[0][0]]: 'tau'}, inplace=True)
                    print 'changing header of het file from ' + alternate_headers_tau[idx_replace[0][0]] + ' to tau'

        return seg_file
