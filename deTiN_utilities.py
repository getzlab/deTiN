import numpy as np
import sys
from scipy.stats import beta
from scipy.stats import fisher_exact
from itertools import compress
import random
import pandas as pd


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def is_member(a, b):
    # type: (nparray, nparray) -> nparray
    # based on the matlab is_member function
    # code from stackoverflow user Joe Kington
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return [bind.get(itm, np.nan) for itm in a]


def chr2num(chr):
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
    distance_centromere = np.zeros([len(positions),len(centromere_positions)])
    distance_telomere = np.zeros([len(positions),len(telomere_positions)])

    for i,centromere in enumerate(centromere_positions):
        distance_centromere[:,i] = np.abs(positions - centromere)
    distance_centromere = np.min(distance_centromere,axis=1)

    for i,telomere in enumerate(telomere_positions):
        distance_telomere[:,i] = np.abs(positions - telomere)
    distance_telomere = np.min(distance_telomere,axis=1)
    het_table = het_table[np.logical_and(distance_centromere>5000000,distance_telomere>5000000)]
    het_table.reset_index(inplace=True,drop=True)
    return het_table

def filter_hets_based_on_coverage(het_table):

    het_table=het_table[np.logical_and(het_table['READ_DEPTH_N']>10,het_table['READ_DEPTH_T']>10)]
    het_table.reset_index(inplace=True,drop=True)
    return het_table

def filter_segments_based_on_size_f_and_tau(seg_table):

    seg_table = seg_table[np.logical_and.reduce(np.array([np.array(seg_table['f']) < 0.40,
                                                seg_table['n_probes'] > 200,seg_table['tau'] > 0]))]
    seg_table.reset_index(inplace=True,drop=True)
    return seg_table

def load_exac(exac_vcf):
    return

def identify_aSCNAs(seg_table,het_table):

    mu_af_n = np.mean(het_table['AF_N'])
    f_detin = np.zeros([len(seg_table),1])
    f_variance = np.zeros([len(seg_table), 1])
    n_snps_above_mu = np.zeros([len(seg_table),1])
    n_snps_below_mu = np.zeros([len(seg_table),1])
    fishers_p_convergent_seg = np.ones([len(seg_table),1])

    for seg_id,seg in seg_table.iterrows():
        seg_hets = het_table[het_table['seg_id'] == seg_id]
        f_detin[seg_id] = mu_af_n - np.mean(np.abs(seg_hets['AF_T']-mu_af_n))
        f_variance[seg_id] = np.var(np.abs(seg_hets['AF_T']-mu_af_n))
        n_snps_above_mu[seg_id] = np.mean([np.sum(seg_hets['AF_T'] > mu_af_n),
                                                        np.sum(seg_hets['AF_N'])> mu_af_n])
        n_snps_below_mu[seg_id] = np.mean([np.sum(seg_hets['AF_T'] <= mu_af_n),
                                                        np.sum(seg_hets['AF_N'])<= mu_af_n])
        fe_tuple = fisher_exact([[np.sum(np.logical_and(seg_hets['AF_T']>mu_af_n,
                                                                                seg_hets['AF_N']>mu_af_n)),
                                                          np.sum(np.logical_and(seg_hets['AF_T']>mu_af_n,
                                                                                seg_hets['AF_N']<=mu_af_n))],
                                                         [np.sum(np.logical_and(seg_hets['AF_T']<=mu_af_n,
                                                                                seg_hets['AF_N']>mu_af_n)),
                                                          np.sum(np.logical_and(seg_hets['AF_T']<=mu_af_n,
                                                                                seg_hets['AF_N']<=mu_af_n))]],'less')
        fishers_p_convergent_seg[seg_id] = fe_tuple[1]
    seg_table['f_detin'] = f_detin
    seg_table['f_variance'] = f_variance
    seg_table['n_snps_above_mu'] = n_snps_above_mu
    seg_table['n_snps_below_mu'] = n_snps_below_mu
    seg_table['fishers_p_convergent_seg'] = fishers_p_convergent_seg
    if any((seg_table['fishers_p_convergent_seg']*len(seg_table))<0.05):
        segs = (seg_table['fishers_p_convergent_seg'] * len(seg_table)) < 0.05
        ix = list(compress(xrange(len(segs)), segs))
        print 'identified convergent aSCNA in normal on chromosomes:' + str(np.unique(seg_table['Chromosome'][ix]+1))
    aSCNAs = seg_table[np.logical_and.reduce(np.array([np.array(seg_table['fishers_p_convergent_seg']*len(seg_table)) >0.05,
                                                                      seg_table['n_snps_above_mu'] > 10,
                                                                      seg_table['n_snps_below_mu'] > 10,
                                                                      seg_table['f_detin']<= 0.4,
                                                                      seg_table['f_variance']< 0.025]))]
    return aSCNAs

def ensure_balanced_hets(seg_table,het_table):
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
            if sum(seg_hets['AF_N'] > 0.5) < sum(seg_hets['AF_N'] <= 0.5) :
                sites = seg_hets['AF_N'] <= 0.5
                index = list(compress(xrange(len(sites)), sites))
                ixs = random.sample(index, (sum(seg_hets['AF_N'] <= 0.5) - sum(seg_hets['AF_N'] > 0.5)))
                seg_hets = seg_hets.drop(seg_hets.index[[ixs]])
                seg_hets.reset_index(inplace=True,drop=True)

            if sum(seg_hets['AF_N'] > 0.5) > sum(seg_hets['AF_N'] <= 0.5):
                sites = seg_hets['AF_N'] > 0.5
                index = list(compress(xrange(len(sites)), sites))
                ixs = random.sample(index, (sum(seg_hets['AF_N'] > 0.5) - sum(seg_hets['AF_N'] <= 0.5)))
                seg_hets = seg_hets.drop(seg_hets.index[[ixs]])
                seg_hets.reset_index(inplace=True, drop=True)
            if len(aSCNA_hets)==0:
                aSCNA_hets = seg_hets
            else :
                aSCNA_hets=pd.concat([aSCNA_hets,seg_hets])
                aSCNA_hets.reset_index(inplace=True,drop=True)
    return aSCNA_hets

def plotting(deTiN_output):
    return


def select_candidate_mutations(call_stats_table):
    failure_reasons = np.array(call_stats_table['failure_reasons'])
    candidate_sites = call_stats_table[np.logical_or.reduce(np.array([np.array(call_stats_table['judgement']) == 'KEEP',
                                                                      failure_reasons == 'normal_lod,alt_allele_in_normal',
                                                                      failure_reasons == 'normal_lod',
                                                                      failure_reasons == 'alt_allele_in_normal']))]
    candidate_sites.reset_index(inplace=True,drop=True)
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
