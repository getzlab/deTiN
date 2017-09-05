import numpy as np
import sys
import gzip


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


def load_exac(exac_vcf):
    return


def identify_aSCNAs(deTiN_input):
    return


def plotting(deTiN_output):
    return


def select_candidate_mutations(call_stats_table):
    failure_reasons = np.array(call_stats_table['failure_reasons'])
    candidate_sites = call_stats_table[np.logical_or.reduce(np.array([np.array(call_stats_table['judgement']) == 'KEEP',
                                                                      failure_reasons == 'normal_lod,alt_allele_in_normal',
                                                                      failure_reasons == 'normal_lod',
                                                                      failure_reasons == 'alt_allele_in_normal']))]
    candidate_sites.reset_index(inplace=True)
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
