## Synopsis

DeTiN estimates tumor in normal (TiN) based on tumor and matched normal sequencing data. The estimate is based on both candidate SSNVs and aSCNAs. DeTiN then applies the joint TiN estimate to reclassify SSNVs and InDels as somatic or germline. Install and run time on standard exome data is about 5 mins. 

## Code Example
Please see github wiki for description of input files. 

python deTiN.py --mutation_data_path example_data/HCC_10_90.call_stats.pon_removed.txt --cn_data_path example_data/HCC-1143_100_T-sim-final.acs.seg --tumor_het_data example_data/HCC_10_90.tumor.hets.tsv --normal_het_data example_data/HCC_10_90.normal.hets.tsv --exac_data_path example_data/exac.pickle_high_af --output_name 10_percent_TiN_simulation --indel_data_path example_data/MuTect2.call_stats.txt --indel_data_type MuTect2 --output_dir example_data/

## Parameter descriptions

See project Wiki for full description of required fields for input data.

Input data:

–-mutation_data_path mutation statistics file (MuTect call stats file (or similar variants file)).

–-cn_data_path allelic copy number segmentation file (GATK4 AllelicCNV seg file).

–-tumor_het_data heterozygous SNP variant counts in the tumor sample. (GATK4 tumor het cov file).

–-normal_het_data heterozygous SNP variant counts in the normal sample. (GATK4 normal het cov file).

–-exac_data_path pickle file of minor allele fraction > 0.01 ExAC sites. 

Parameters:

–-output_name sample name

–-output_dir output directory

Optional parameters:

--TiN_prior (default = 0.5)
0.5 is a null prior. If users wish to require more evidence for TiN > 0 this can be lowered. 

--mutation_prior (default = 0.15)
The ratio of sites expected to be mutated somatically to rare germline events (0.15 corresponds to ~2 mutations per megabase)

--ascna_probe_number_filter (default = 200)
We require 200 probes based on empirical results using GATK4CNV that segments smaller than this tend to be enriched for artifacts. For WGS this parameter can be set to 0.

--ascna_SNP_number_filter (default = 20)
We require 20 SNPs based on empirical results using GATK4CNV that segments smaller than this tend to be enriched for artifacts. 

--coverage_threshold (default = 15)
Number of reads required to use a variant for estimation of TiN. We require 15x coverage since low coverage sites tend to be enriched for artifacts. (NOTE: all sites are considered for somatic recovery)

--SSNV_af_threshold (default = 0.15)
Fraction of alternate allele required for site to be used in SSNV based estimation of TiN. We require 15% since low af sites tend to be enriched for artifacts. If users are using more deeply sequenced data this should be set to a lower value.
(NOTE: all sites are considered for somatic recovery)

--aSCNA_threshold (default = 0.1)
Fraction of allele shift required to use a segment for TiN estimation. Lower this value for extremely well covered samples (e.g. 500x). 

--aSCNA_variance_threshold (default = 0.025)
Variance tolerated in allele shift of a segment before removal. This filter helps to remove regions enriched for artifact sites such as centromeres/telomeres and low mapping regions. 

--cancer_hot_spots
Optional BED file of cancer hot spot mutations which the user has a stronger prior on being somatic e.g. BRAF v600E mutations.

## Motivation
Genomic characterization is vital to the understanding and treatment of cancer.  Detection of somatic mutations is a critical component of this process. A key step in sensitive and specific somatic mutation detection is comparison of the tumor sample to a matched germline control. Sensitivity to detect somatic variants is greatly reduced when the matched normal sample is contaminated with tumor cells. To overcome this limitation, we developed deTiN, a method that estimates tumor-in-normal contamination (TiN), and improves detection sensitivity when using a contaminated normal. 

## Installation

deTiN requires Numpy, Pandas, Scipy and Python 2.7. 

pip install deTiN

or

git clone https://github.com/broadinstitute/deTiN.git

cd deTiN/

python setup install


## Example Data

The above code example will run using the included data from an artifically mixed 10% contaminated normal. Input files were generated using MuTect and GATK4ACNV. 


## License

Copyright 2017 Amaro Taylor-Weiner

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
