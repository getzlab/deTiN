## Synopsis

DeTiN estimates tumor in normal (TiN) based on tumor and matched normal sequencing data. The estimate is based on both candidate SSNVs and aSCNAs. DeTiN then applies the joint TiN estimate to reclassify SSNVs and InDels as somatic or germline. 

## Code Example

python deTiN.py --mutation_data_path example_data/HCC_10_90.call_stats.pon_removed.txt --cn_data_path example_data/HCC-1143_100_T-sim-final.acs.seg --tumor_het_data example_data/HCC_10_90.tumor.hets.tsv --normal_het_data example_data/HCC_10_90.normal.hets.tsv --exac_data_path example_data/exac.pickle_high_af --output_name 10_percent_TiN_simulation --indel_data_path example_data/MuTect2.call_stats.txt --indel_data_type MuTect2 --output_dir example_data/
## Motivation

Genomic characterization is vital to the understanding and treatment of cancer.  Detection of somatic mutations is a critical component of this process. A key step in sensitive and specific somatic mutation detection is comparison of the tumor sample to a matched germline control. Sensitivity to detect somatic variants is greatly reduced when the matched normal sample is contaminated with tumor cells. To overcome this limitation, we developed deTiN, a method that estimates tumor-in-normal contamination (TiN), and improves detection sensitivity when using a contaminated normal. 

## Installation

Scripts are standalone but require Numpy, Pandas, Scipy and Python 2.7. 

git clone https://github.com/broadinstitute/deTiN.git

## Example Data

The above code example will run using the included data from an artifically mixed 10% contaminated normal. Input files were generated using MuTect and GATK4ACNV. 

## License

Copyright 2017 Amaro Taylor-Weiner

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
