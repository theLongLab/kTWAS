
## Users’ Manual of kTWAS (Version 1.0.0)
### Preamble

The kTWAS program leverages TWAS-like feature selection (Elastic Net weights) followed by a SKAT-like kernel-based score test, to combine advantages from both approaches.

### Installation
kTWAS is a batteries-included JAR executable. All needed external jar packages are included in the downloadable, kTWAS.jar. To download all necessary files, users can use the command 
`git clone https://github.com/theLongLab/kTWAS.git`
As we used an R package SKAT, the users have to install R and SKAT (https://cran.r-project.org/web/packages/SKAT/index.html. The versions of R and R package L0Learn that we have used on our platform are: version 2.0.0 for SKAT and version 3.5.1 for R. Other versions are not tested, although they may work. Users are also expected to have java (version: 1.8) on their platform. plink (v1.07, https://zzz.bwh.harvard.edu/plink/download.shtml) should also be installed.

Usage:

`java -jar PoolHapX.jar kTWAS -format CSV|PLINK -input_genotype INPUT_GENOTYPE_FILE -input_phenotype INPUT_PHENOTYPE_FILE -input_phenotype_column INPUT_PHENOTYPE_COLUMN_START:2-en_info_path INPUT_ELASTICNET_INFORMATION_FILE -gene INPUT_ENSEMBL_GENE_ID -output_folder OUTPUT_FOLDER_PATH`


### Copyright License (MIT Open Source)
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
