# tcgaBrca
The Johan Staaf-lab repository containing scripts for downloading and preprocessing TCGA BRCA data from the [Genomic Data Commons](https://portal.gdc.cancer.gov/) repository. 

## Description
The pipeline uses the included manifest and annotation files to perform automated download and processing of gene expression, DNA methylation, copy-number, and somatic mutation data for all available tumor samples. A core data set with available data on all levels is produced and filtered against TCGA sample blacklists. The scripts produce a final data set comprising 630 unique breast cancer cases. 

## Data types

### RNA-seq
Downloaded "as is". Avaliable as counts, fpkm and upper-quantile normalized data. 

### SNP6-Copy number data
Downloaded "as is" with segmentation data on gene- and genome level.

### DNA methylation
Downloads raw Illumina iDat-files which are subsequently normlaized using the [minfi](https://bioconductor.org/packages/release/bioc/html/minfi.html) package followed by Infinium I/II scaling. Poorly perfomring probes are filtered using a CpG blacklist by [Zhou et al.](http://zwdzwd.github.io/InfiniumAnnotation). Annotations for Illumina Infinium CpGs are generated using custom scripts.

### Whole exome sequencing
Downloaded "as is" for the MUSE, SomaticSniper, Mutect, and varscan pipelines. Mutations called using at least two independent methods were kept in a filtering step, resulting in 58 973 SNVs across 630 samples.

## Usage
The scripts contained perform data download, preprocessing and annotation.

