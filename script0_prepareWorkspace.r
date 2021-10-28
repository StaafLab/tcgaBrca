######================================================######
### script 0 - prepare github workspace                  ###
######================================================######

##1. create directory structure for TCGA/GDC download manifests
##2. add manually created manifests by cancer and data type
	#1. Harmonized hg38 data: https://portal.gdc.cancer.gov/
	#2. Legacy hg19 data: https://portal.gdc.cancer.gov/legacy-archive/search/f
##3. build script pipeline for download and processing of data

##The full workspace can be cloned from github using the command
'git clone https://github.com/StaafLab/tcgaBrca.git'

##no need to run commands below if you use the gitHub-link

############################################################
##Main folder
dir.create("tcgaBrca")
setwd("tcgaBrca")

############################################################
##create work folders in "."
dir.create("manifest/brca",recursive=T)

##subfolders
for (i in c("rnaseq",
            "me450",
            "snp6/ascat2","snp6/dnacopy",
            "wes",
            "image/dx","image/tissue",
            "rnaseq_hg38")) {
  dir.create(paste0("manifest/brca/",i),recursive=T)
} ; rm(i)

############################################################
##Pan-immune subtype annotations https://gdc.cancer.gov/about-data/publications/panimmune

##container for processed manifests
dir.create("manifest/brca/annotations",recursive=T)

download.file("https://api.gdc.cancer.gov/data/0f31b768-7f67-4fc4-abc3-06ac5bd90bf0","manifest/brca/annotations/TCGASubtype.20170308.tsv",mode="wb")

#dir.create("manifest/luad/annotations",recursive=T)
#dir.create("manifest/lusc/annotations",recursive=T)

############################################################
##TCGA 2018 Science paper - ATAC data manifest
dir.create("manifest/atac",recursive=T)

##get ATAC manifest
#download.file("https://gdc.cancer.gov/node/1060/","manifest/atac/dataManifest_tcgaAtac_2018.txt",mode="wb")
download.file("https://gdc.cancer.gov/files/public/file/ATACseq-AWG_Open_GDC-Manifest.txt","manifest/atac/ATACseq-AWG_Open_GDC-Manifest.txt",mode="wb")

############################################################
##TCGA 2018 Pan-can annotations
dir.create("manifest/pancan",recursive=T)

##get sample purity/ploidy from ABSOLUTE
download.file("http://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5","manifest/pancan/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",mode="wb")

##get recommended survival data file
download.file("https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81","manifest/pancan/TCGA-CDR-SupplementalTableS1.xlsx",mode="wb")

##get recommended clinical data file
download.file("http://api.gdc.cancer.gov/data/0fc78496-818b-4896-bd83-52db1f533c5c","manifest/pancan/clinical_PANCAN_patient_with_followup.tsv",mode="wb")

##get sample quality annotations
download.file("http://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf","manifest/pancan/merged_sample_quality_annotations.tsv",mode="wb"))

############################################################
##TCGA/GDC data transfer tool

dir.create("gdcTransferTool",recursive=T)

##prepare GDC downloader
#https://gdc.cancer.gov/files/public/file/gdc-client_v1.5.0_Windows_x64_0.zip
#https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.5.0_Ubuntu_x64.zip
#https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.5.0_OSX_x64.zip

##legacy
#https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.4.0_OSX_x64_10.12.6.zip
#https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.4.0_Windows_x64.zip
#https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.4.0_Ubuntu_x64.zip

##get GDC-client - UBUNTU
download.file("https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.5.0_Ubuntu_x64.zip","gdcTransferTool/gdc-client_v1.5.0_Ubuntu_x64.zip",mode="wb")
unzip("gdcTransferTool/gdc-client_v1.5.0_Ubuntu_x64.zip",exdir="gdcTransferTool")
unlink("gdcTransferTool/gdc-client_v1.5.0_Ubuntu_x64.zip")
#system('chmod 777 gdcTransferTool/gdc-client')

##get GDC-client - Mac
#download.file("https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.5.0_OSX_x64.zip","gdcTransferTool/gdc-client_v1.5.0_OSX_x64.zip",mode="wb")
#unzip("gdcTransferTool/gdc-client_v1.5.0_OSX_x64.zip",exdir="gdcTransferTool/tmp")
unzip("gdcTransferTool/gdc-client_v1.5.0_OSX_x64.zip",exdir="gdcTransferTool")
unlink("gdcTransferTool/gdc-client_v1.5.0_OSX_x64.zip")
#system('chmod 777 gdcTransferTool/gdc-client')

##get GDC-client - Win
download.file("https://gdc.cancer.gov/files/public/file/gdc-client_v1.5.0_Windows_x64_0.zip","gdcTransferTool/gdc-client_v1.5.0_Windows_x64.zip",mode="wb")
unzip("gdcTransferTool/gdc-client_v1.5.0_Windows_x64.zip",exdir="gdcTransferTool")
unlink("gdcTransferTool/gdc-client_v1.5.0_Windows_x64.zip")

q("no")
###END
