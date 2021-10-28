#####======================================================================#####
### Create final workspace
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncoloy and Pathology

################################################################################
##Set home directory

##set/create own home directory below:

##work
HOME<-"~/hdd1/tcgaBrca"
MANIFEST<-"~/Documents/tcgaBrca/manifest"
##home
HOME<-"I:/data/tcgaBrca"
MANIFEST<-"F:/gitProjects/tcgaBrca/manifest"

##tumor type
TUMOR_TYPE<-"brca"

list.files(MANIFEST,full.names=T)
#[1] "/home/med-mai/Documents/tcgaBrca/manifest/atac"  
#[2] "/home/med-mai/Documents/tcgaBrca/manifest/brca"  
#[5] "/home/med-mai/Documents/tcgaBrca/manifest/pancan"

##create data directories
#dir.create(paste0(HOME,"/","me/norm"),recursive=TRUE)

################################################################################
##load required packages

if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer") }

library("RColorBrewer")

if(!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap") }

library("pheatmap")

if(!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr") }

library("dplyr")

if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges") }

library("GenomicRanges")

################################################################################
##Get core tumor set from methylation data

load(file=paste0(HOME,"/","workspace_blacklistFiltered_atacCnGexMeWes_withSampleAnnotations.RData"))

ls()
#  [1] "betaOrig"   "dataAtac"   "dataCn"     "dataMut"    "dataSeg"   
#  [6] "gexCounts"  "gexFpkm"    "gexUq"      "HOME"       "MANIFEST"  
# [11] "sampleAnno" "TUMOR_TYPE"

#object - "betaAdj" - purity adjusted betas
load(list.files(paste0(HOME,"/me/"),pattern="data450k_.+_minfiNormalized_ringnerAdjusted_purityAdjusted_tumorBetaValues.RData",full.names=TRUE) )

#object - "betaNorm" - inferred normal methylation
load(list.files(paste0(HOME,"/me/"),pattern="data450k_.+_minfiNormalized_ringnerAdjusted_purityAdjusted_normalBetaValues.RData",full.names=TRUE) )

#object - "annoObj" - custom CpG annotations
tmp<-list.files(paste0(HOME,"/annotateFeatures/"),pattern="object_noChromFilter_.+_expandedAnnotations.RData",full.names=TRUE)
tmp<-tmp[ length(tmp) ]
load( tmp  ) ; rm(tmp)

ls()

################################################################################
### Gather final object

ls()

save(annoObj,betaAdj,betaNorm,betaOrig,
	dataAtac,
	dataCn,dataSeg,
	dataMut,
	gexCounts,gexFpkm,gexUq,
	sampleAnno,
	HOME,MANIFEST,TUMOR_TYPE,
	file=paste0(HOME,"/","finalWorkspace_atacCnGexMeWes_withAnnotations.RData")
)

q("no")
###END