#####======================================================================#####
###  Filter TCGA BRCA against blacklists
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

load(file=paste0(HOME,"/","coreData_gexCnWesMeAtac_unique_matched_samples.RData"))

ls()
#  [1] "data.atac.core"      "data.cn.core"        "data.counts.core"   
#  [4] "data.fpkm.core"      "data.mutations.core" "data.seg.core"      
#  [7] "data.uq.core"        "HOME"                "MANIFEST"           
# [10] "mutations.core"      "samples.me.core"     "TUMOR_TYPE"        

##load betaFinal,sampleMat and probesKeep
load(file=paste0(HOME,"/me/","workspace_minfiNormInfiniumAdjFinalBetas.RData"))

ls()
#  [1] "betaFinal"           "data.atac.core"      "data.cn.core"       
#  [4] "data.counts.core"    "data.fpkm.core"      "data.mutations.core"
#  [7] "data.seg.core"       "data.uq.core"        "HOME"               
# [10] "MANIFEST"            "mutations.core"      "probesKeep"         
# [13] "sampleMat"           "samples.me.core"     "TUMOR_TYPE"         

################################################################################
##Filter vs sample exclude list and results from beta adjustment

all(sampleMat$Sample_ID==colnames(data.counts.core))
#[1] TRUE

##load TCGA blacklist
blacklist<-read.delim(paste0(HOME,"/sampleAnno/merged_sample_quality_annotations.tsv"),header=TRUE,as.is=TRUE)

str(blacklist)
# 'data.frame':   79286 obs. of  12 variables:
#  $ patient_barcode                  : chr  "TCGA-01-0628" "TCGA-01-0628" "TCGA-01-0630" "TCGA-01-0630" ...
#  $ aliquot_barcode                  : chr  "TCGA-01-0628-11A-01D-0356-01" "TCGA-01-0628-11A-01D-0383-05" "TCGA-01-0630-11A-01D-0356-01" "TCGA-01-0630-11A-01D-0383-05" ...
#  $ cancer.type                      : chr  "OV" "OV" "OV" "OV" ...
#  $ platform                         : chr  "Genome_Wide_SNP_6" "HumanMethylation27" "Genome_Wide_SNP_6" "HumanMethylation27" ...
#  $ patient_annotation               : chr  "Organ-Specific Control" "Organ-Specific Control" "Organ-Specific Control" "Organ-Specific Control" ...
#  $ sample_annotation                : chr  "" "" "" "" ...
#  $ aliquot_annotation               : chr  "" "" "" "" ...
#  $ aliquot_annotation_updated       : chr  "" "" "" "" ...
#  $ AWG_excluded_because_of_pathology: num  0 0 0 0 0 0 0 0 0 0 ...
#  $ AWG_pathology_exclusion_reason   : chr  "" "" "" "" ...
#  $ Reviewed_by_EPC                  : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ Do_not_use                       : chr  "False" "False" "False" "False" ...

blacklist$sampleID<-sub("-...-....-..","",blacklist$aliquot_barcode)

table(blacklist$platform)
 #                            Genome_Wide_SNP_6     HumanMethylation27    HumanMethylation450    IlluminaGA_miRNASeq    IlluminaGA_RNASeqV2 
 #                  1013                  22632                   2734                   9854                   2713                    646 
 # IlluminaHiSeq_DNASeqC IlluminaHiSeq_miRNASeq IlluminaHiSeq_RNASeqV2          MDA_RPPA_Core 
 #                  2347                  18445                  10439                   8463 

table(unique(blacklist$sampleID) %in% sampleMat$Sample_ID)
# FALSE  TRUE 
# 40862   669 

table(sampleMat$Sample_ID %in% unique(blacklist$sampleID) )
# TRUE 
#  669 

blacklist<-blacklist[blacklist$sampleID %in% sampleMat$Sample_ID,]

table(blacklist$platform)
#                             Genome_Wide_SNP_6    HumanMethylation450    IlluminaGA_miRNASeq  IlluminaHiSeq_DNASeqC IlluminaHiSeq_miRNASeq 
#                      1                    671                    671                    125                     17                    544 
# IlluminaHiSeq_RNASeqV2 
#                    669 

blacklist<-blacklist[blacklist$platform %in% c("Genome_Wide_SNP_6","HumanMethylation450","IlluminaHiSeq_RNASeqV2"),]

table(blacklist$Do_not_use)
# False  True 
#  1996    15 

##for a few there are 2 or more samples annotated. Choose 1..
samplesMe<-blacklist[blacklist$platform %in% c("HumanMethylation450"),]

samplesGex<-blacklist[blacklist$platform %in% c("IlluminaHiSeq_RNASeqV2"),]

samplesSnp<-blacklist[blacklist$platform %in% c("Genome_Wide_SNP_6"),]

##only care about "Do not use"
samplesMe<-unlist(lapply(split(samplesMe$Do_not_use,samplesMe$sampleID),function(x) { if( any(x=="True")) { return(FALSE)  } else { return(TRUE) } } ))

samplesGex<-unlist(lapply(split(samplesGex$Do_not_use,samplesGex$sampleID),function(x) { if( any(x=="True")) { return(FALSE)  } else { return(TRUE) } } ))

samplesSnp<-unlist(lapply(split(samplesSnp$Do_not_use,samplesSnp$sampleID),function(x) { if( any(x=="True")) { return(FALSE)  } else { return(TRUE) } } ))

table(samplesMe,samplesGex)
#          samplesGex
# samplesMe FALSE TRUE
#     FALSE     5    0
#     TRUE      0  664
table(samplesMe,samplesSnp)
#          samplesSnp
# samplesMe FALSE TRUE
#     FALSE     5    0
#     TRUE      0  664
table(samplesGex,samplesSnp)
#           samplesSnp
# samplesGex FALSE TRUE
#      FALSE     5    0
#      TRUE      0  664

identical(names(samplesMe),names(samplesGex))
#[1] TRUE

identical(names(samplesMe),names(samplesSnp))
#[1] TRUE

samplesUse<-intersect( names(samplesMe)[samplesMe] , intersect( names(samplesGex)[samplesGex],names(samplesSnp)[samplesSnp] ) )

length(samplesUse)
#[1] 664

##664 potentially ok to use

################################################################################
### Name swap

##clean Methylation
identical(paste0(sampleMat$Slide,"_",sampleMat$Array),colnames(betaFinal))
#[1] TRUE

##add rownames to smaple anno obj
rownames(sampleMat)<-paste0(sampleMat$Slide,"_",sampleMat$Array)

##swap names to TCGA-id
rownames(sampleMat)<-sampleMat$Sample_ID
colnames(betaFinal)<-sampleMat$Sample_ID

identical(rownames(sampleMat),colnames(betaFinal))
#[1] TRUE
identical(colnames(betaFinal),colnames(data.fpkm.core))
#[1] TRUE

################################################################################
### Additional sample removal?

if( file.exists( paste0(HOME,"/me/","methylationSampleBlacklist.txt") ) ) {

  ##load methylation normalization blacklist
  blacklist.me<-read.table(paste0(HOME,"/me/","methylationSampleBlacklist.txt"),sep="\t",as.is=T)

  ##Infinium intensity filter and Ringner-correction failures could be filtered
    #infinium intensity filter seldom used - no info on impact on data
    #Ringner (infinium I/II) adjustment failures can become quite skewed

  flaggedIDs<-NULL

  ##Filter Ringner
  flaggedIDs<-sub(".+\\t","",readLines(con=paste0(HOME,"/me/","methylationSampleBlacklist.txt")))

  print(flaggedIDs)
  #  [1] "6042324072_R01C02" "9283265077_R02C02" "9283265014_R01C01"
  #  [4] "9283265014_R01C02" "9993943005_R01C01" "6004791006_R06C02"
  #  [7] "9283265014_R06C02" "8784241047_R01C01" "6042324070_R03C02"
  # [10] "6042324070_R02C02" "6042316130_R02C01" "6026818055_R04C02"
  # [13] "9422492113_R01C01" "9993943017_R06C02"

  print(length(samplesUse))
  #[1] 664

  if( !is.null(flaggedIDs) & any(flaggedIDs %in% paste0(sampleMat$Slide,"_",sampleMat$Array)) ) {
    flaggedIDs<-rownames(sampleMat)[ which( paste0(sampleMat$Slide,"_",sampleMat$Array) %in% flaggedIDs ) ]
    if( any(flaggedIDs %in% samplesUse)  ) {
      samplesUse<-samplesUse[ -which(samplesUse %in% flaggedIDs) ]
    }
  }

  print(length(samplesUse))
  #[1] 650

  rm(flaggedIDs)

}

################################################################################
### Prepare beta adjustment object

## Get tumor purities

##from TCGA pancan companion page
sample.purity<-read.table(paste0(MANIFEST,"/pancan/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"),sep="\t",header=TRUE,as.is=TRUE)

str(sample.purity)
# 'data.frame':   10786 obs. of  10 variables:
#  $ array                    : chr  "TCGA-OR-A5J1-01" "TCGA-OR-A5J2-01" "TCGA-OR-A5J3-01" "TCGA-OR-A5J4-01" ...
#  $ sample                   : chr  "TCGA-OR-A5J1-01A-11D-A29H-01" "TCGA-OR-A5J2-01A-11D-A29H-01" "TCGA-OR-A5J3-01A-11D-A29H-01" "TCGA-OR-A5J4-01A-11D-A29H-01" ...
#  $ call.status              : chr  "called" "called" "called" "called" ...
#  $ purity                   : num  0.9 0.89 0.93 0.87 0.93 0.69 0.84 0.76 0.84 0.75 ...
#  $ ploidy                   : num  2 1.3 1.27 2.6 2.79 3.34 2.6 1.23 2.61 5.52 ...
#  $ Genome.doublings         : num  0 0 0 1 1 1 1 0 1 2 ...
#  $ Coverage.for.80..power   : num  9 6 5 12 12 17 12 7 12 32 ...
#  $ Cancer.DNA.fraction      : num  0.9 0.84 0.89 0.89 0.95 0.79 0.87 0.66 0.87 0.89 ...
#  $ Subclonal.genome.fraction: num  0.02 0.16 0.11 0.08 0.15 0.06 0 0.52 0.23 0.12 ...
#  $ solution                 : chr  "new" "new" "new" "new" ...

##use cancer DNA fraction as this is what should correspond best to what is seen on betas
str(samples.me.core)
# 'data.frame':   669 obs. of  2 variables:
#  $ tcga_id: chr  "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
#  $ array  : chr  "9993943013_R04C01" "9993943013_R01C02" "9993943005_R02C02" "9993943017_R01C01" ...

length(unique(sample.purity$array)) == length((sample.purity$array))
#[1] TRUE

length(intersect(sample.purity$array,sub(".$","",samplesUse)))
#[1] 645

uids<-intersect(sample.purity$array,sub(".$","",samplesUse))

all(sample.purity$array[match(uids,sample.purity$array)]==uids)
#[1] TRUE

sample.purity<-sample.purity[match(uids,sample.purity$array),]

head(sample.purity)
#               array                       sample call.status purity ploidy Genome.doublings
# 492 TCGA-3C-AAAU-01 TCGA-3C-AAAU-01A-11D-A41E-01      called   0.81   2.56                1
# 493 TCGA-3C-AALI-01 TCGA-3C-AALI-01A-11D-A41E-01      called   0.64   3.05                1
# 494 TCGA-3C-AALJ-01 TCGA-3C-AALJ-01A-31D-A41E-01      called   0.53   3.90                1
# 495 TCGA-3C-AALK-01 TCGA-3C-AALK-01A-11D-A41E-01      called   0.61   2.05                0
# 496 TCGA-4H-AAAK-01 TCGA-4H-AAAK-01A-12D-A41E-01      called   0.47   3.99                1
# 497 TCGA-5L-AAT0-01 TCGA-5L-AAT0-01A-12D-A41E-01      called   0.53   2.00                0
#     Coverage.for.80..power Cancer.DNA.fraction Subclonal.genome.fraction solution
# 492                     12                0.85                      0.05      new
# 493                     17                0.73                      0.00      new
# 494                     28                0.69                      0.13      new
# 495                     13                0.62                      0.03      new
# 496                     32                0.64                      0.00      new
# 497                     15                0.53                      0.00      new

table( colSums(apply(sample.purity,1,is.na)) )
#   0   1   6 
# 630   2  13 

sample.purity[colSums(apply(sample.purity[,c("purity","ploidy","Cancer.DNA.fraction","Subclonal.genome.fraction")],1,is.na))>0,]
#                 array                               sample call.status purity ploidy Genome.doublings
# 9859  TCGA-A1-A0SO-01         TCGA-A1-A0SO-01A-22D-A087-01                 NA     NA               NA
# 9860  TCGA-A2-A0EP-01         TCGA-A2-A0EP-01A-52D-A22W-01                 NA     NA               NA
# 9861  TCGA-A2-A0ST-01         TCGA-A2-A0ST-01A-12D-A087-01                 NA     NA               NA
# 9862  TCGA-A2-A1G6-01         TCGA-A2-A1G6-01A-11D-A13J-01                 NA     NA               NA
# 9863  TCGA-A2-A25B-01         TCGA-A2-A25B-01A-11D-A166-01                 NA     NA               NA
# 9864  TCGA-A2-A25F-01         TCGA-A2-A25F-01A-11D-A166-01                 NA     NA               NA
# 9865  TCGA-A2-A3XW-01         TCGA-A2-A3XW-01A-11D-A238-01                 NA     NA               NA
# 9869  TCGA-AO-A0JC-01         TCGA-AO-A0JC-01A-11D-A059-01                 NA     NA               NA
# 9872  TCGA-BH-A1FR-01         TCGA-BH-A1FR-01A-11D-A13J-01                 NA     NA               NA
# 9873  TCGA-BH-A42U-01         TCGA-BH-A42U-01A-12D-A242-01                 NA     NA               NA
# 9876  TCGA-E2-A1B6-01         TCGA-E2-A1B6-01A-31D-A12N-01                 NA     NA               NA
# 9878  TCGA-LD-A9QF-01         TCGA-LD-A9QF-01A-32D-A41E-01                 NA     NA               NA
# 9879  TCGA-LL-A5YO-01         TCGA-LL-A5YO-01A-21D-A28A-01                 NA     NA               NA
# 10545 TCGA-B6-A0X4-01 WU-BRCA-TCGA-B6-A0X4-01A-11D-A10G-09    maf_call   0.81   2.03                0
# 10730 TCGA-E9-A1RD-01 WU-BRCA-TCGA-E9-A1RD-01A-11D-A159-09    maf_call   0.65   3.15                1
#       Coverage.for.80..power Cancer.DNA.fraction Subclonal.genome.fraction solution
# 9859                      NA                  NA                        NA      new
# 9860                      NA                  NA                        NA      new
# 9861                      NA                  NA                        NA      new
# 9862                      NA                  NA                        NA      new
# 9863                      NA                  NA                        NA      new
# 9864                      NA                  NA                        NA      new
# 9865                      NA                  NA                        NA      new
# 9869                      NA                  NA                        NA      new
# 9872                      NA                  NA                        NA      new
# 9873                      NA                  NA                        NA      new
# 9876                      NA                  NA                        NA      new
# 9878                      NA                  NA                        NA      new
# 9879                      NA                  NA                        NA      new
# 10545                     10                0.81                        NA      old
# 10730                     17                0.74                        NA      old

sample.purity<-sample.purity[colSums(apply(sample.purity[,c("purity","ploidy","Cancer.DNA.fraction","Subclonal.genome.fraction")],1,is.na))==0,]

dim(sample.purity)
#[1] 630  10

uids<-intersect(sample.purity$array,uids)
samplesUse<-samplesUse[which(sub(".$","",samplesUse) %in% sample.purity$array)]

sample.purity<-sample.purity[match(uids,sample.purity$array),]
samplesUse<-samplesUse[match(uids,sub(".$","",samplesUse))]

identical(sample.purity$array,uids)
#[1] TRUE
identical(sample.purity$array,sub(".$","",samplesUse))
#[1] TRUE

##get ascat tumor
fracA<-sample.purity$purity

##get battenberg tumor
fracB<-sample.purity$Cancer.DNA.fraction

quantile(fracA,na.rm=TRUE)
#     0%    25%    50%    75%   100% 
# 0.0900 0.4400 0.5700 0.7175 1.0000 

quantile(fracB,na.rm=TRUE)
#   0%  25%  50%  75% 100% 
# 0.14 0.50 0.63 0.75 1.00 

cor.test(fracA,fracB,na.rm=T)
#         Pearson's product-moment correlation

# data:  fracA and fracB
# t = 70.731, df = 628, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9331838 0.9507025
# sample estimates:
#       cor 
# 0.9425881 

length(uids)
#[1] 630

##Save 
write.table(cbind(sample=sample.purity$array,tumor_purity=sample.purity$purity),
  file=paste0(HOME,"/me/samplePurityVector.txt"),
  sep="\t",quote=F,row.names=FALSE)

write.table(cbind(sample=sample.purity$array,tumor_dna_fraction=sample.purity$Cancer.DNA.fraction),
  file=paste0(HOME,"/me/sampleTumorDnaFractionVector.txt"),
  sep="\t",quote=F,row.names=FALSE)

write.table(cbind(sample=sample.purity$array,tumor_ploidy=sample.purity$ploidy),
  file=paste0(HOME,"/me/samplePloidyVector.txt"),
  sep="\t",quote=F,row.names=FALSE)

write.table(cbind(sample=sample.purity$array,tumor_subclonal_fraction=sample.purity$Subclonal.genome.fraction),
  file=paste0(HOME,"/me/sampleSubclonalFractionVector.txt"),
  sep="\t",quote=F,row.names=FALSE)

#make vector
fracTum<-sample.purity$Cancer.DNA.fraction
rm(fracA,fracB)

length(fracTum)
#[1] 630

################################################################################
###Make object for flexmix correction

identical(sample.purity$array,uids)
#[1] TRUE
identical(sample.purity$array,sub(".$","",samplesUse))
#[1] TRUE

length( intersect(samplesUse,colnames(betaFinal)) )
#[1] 630

betaNew<-betaFinal[,samplesUse]

str(betaNew)
 # num [1:421368, 1:630] 0.718 0.766 0.189 0.925 0.868 0.642 0.95 1 0.986 1 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:421368] "cg21870274" "cg08258224" "cg16619049" "cg18147296" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

################################################################################
###clean workspace - only need beta object and purities for adjustment

##source correction function - flexmix loaded on source
source(paste0(sub("/manifest","",MANIFEST),"/function_correctBetas.r"))

save(list=ls()[  ls() %in% c("fracTum","betaNew","adjustBeta","uids","MANIFEST","HOME") ],
  file=paste0(HOME,"/me/object_flexmixAdjustmentObject.RData")
  )

ls()[  ls() %in% c("fracTum","betaNew","adjustBeta","uids","MANIFEST","HOME") ]
#[1] "adjustBeta" "betaNew"    "fracTum"    "HOME"       "MANIFEST"   "uids"      

gc()

rm(betaNew,adjustBeta,fracTum)

################################################################################
### Do layer intersects

##Methylation 

###Filter against blacklist
table(rownames(sampleMat) %in% samplesUse,colnames(betaFinal) %in% samplesUse)
  #       FALSE TRUE
  # FALSE    39    0
  # TRUE      0  630

##lose a few more..
table(rownames(sampleMat) %in% samplesUse )
# FALSE  TRUE 
#    39   630 

##check intersect
length(intersect(rownames(sampleMat),samplesUse))
#[1] 630

samplesUseFinal<-intersect(samplesUse,rownames(sampleMat))

##GEX
str(data.counts.core)
 # int [1:60483, 1:669] 0 0 839 0 320 4151 0 3 22081 1 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60483] "ENSG00000242268.2" "ENSG00000270112.3" "ENSG00000167578.15" "ENSG00000273842.1" ...
 #  ..$ : chr [1:669] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

str(data.fpkm.core)
 # num [1:60483, 1:669] 0 0 3.68 0 1.83 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60483] "ENSG00000242268.2" "ENSG00000270112.3" "ENSG00000167578.15" "ENSG00000273842.1" ...
 #  ..$ : chr [1:669] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
 
str(data.uq.core)
 # num [1:60483, 1:669] 0 0 76028 0 37765 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60483] "ENSG00000242268.2" "ENSG00000270112.3" "ENSG00000167578.15" "ENSG00000273842.1" ...
 #  ..$ : chr [1:669] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

length(intersect(colnames(data.counts.core),samplesUse))
#[1] 630

samplesUseFinal<-intersect(samplesUseFinal,colnames(data.counts.core))

##Mutations
str(mutations.core)
# 'data.frame':   63767 obs. of  18 variables:
#  $ Tumor_Sample_Barcode  : chr  "TCGA-3C-AAAU-01A" "TCGA-3C-AAAU-01A" "TCGA-3C-AAAU-01A" "TCGA-3C-AAAU-01A" ...
#  $ Hugo_Symbol           : chr  "MST1L" "DDAH1" "CELSR2" "PLEKHO1" ...
#  $ Entrez_Gene_Id        : chr  "11223" "23576" "1952" "51177" ...
#  $ Gene                  : chr  "ENSG00000186715" "ENSG00000153904" "ENSG00000143126" "ENSG00000023902" ...
#  $ Chromosome            : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ Start_Position        : chr  "16757604" "85324885" "109274629" "150155993" ...
#  $ End_Position          : chr  "16757604" "85324885" "109274629" "150155993" ...
#  $ Strand                : chr  "+" "+" "+" "+" ...
#  $ Reference_Allele      : chr  "G" "T" "C" "G" ...
#  $ Tumor_Seq_Allele1     : chr  "G" "T" "C" "G" ...
#  $ Tumor_Seq_Allele2     : chr  "A" "A" "G" "T" ...
#  $ IMPACT                : chr  "MODIFIER" "HIGH" "MODIFIER" "MODIFIER" ...
#  $ Variant_Classification: chr  "RNA" "Splice_Site" "3'UTR" "Intron" ...
#  $ Variant_Type          : chr  "SNP" "SNP" "SNP" "SNP" ...
#  $ muse_row              : int  NA 89503 89504 89506 89507 89508 89509 89510 89511 89512 ...
#  $ mutect_row            : int  106929 106931 106932 106934 106937 106938 106939 106940 106941 106943 ...
#  $ somsniper_row         : int  NA 33141 33142 33143 33144 NA 33145 33146 33147 NA ...
#  $ varscan_row           : int  52748 52749 52750 52751 52752 52753 NA 52754 52755 NA ...

length(intersect(mutations.core$Tumor_Sample_Barcode,samplesUse))
#[1] 630

samplesUseFinal<-intersect(samplesUseFinal,mutations.core$Tumor_Sample_Barcode)

##SNP
dim(data.cn.core)
#[1] 60446   669

length(unique(colnames(data.cn.core)))
#[1] 669

dim(data.seg.core)
#[1] 133229      7

length(unique(data.seg.core$Sample_ID))
#[1] 669

length(intersect(colnames(data.cn.core),samplesUse))
#[1] 630

length(intersect(data.seg.core$Sample_ID,samplesUse))
#[1] 630

samplesUseFinal<-intersect(samplesUseFinal,colnames(data.cn.core))
samplesUseFinal<-intersect(samplesUseFinal,data.seg.core$Sample_ID)

length(samplesUseFinal)
#[1] 630

table(colnames(data.atac.core) %in% samplesUseFinal)
# TRUE 
#   39 

################################################################################
### Filter to final overlaps

identical(rownames(sampleMat),colnames(betaFinal))
#[1] TRUE
identical(colnames(betaFinal),colnames(data.counts.core))
#[1] TRUE

##Meth
sampleMat<-sampleMat[samplesUseFinal,]
betaFinal<-betaFinal[,samplesUseFinal]

##Gex
data.counts.core<-data.counts.core[,samplesUseFinal]
data.fpkm.core<-data.fpkm.core[,samplesUseFinal]
data.uq.core<-data.uq.core[,samplesUseFinal]

##Mutations
mutations.core<-mutations.core[mutations.core$Tumor_Sample_Barcode %in% samplesUseFinal,]

##Snp6
data.cn.core<-data.cn.core[,samplesUseFinal]

data.seg.core<-data.seg.core[data.seg.core$Sample_ID %in% samplesUseFinal,]

identical(rownames(sampleMat),colnames(betaFinal))
#[1] TRUE

identical(colnames(data.counts.core),colnames(data.fpkm.core))
#[1] TRUE
identical(colnames(data.counts.core),colnames(data.uq.core))
#[1] TRUE

all(colnames(data.cn.core) %in% data.seg.core$Sample_ID)
#[1] TRUE
all(data.seg.core$Sample_ID %in% colnames(data.cn.core))
#[1] TRUE

identical(rownames(sampleMat),colnames(data.counts.core))
#[1] TRUE
identical(rownames(sampleMat),colnames(data.cn.core))
#[1] TRUE

##also fix ATAC
str(data.atac.core)
 # num [1:215920, 1:39] 1.031 1.667 1.839 1.375 0.769 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:215920] "chr1:17234-17733|BRCA_2" "chr1:180634-181133|BRCA_3" "chr1:181207-181706|BRCA_4" "chr1:183557-184056|BRCA_5" ...
 #  ..$ : chr [1:39] "TCGA-3C-AALJ-01A" "TCGA-4H-AAAK-01A" "TCGA-A2-A0SV-01A" "TCGA-A2-A0SW-01A" ...

data.atac.core<-data.atac.core[,intersect(colnames(data.atac.core),samplesUseFinal)]

################################################################################
##Get clinical TNBC from pancan annotations

tmp<-readLines(con=paste0(MANIFEST,"/pancan/","clinical_PANCAN_patient_with_followup.tsv"))

tmp<-lapply(tmp,function(x) unlist(strsplit(x,"\t")))

table( unlist(lapply(tmp,function(x) length(x))) )
#   745   746 
# 10945    12 

tmp<-lapply(tmp,function(x) { c(x,rep("",746-length(x))) })

table( unlist(lapply(tmp,function(x) length(x))) )
#   746 
# 10957 

tmpH<-tmp[[1]]
tmp<-tmp[-1]

tmp<-do.call("rbind",tmp)

length(unique(tmp[,2]))==nrow(tmp)
#[1] TRUE

colnames(tmp)<-tmpH
rm(tmpH)

rownames(tmp)<-tmp[,2]

table( sub("-...$","",samplesUseFinal) %in% rownames(tmp) )
# TRUE 
#  630 

tmp<-tmp[sub("-...$","",samplesUseFinal),]

str(tmp)
 # chr [1:630, 1:746] "6E7D5EC6-A469-467C-B748-237353C23416" "55262FCB-1B01-4480-B322-36570430C917" ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:630] "TCGA-3C-AAAU" "TCGA-3C-AALI" "TCGA-3C-AALJ" "TCGA-3C-AALK" ...
 #  ..$ : chr [1:746] "bcr_patient_uuid" "bcr_patient_barcode" "acronym" "gender" ...

table( sub("-...$","",samplesUseFinal) == rownames(tmp) )
# TRUE 
#  630 

##get receptor status when possible
colnames(tmp)[grep("estrogen",colnames(tmp))]
# [1] "pos_finding_metastatic_breast_carcinoma_estrogen_receptor_other_measuremenet_scale_text"
# [2] "breast_carcinoma_estrogen_receptor_status"                                              
# [3] "positive_finding_estrogen_receptor_other_measurement_scale_text"                        
# [4] "metastatic_breast_carcinoma_estrogen_receptor_detection_method_text"                    
# [5] "metastatic_breast_carcinoma_estrogen_receptor_level_cell_percent_category"              
# [6] "metastatic_breast_carcinoma_estrogen_receptor_status"          

colnames(tmp)[grep("progesterone",colnames(tmp),ignore.case=T)]
# [1] "breast_carcinoma_immunohistochemistry_progesterone_receptor_pos_finding_scale"         
# [2] "breast_carcinoma_progesterone_receptor_status"                                         
# [3] "pos_finding_progesterone_receptor_other_measurement_scale_text"                        
# [4] "metastatic_breast_carcinoma_progesterone_receptor_level_cell_percent_category"         
# [5] "metastatic_breast_carcinoma_pos_finding_progesterone_receptor_other_measure_scale_text"
# [6] "metastatic_breast_carcinoma_progesterone_receptor_status"                              
# [7] "progesterone_receptor_level_cell_percent_category"                                     
# [8] "metastatic_breast_carcinoma_progesterone_receptor_detection_method_text"               

colnames(tmp)[grep("her2|erbb",colnames(tmp),ignore.case=T)]
#  [1] "her2_erbb_method_calculation_method_text"                                                                    
#  [2] "her2_immunohistochemistry_level_result"                                                                      
#  [3] "her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text"                            
#  [4] "her2_erbb_pos_finding_cell_percent_category"                                                                 
#  [5] "her2_and_centromere_17_positive_finding_other_measurement_scale_text"                                        
#  [6] "metastatic_breast_carcinoma_lab_proc_her2_neu_immunohistochemistry_receptor_status"                          
#  [7] "metastatic_breast_carcinoma_her2_neu_chromosone_17_signal_ratio_value"                                       
#  [8] "metastatic_breast_carcinoma_her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text"
#  [9] "metastatic_breast_carcinoma_her2_erbb_pos_finding_cell_percent_category"                                     
# [10] "metastatic_breast_carcinoma_her2_erbb_method_calculation_method_text"                                        
# [11] "her2_neu_and_centromere_17_copy_number_analysis_input_total_number_count"                                    
# [12] "pos_finding_her2_erbb2_other_measurement_scale_text"                                                         
# [13] "metastatic_breast_carcinoma_erbb2_immunohistochemistry_level_result"                                         
# [14] "her2_neu_chromosone_17_signal_ratio_value"                                                                   
# [15] "metastatic_breast_carcinoma_lab_proc_her2_neu_in_situ_hybridization_outcome_type"                            
# [16] "metastatic_breast_carcinoma_pos_finding_her2_erbb2_other_measure_scale_text"                                 
# [17] "lab_proc_her2_neu_immunohistochemistry_receptor_status"                                                      
# [18] "lab_procedure_her2_neu_in_situ_hybrid_outcome_type"                                                          
# [19] "her2_neu_breast_carcinoma_copy_analysis_input_total_number"                                                  
# [20] "her2_neu_and_centromere_17_copy_number_metastatic_breast_carcinoma_analysis_input_total_number_count"        
# [21] "her2_neu_metastatic_breast_carcinoma_copy_analysis_input_total_number"                                       

getVars<-c(
  "bcr_patient_barcode",
  "diagnosis",
  "gender",
  "age_at_initial_pathologic_diagnosis",
  "race",
  "ethnicity",
  "vital_status",
  "days_to_birth",
  "days_to_death",
  "days_to_last_followup",
  "year_of_initial_pathologic_diagnosis",
  "radiation_therapy",
  "histological_type",
  "histological_type_other",
  "anatomic_neoplasm_subdivision",
  "pathologic_T",
  "pathologic_N",
  "pathologic_M",
  "pathologic_stage",
  "clinical_T",
  "clinical_N",
  "clinical_M",
  "clinical_stage",
  "number_of_lymphnodes_positive_by_ihc",
  "number_of_lymphnodes_positive_by_he",
  "postoperative_rx_tx",
  "menopause_status",
  "eastern_cancer_oncology_group",
  "karnofsky_performance_score",
  "history_of_neoadjuvant_treatment",
  ##tumor type specific
  "breast_carcinoma_progesterone_receptor_status",
  "breast_carcinoma_estrogen_receptor_status",
  "lab_proc_her2_neu_immunohistochemistry_receptor_status",
  "lab_procedure_her2_neu_in_situ_hybrid_outcome_type"
)

##save variables - all 
sampleAnno<-as.data.frame(tmp[,getVars],stringsAsFactors=FALSE)
rownames(sampleAnno)<-samplesUseFinal

str(sampleAnno)
# 'data.frame':   630 obs. of  34 variables:
#  $ bcr_patient_barcode                                   : chr  "TCGA-3C-AAAU" "TCGA-3C-AALI" "TCGA-3C-AALJ" "TCGA-3C-AALK" ...
#  $ diagnosis                                             : chr  "" "" "" "" ...
#  $ gender                                                : chr  "FEMALE" "FEMALE" "FEMALE" "FEMALE" ...
#  $ age_at_initial_pathologic_diagnosis                   : chr  "55" "50" "62" "52" ...
#  $ race                                                  : chr  "WHITE" "BLACK OR AFRICAN AMERICAN" "BLACK OR AFRICAN AMERICAN" "BLACK OR AFRICAN AMERICAN" ...
#  $ ethnicity                                             : chr  "NOT HISPANIC OR LATINO" "NOT HISPANIC OR LATINO" "NOT HISPANIC OR LATINO" "NOT HISPANIC OR LATINO" ...
#  $ vital_status                                          : chr  "Alive" "Alive" "Alive" "Alive" ...
#  $ days_to_birth                                         : chr  "-20211" "-18538" "-22848" "-19074" ...
#  $ days_to_death                                         : chr  "[Not Applicable]" "[Not Applicable]" "[Not Applicable]" "[Not Applicable]" ...
#  $ days_to_last_followup                                 : chr  "4047.0" "4005.0" "1474.0" "1448.0" ...
#  $ year_of_initial_pathologic_diagnosis                  : chr  "2004" "2003" "2011" "2011" ...
#  $ radiation_therapy                                     : chr  "NO" "YES" "NO" "NO" ...
#  $ histological_type                                     : chr  "Infiltrating Lobular Carcinoma" "Infiltrating Ductal Carcinoma" "Infiltrating Ductal Carcinoma" "Infiltrating Ductal Carcinoma" ...
#  $ histological_type_other                               : chr  "[Not Applicable]" "[Not Applicable]" "[Not Applicable]" "[Not Applicable]" ...
#  $ anatomic_neoplasm_subdivision                         : chr  "Left Lower Outer Quadrant" "Right Upper Outer Quadrant" "Right" "Right" ...
#  $ pathologic_T                                          : chr  "TX" "T2" "T2" "T1c" ...
#  $ pathologic_N                                          : chr  "NX" "N1a" "N1a" "N0 (i+)" ...
#  $ pathologic_M                                          : chr  "MX" "M0" "M0" "M0" ...
#  $ pathologic_stage                                      : chr  "Stage X" "Stage IIB" "Stage IIB" "Stage IA" ...
#  $ clinical_T                                            : chr  "[Not Available]" "[Not Available]" "[Not Available]" "[Not Available]" ...
#  $ clinical_N                                            : chr  "[Not Applicable]" "[Not Applicable]" "[Not Applicable]" "[Not Applicable]" ...
#  $ clinical_M                                            : chr  "[Not Applicable]" "[Not Applicable]" "[Not Applicable]" "[Not Applicable]" ...
#  $ clinical_stage                                        : chr  "[Not Applicable]" "[Not Applicable]" "[Not Applicable]" "[Not Applicable]" ...
#  $ number_of_lymphnodes_positive_by_ihc                  : chr  "[Not Available]" "[Not Available]" "[Not Available]" "1" ...
#  $ number_of_lymphnodes_positive_by_he                   : chr  "4" "1" "1" "0" ...
#  $ postoperative_rx_tx                                   : chr  "YES" "YES" "YES" "YES" ...
#  $ menopause_status                                      : chr  "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)" "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)" "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)" "[Unknown]" ...
#  $ eastern_cancer_oncology_group                         : chr  "" "" "" "" ...
#  $ karnofsky_performance_score                           : chr  "" "" "" "" ...
#  $ history_of_neoadjuvant_treatment                      : chr  "No" "No" "No" "No" ...
#  $ breast_carcinoma_progesterone_receptor_status         : chr  "Positive" "Positive" "Positive" "Positive" ...
#  $ breast_carcinoma_estrogen_receptor_status             : chr  "Positive" "Positive" "Positive" "Positive" ...
#  $ lab_proc_her2_neu_immunohistochemistry_receptor_status: chr  "Negative" "Positive" "Indeterminate" "Positive" ...
#  $ lab_procedure_her2_neu_in_situ_hybrid_outcome_type    : chr  "[Not Evaluated]" "[Not Evaluated]" "[Not Evaluated]" "[Not Evaluated]" ...
 
table(tmp[,"breast_carcinoma_progesterone_receptor_status"])
# [Not Evaluated]   Indeterminate        Negative        Positive 
#              40               2             186             402 

table(tmp[,"breast_carcinoma_estrogen_receptor_status"	])
# [Not Evaluated]        Negative        Positive 
#              39             134             457 

table(tmp[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"])
# [Not Available] [Not Evaluated]       Equivocal   Indeterminate        Negative        Positive 
#               3             117             109              11             317              73 

table(tmp[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"],
  tmp[,"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"]
  )
  #                 [Not Available] [Not Evaluated] Equivocal Indeterminate Negative Positive
  # [Not Available]               1               0         0             0        2        0
  # [Not Evaluated]              13              65         0             0       34        5
  # Equivocal                     2              10         2             0       80       15
  # Indeterminate                 0               3         0             0        6        2
  # Negative                     15             261         0             3       37        1
  # Positive                      5              55         0             0        6        7

##make hierarchy
her2<- rep("",nrow(tmp))
her2[ tmp[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"] == "[Not Evaluated]" | 
  tmp[,"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"] == "[Not Evaluated]" ] <- "NA"
her2[ tmp[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"] == "[Not Available]" | 
  tmp[,"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"] == "[Not Available]" ] <- "NA"
her2[ tmp[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"] == "Indeterminate" | 
  tmp[,"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"] == "Indeterminate" ] <- "Indeterminate"
her2[ tmp[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"] == "Equivocal" | 
  tmp[,"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"] == "Equivocal" ] <- "Equivocal"
her2[ tmp[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"] == "Negative" | 
  tmp[,"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"] == "Negative" ] <- "Negative"
her2[ tmp[,"lab_proc_her2_neu_immunohistochemistry_receptor_status"] == "Positive" | 
  tmp[,"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"] == "Positive" ] <- "Positive"

table(her2)
# her2
#     Equivocal Indeterminate            NA      Negative      Positive 
#            14             3            79           438            96 

combined<-paste0(paste("ER:",tmp[,"breast_carcinoma_estrogen_receptor_status" ]),
	paste(" PR:",tmp[,"breast_carcinoma_progesterone_receptor_status"]),
	paste(" HER2:",her2)
  )
names(combined)<-rownames(tmp)

length(combined)
#[1] 630

sampleAnno$receptorCombined<-combined

rm(combined)

##save variables - all 
sampleAnno$ER=tmp[,"breast_carcinoma_estrogen_receptor_status"]
sampleAnno$PR=tmp[,"breast_carcinoma_progesterone_receptor_status"]
sampleAnno$HER2=her2

sampleAnno$TNBC<-sampleAnno$ER=="Negative" &
	sampleAnno$PR=="Negative" &
	sampleAnno$HER2=="Negative"

head(sampleAnno[,-1:-34])
#                                               receptorCombined       ER       PR          HER2  TNBC
# TCGA-3C-AAAU-01A      ER: Positive PR: Positive HER2: Negative Positive Positive      Negative FALSE
# TCGA-3C-AALI-01A      ER: Positive PR: Positive HER2: Positive Positive Positive      Positive FALSE
# TCGA-3C-AALJ-01A ER: Positive PR: Positive HER2: Indeterminate Positive Positive Indeterminate FALSE
# TCGA-3C-AALK-01A      ER: Positive PR: Positive HER2: Positive Positive Positive      Positive FALSE
# TCGA-4H-AAAK-01A     ER: Positive PR: Positive HER2: Equivocal Positive Positive     Equivocal FALSE
# TCGA-5L-AAT0-01A      ER: Positive PR: Positive HER2: Negative Positive Positive      Negative FALSE

sampleAnno$erPositive<-sampleAnno$ER=="Positive"
sampleAnno$prPositive<-sampleAnno$PR=="Positive"
sampleAnno$her2Positive<-sampleAnno$HER2=="Positive"

##load pan-immune paper PAM50 subtypes (from https://gdc.cancer.gov/about-data/publications/panimmune)
tmp<-read.table(file=paste0(MANIFEST,"/brca/annotations/","TCGASubtype.20170308.tsv"),sep="\t",header=T,row.names=1,as.is=T)
tmp<-tmp[tmp$cancer.type=="BRCA",]

table( sub("TCGA-..-....-(..).-.+","\\1",rownames(tmp)) )
#   01   06   11 
# 1097    7  114 

length( unique(  sub("(TCGA-..-....)-...-.+","\\1",rownames(tmp))  ) )
#[1] 1097

length( unique(  sub("(TCGA-..-....)-...-.+","\\1",rownames(tmp))[sub("TCGA-..-....-(..).-.+","\\1",rownames(tmp)) == "01"]  ) )
#[1] 1097

tmp<-tmp[sub("TCGA-..-....-(..).-.+","\\1",rownames(tmp)) =="01",]
rownames(tmp)<-sub("(TCGA-..-....-...)-.+","\\1",rownames(tmp))

all(rownames(sampleAnno) %in% rownames(tmp))
#[1] TRUE

tmp<-tmp[ match(rownames(sampleAnno),rownames(tmp)),]

all(rownames(sampleAnno) == rownames(tmp))
#[1] TRUE

sampleAnno$subtype.canonical<-tmp$Subtype_Selected
sampleAnno$subtype.full<-tmp$Subtype_mRNA

identical(rownames(sampleAnno),colnames(betaFinal))
#[1] TRUE

table(sampleAnno[colnames(data.atac.core),"TNBC"],sampleAnno[colnames(data.atac.core),"pam50.full"])
  #       Basal Her2 LumA LumB Normal
  # FALSE     5    1   17   10      3
  # TRUE      3    0    0    0      0

rm(tmp,her2,getVars)

###Add purity info
str(sample.purity)
# 'data.frame':   630 obs. of  10 variables:
#  $ array                    : chr  "TCGA-3C-AAAU-01" "TCGA-3C-AALI-01" "TCGA-3C-AALJ-01" "TCGA-3C-AALK-01" ...
#  $ sample                   : chr  "TCGA-3C-AAAU-01A-11D-A41E-01" "TCGA-3C-AALI-01A-11D-A41E-01" "TCGA-3C-AALJ-01A-31D-A41E-01" "TCGA-3C-AALK-01A-11D-A41E-01" ...
#  $ call.status              : chr  "called" "called" "called" "called" ...
#  $ purity                   : num  0.81 0.64 0.53 0.61 0.47 0.53 0.83 0.47 0.44 0.64 ...
#  $ ploidy                   : num  2.56 3.05 3.9 2.05 3.99 2 3.12 2.04 3.26 2.1 ...
#  $ Genome.doublings         : num  1 1 1 0 1 0 1 0 1 0 ...
#  $ Coverage.for.80..power   : num  12 17 28 13 32 15 14 18 29 13 ...
#  $ Cancer.DNA.fraction      : num  0.85 0.73 0.69 0.62 0.64 0.53 0.89 0.47 0.56 0.65 ...
#  $ Subclonal.genome.fraction: num  0.05 0 0.13 0.03 0 0 0.1 0.11 0.02 0.05 ...
#  $ solution                 : chr  "new" "new" "new" "new" ...

identical(sampleAnno$bcr_patient_barcode,sub("-..$","",sample.purity$array))
#[1] TRUE

sampleAnno$purity<-sample.purity$purity
sampleAnno$ploidy<-sample.purity$ploidy
sampleAnno$Cancer.DNA.fraction<-sample.purity$Cancer.DNA.fraction
sampleAnno$Subclonal.genome.fraction<-sample.purity$Subclonal.genome.fraction

################################################################################
### Gather objects

##change matrix name to indicate that it is uncorrected betas.
betaOrig<-betaFinal

##change some other names
dataAtac<-data.atac.core

dataCn<-data.cn.core
dataSeg<-data.seg.core

dataMut<-mutations.core

gexCounts<-data.counts.core
gexFpkm<-data.fpkm.core
gexUq<-data.uq.core

save(betaOrig,
	dataAtac,
	dataCn,dataSeg,
	dataMut,
	gexCounts,gexFpkm,gexUq,
	sampleAnno,
	HOME,MANIFEST,TUMOR_TYPE,
	file=paste0(HOME,"/","workspace_blacklistFiltered_atacCnGexMeWes_withSampleAnnotations.RData")
)

q("no")
###END