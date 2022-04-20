#####======================================================================#####
### Normalize methylation for TCGA BRCA
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncoloy and Pathology

################################################################################
##Set home directory

##downloads will require substantial amount of space, e.g. 1-2Tb for TCGA-BRCA

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

if(!requireNamespace("tools", quietly = TRUE)) {
  install.packages("tools") }

library("tools")

if(!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr") }

library("dplyr")

if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") }

library("BiocManager")

##swamp+dependencies
if (!requireNamespace("impute", quietly = TRUE)) {
    BiocManager::install("impute") }

library("impute")

if(!requireNamespace("swamp", quietly = TRUE)) {
  install.packages("swamp") }

library("swamp")

##granges
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    BiocManager::install("GenomicRanges") }

library("GenomicRanges")

##minfi+dependencies
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
    BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19") }

library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

if (!requireNamespace("IlluminaHumanMethylation450kmanifest", quietly = TRUE)) {
    BiocManager::install("IlluminaHumanMethylation450kmanifest") }

library("IlluminaHumanMethylation450kmanifest")

if (!requireNamespace("minfi", quietly = TRUE)) {
    BiocManager::install("minfi") }

library(minfi)

##850k
# if (!requireNamespace("IlluminaHumanMethylationEPICmanifest", quietly = TRUE)) {
#     BiocManager::install("IlluminaHumanMethylationEPICmanifest") }

# library("IlluminaHumanMethylationEPICmanifest")

# if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b2.hg19", quietly = TRUE)) {
#     BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19") }

# library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

################################################################################
##Get core tumor set from methylation data

load(file=paste0(HOME,"/","coreData_gexCnWesMeAtac_unique_matched_samples.RData"))

ls()
#  [1] "data.atac.core"      "data.cn.core"        "data.counts.core"   
#  [4] "data.fpkm.core"      "data.mutations.core" "data.seg.core"      
#  [7] "data.uq.core"        "HOME"                "MANIFEST"           
# [10] "mutations.core"      "samples.me.core"     "TUMOR_TYPE"        

str(data.counts.core)
 # int [1:60483, 1:669] 0 0 839 0 320 4151 0 3 22081 1 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60483] "ENSG00000242268.2" "ENSG00000270112.3" "ENSG00000167578.15" "ENSG00000273842.1" ...
 #  ..$ : chr [1:669] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

str(data.mutations.core)
 # int [1:60483, 1:669] 0 0 0 0 0 0 0 0 0 0 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60483] "ENSG00000242268" "ENSG00000270112" "ENSG00000167578" "ENSG00000273842" ...
 #  ..$ : chr [1:669] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

length(data.cn.core)
#[1] 669

data.cn.core[1:2,1:5]
#                   TCGA-3C-AAAU-01A TCGA-3C-AALI-01A TCGA-3C-AALJ-01A
# ENSG00000223972.5               NA               NA               NA
# ENSG00000227232.5               NA               NA               NA
#                   TCGA-3C-AALK-01A TCGA-4H-AAAK-01A
# ENSG00000223972.5               NA               NA
# ENSG00000227232.5               NA               NA

str(data.seg.core)
# 'data.frame': 133229 obs. of  7 variables:
#  $ Sample_ID        : chr  "TCGA-3C-AAAU-01A" "TCGA-3C-AAAU-01A" "TCGA-3C-AAAU-01A" "TCGA-3C-AAAU-01A" ...
#  $ Chromosome       : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ Start            : int  62920 72302815 72346221 85151891 85426191 149927059 150522687 151201646 151849527 153290365 ...
#  $ End              : int  72302735 72345465 85150197 85426026 149918981 150520671 151198313 151847886 153282939 153983915 ...
#  $ Copy_Number      : int  3 2 3 5 3 18 11 3 4 3 ...
#  $ Major_Copy_Number: int  2 2 2 3 2 17 10 2 3 2 ...
#  $ Minor_Copy_Number: int  1 0 1 2 1 1 1 1 1 1 ...

str(samples.me.core)
# 'data.frame': 669 obs. of  2 variables:
#  $ tcga_id: chr  "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
#  $ array  : chr  "9993943013_R04C01" "9993943013_R01C02" "9993943005_R02C02" "9993943017_R01C01" ...

##load raw methylation data
load(file=paste0(HOME,"/me/","me_minfi_RGset_unique_samples.RData"))

methSet
# class: RGChannelSet 
# dim: 622399 893 
# metadata(0):
# assays(2): Green Red
# rownames(622399): 10600313 10600322 ... 74810490 74810492
# rowData names(0):
# colnames(893): 9993943013_R04C01 9993943013_R01C02 ...
#   9993943017_R06C01 3999997079_R01C02
# colData names(6): Sample_ID Sample_Type ... Basename filenames
# Annotation
#   array: IlluminaHumanMethylation450k
#   annotation: ilmn12.hg19

str(samples.me.core)
# 'data.frame':	669 obs. of  2 variables:
#  $ tcga_id: chr  "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
#  $ array  : chr  "9993943013_R04C01" "9993943013_R01C02" "9993943005_R02C02" "9993943017_R01C01" ...

methSet<-methSet[,samples.me.core$array]

methSet
# class: RGChannelSet 
# dim: 622399 669 
# metadata(0):
# assays(2): Green Red
# rownames(622399): 10600313 10600322 ... 74810490 74810492
# rowData names(0):
# colnames(669): 9993943013_R04C01 9993943013_R01C02 ...
#   9993943017_R06C01 3999997079_R01C02
# colData names(6): Sample_ID Sample_Type ... Basename filenames
# Annotation
#   array: IlluminaHumanMethylation450k
#   annotation: ilmn12.hg19

################################################################################
##Do raw preprocess..

mSet <- preprocessRaw(methSet)

mSet
# class: MethylSet 
# dim: 485512 669 
# metadata(0):
# assays(2): Meth Unmeth
# rownames(485512): cg00050873 cg00212031 ... ch.22.47579720R
#   ch.22.48274842R
# rowData names(0):
# colnames(669): 9993943013_R04C01 9993943013_R01C02 ...
#   9993943017_R06C01 3999997079_R01C02
# colData names(6): Sample_ID Sample_Type ... Basename filenames
# Annotation
#   array: IlluminaHumanMethylation450k
#   annotation: ilmn12.hg19
# Preprocessing
#   Method: Raw (no normalization or bg correction)
#   minfi version: 1.34.0
#   Manifest version: 0.4.0

##check median array intensities
qc <- getQC(mSet)

plotQC(qc)
savePlot(filename=paste0(HOME,"/me/","minfiIntQCplot.tif"),type="tif")

dev.off()

save(qc,file=paste0(HOME,"/me/","minfi_mSetIntensityQc.RData"))

##poor quality sample/s to remove before normalization
  ##define which samples to remove
all(rownames(qc)==colnames(methSet))
#[1] TRUE
all(rownames(qc)==samples.me.core$array)
#[1] TRUE

##cutoff that separates good from bad samples is 
which( (qc$mMed + qc$uMed)/2 <10.5 )
#[1]  79 635 636 637 641

##Record bad samples, remove after normalization - will also have to filter against TCGA blacklist..
write(apply(samples.me.core[which( (qc$mMed + qc$uMed)/2 <10.5 ),],1,paste,collapse="\t"),file=paste0(HOME,"/me/","methylationSampleBlacklist.txt"),append=T)

rm(qc,mSet)

ls()
#  [1] "data.atac.core"      "data.cn.core"        "data.counts.core"   
#  [4] "data.fpkm.core"      "data.mutations.core" "data.seg.core"      
#  [7] "data.uq.core"        "HOME"                "MANIFEST"           
# [10] "methSet"             "mutations.core"      "samples.me.core"    
# [13] "TUMOR_TYPE"      

################################################################################
##check reported vs data derived sex for patients

GRset <- mapToGenome(methSet)
predictedSex <- getSex(GRset, cutoff = -2)
GRset <-addSex(GRset, sex = predictedSex)

pdf(paste0(HOME,"/me/","sexPredictionFromIntensitiesPlot.pdf"),width=12,height=12,useDingbats=F)
par(mfrow=c(3,1),mar=c(4.1,4.1,0.5,2.1),font=2,font.axis=2,font.lab=2)
plotSex(GRset) 
plot(predictedSex$xMed,predictedSex$yMed,col=1,pch=16,xlab="chrX int",ylab="chrY int",cex=1.2,ylim=c(6,14))
legend("bottomleft",legend=c("raw ints"),col=1:2,pch=16,bty="n")
abline(a=0,b=8/10)  
dev.off()

table(predictedSex$predictedSex)
#   F   M 
# 661   8 

all(samples.me.core$array==rownames(predictedSex))
#[1] TRUE

##Record male samples, remove after normalization - will also have to filter against TCGA blacklist later..
write(apply(samples.me.core[predictedSex$predictedSex=="M",],1,paste,collapse="\t"),file=paste0(HOME,"/me/","methylationSampleBlacklist.txt"),append=T)

##verified in patient table that 8 samples are male

rm(GRset)
gc()

################################################################################

##perform functional normalization as in Fortin et al. 2014 (Gen Biol)
  ##default parameters
    ##noob (bg+dye) correction + functional normalization
    ##pretty RAM-intensive, do not attempt with less than 32GB
 
all.equal(rownames(predictedSex),colnames(methSet))
#[1] TRUE

##normalize with predicted sex..
rgSetFnorm <- preprocessFunnorm(methSet,sex=predictedSex$predictedSex,keepCN=FALSE)
#[preprocessFunnorm] Background and dye bias correction with noob
#[preprocessFunnorm] Mapping to genome
#[preprocessFunnorm] Quantile extraction
#[preprocessFunnorm] Normalization

rgSetFnorm
# class: GenomicRatioSet 
# dim: 485512 669 
# metadata(0):
# assays(2): Beta CN
# rownames(485512): cg13869341 cg14008030 ... cg08265308 cg14273923
# rowData names(0):
# colnames(669): 9993943013_R04C01 9993943013_R01C02 ...
#   9993943017_R06C01 3999997079_R01C02
# colData names(6): Sample_ID Sample_Type ... Basename filenames
# Annotation
#   array: IlluminaHumanMethylation450k
#   annotation: ilmn12.hg19
# Preprocessing
#   Method: NA
#   minfi version: NA
#   Manifest version: NA

save(rgSetFnorm,file=paste0(HOME,"/me/","object_minfi_rgSetFnorm.RData"))

save(predictedSex,file=paste0(HOME,"/me/","object_minfi_predictedSex.RData"))

rm(predictedSex)

gc()

################################################################################
##Remove unreliable probes

##add annotations for 450K from Zhou et al. 2016 NAR
  ##https://doi.org/10.1093/nar/gkw967
    ##this publication may eventually lead to some changes to the probe-sets used in noob+funnorm
      ##for the moment these are best practices though
  ##publication also provides updated probe mapping and exclude-infomration
    ##available from: http://zwdzwd.github.io/InfiniumAnnotation

##https://github.com/zhou-lab/InfiniumAnnotation
  ##URL keeps changing. Has been hosted on Amazon AWS
  ##15 Nov 2021 - AWS no longer working? 
    ##try: https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.manifest.rds
    ##try: https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/HM450/HM450.hg38.manifest.rds

##Updated annotations for 450K hg19/38 or 850K hg38 - use frozen version 20180909

#download.file("http://zwdzwd.io/InfiniumAnnotation/20180909/HM450/HM450.hg19.manifest.rds",paste0(HOME,"/me/","HM450.hg19.manifest.rds"))
download.file("http://zwdzwd.io/InfiniumAnnotation/20180909/HM450/HM450.hg38.manifest.rds",paste0(HOME,"/me/","HM450.hg38.manifest.rds"))

md5sum(paste0(HOME,"/me/","HM450.hg38.manifest.rds"))
#"c3086f64841512c2a4f1c1200e723ea2" 

##450K 2021 - resource migrated..
#download.file("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/HM450/HM450.hg19.manifest.tsv.gz",paste0(HOME,"/me/","HM450.hg19.manifest.tsv.gz"))
#download.file("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/HM450/HM450.hg38.manifest.tsv.gz",paste0(HOME,"/me/","HM450.hg38.manifest.tsv.gz"))

##EPIC
#download.file("http://zwdzwd.io/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.manifest.rds",paste0(HOME,"/me/","EPIC.hg38.manifest.rds"))

##EPIC 2021 - resource migrated..
#download.file("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.manifest.tsv.gz",paste0(HOME,"/me/","EPIC.hg38.manifest.tsv.gz"))
#download.file("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/EPIC/EPIC.hg19.manifest.tsv.gz",paste0(HOME,"/me/","EPIC.hg19.manifest.tsv.gz"))

##Annotations:
#seqnames, start and end - the location of the target (1-based coordinates, 2 nucleotides for CpG probes, or 1 nucleotide for CpH and SNP probes). strand is left as "*" always. Some erroneous CpH probe coordinates mapping information in the manufacturer's manifest have been corrected. SNP probe coordinates are provided.
#addressA and addressB - addresses of probe A and B on the chip designated by the original manifest.
#channel - "Both" for type II probes and "Grn"/"Red" for type I probes.
#designType - either "I" or "II".
#nextBase - the actual extension base (on the probe strand) after bisulfite conversion ("A" or "C" or "T"). Unmapped probe has extension base labeled in the original manifest.
#nextBaseRef - the extension base (on the hybridized DNA) before bisulfite conversion ("A", "C", "G" or "T"). Unmapped probe has "NA".
#probeType - either "cg", "ch" or "rs".
#orientation - either "up" or "down" specifying whether the probe is positioned upstream or downstream the target.
#probeCpGcnt - the number of CpG in the probe.
#context35 - the number of CpG in the [-35bp, +35bp] window.
#probeStart and probeEnd - the mapped start and end position of the probe, it is always 50bp long.
#ProbeSeqA and ProbeSeqB - the probe sequence for allele A and B.
#chrmA, begA, flag.A, mapQ.A, cigarA - the mapping info for probe A. mapQ=mapping quality score, 0-60, with 60 being the best
#chrmB, begB, flag.B, mapQ.B, cigarB - the mapping info for probe B.
#posMatch - whether the mapping matches the original manifest.

##Masking categories:
#MASK.mapping - whether the probe is masked for mapping reason. Probes retained should have high quality (>10 on 0-60 scale) consistent (with designed MAPINFO) mapping (for both in the case of type I) without INDELs.
#MASK.typeINextBaseSwitch - whether the probe has a SNP in the extension base that causes a color channel switch from the official annotation (described as color-channel-switching, or CCS SNP in the reference). These probes should be processed differently than designed (by summing up both color channels instead of just the annotated color channel).
#MASK.rmsk15 - whether the 15bp 3'-subsequence of the probe overlap with repeat masker, this MASK is NOT recommended.
#MASK.sub25.copy, MASK.sub30.copy, MASK.sub35.copy and MASK.sub40.copy - whether the 25bp, 30bp, 35bp and 40bp 3'-subsequence of the probe is non-unique.
#MASK.snp5.common - whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the common SNPs from dbSNP (global MAF can be under 1%).
#MASK.snp5.GMAF1p - whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the SNPs with global MAF >1%.
#MASK.extBase - probes masked for inconsistent extension base based on mapping.
#MASK.general - recommended general purpose masking merged from "MASK.sub30.copy", "MASK.mapping", "MASK.extBase", "MASK.typeINextBaseSwitch" and "MASK.snp5.GMAF1p".

##get updated annotations
hm.manifest<-readRDS( paste0(HOME,"/me/","HM450.hg38.manifest.rds") )

ls()
#  [1] "data.atac.core"      "data.cn.core"        "data.counts.core"   
#  [4] "data.fpkm.core"      "data.mutations.core" "data.seg.core"      
#  [7] "data.uq.core"        "hm.manifest"         "HOME"               
# [10] "MANIFEST"            "methSet"             "mutations.core"     
# [13] "rgSetFnorm"          "samples.me.core"     "TUMOR_TYPE"       

##fairly large number of probes recommended for exclusion
table(hm.manifest$MASK_general)
#  FALSE   TRUE 
# 421433  64144 

##infinium I vs II
table(hm.manifest$designType,hm.manifest$MASK_general)
  #     FALSE   TRUE
  # I  115147  20354
  # II 306286  43790

table(hm.manifest$designType,hm.manifest$probeType)
  #        cg     ch     rs
  # I  135476      0     25
  # II 346945   3091     40

table(hm.manifest$MASK_general,hm.manifest$probeType)
  #           cg     ch     rs
  # FALSE 418858   2510     65
  # TRUE   63563    581      0

##figure out how to treat CpH-probes later...

################################################################################
##fix data

##create sample annotation object
sampleMat<-pData(rgSetFnorm)

sampleMat<-as.data.frame(sampleMat)

str(sampleMat)
# 'data.frame':	669 obs. of  6 variables:
#  $ Sample_ID  : chr  "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
#  $ Sample_Type: chr  "Primary Solid Tumor" "Primary Solid Tumor" "Primary Solid Tumor" "Primary Solid Tumor" ...
#  $ Slide      : chr  "9993943013" "9993943013" "9993943005" "9993943017" ...
#  $ Array      : chr  "R04C01" "R01C02" "R02C02" "R01C01" ...
#  $ Basename   : chr  "/tmp/Rtmp505mX4/arrays/9993943013/9993943013_R04C01" "/tmp/Rtmp505mX4/arrays/9993943013/9993943013_R01C02" "/tmp/Rtmp505mX4/arrays/9993943005/9993943005_R02C02" "/tmp/Rtmp505mX4/arrays/9993943017/9993943017_R01C01" ...
#  $ filenames  : chr  "/tmp/Rtmp505mX4/arrays/9993943013/9993943013_R04C01" "/tmp/Rtmp505mX4/arrays/9993943013/9993943013_R01C02" "/tmp/Rtmp505mX4/arrays/9993943005/9993943005_R02C02" "/tmp/Rtmp505mX4/arrays/9993943017/9993943017_R01C01" ...

##get Beta and M-values for data set
betaData <- getBeta(rgSetFnorm)

str( betaData )
 # num [1:485512, 1:669] 0.7205 0.7713 0.0351 0.1073 0.0277 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:485512] "cg13869341" "cg14008030" "cg12045430" "cg20826792" ...
 #  ..$ : chr [1:669] "9993943013_R04C01" "9993943013_R01C02" "9993943005_R02C02" "9993943017_R01C01" ...

##remove "unreliable" probes
probesKeep<-cbind(id=names(hm.manifest),
	chr=as.character(seqnames(hm.manifest)),
	start=as.integer(start(hm.manifest)),
	end=as.integer(end(hm.manifest)),
  isBad=hm.manifest$MASK_general,
  probeType=hm.manifest$probeType,
  designType=hm.manifest$designType)

head(probesKeep)
#      id           chr    start   end     isBad  probeType designType
# [1,] "cg13869341" "chr1" "15865" "15866" "TRUE" "cg"      "I"       
# [2,] "cg14008030" "chr1" "18827" "18828" "TRUE" "cg"      "II"      
# [3,] "cg12045430" "chr1" "29407" "29408" "TRUE" "cg"      "I"       
# [4,] "cg20826792" "chr1" "29425" "29426" "TRUE" "cg"      "I"       
# [5,] "cg00381604" "chr1" "29435" "29436" "TRUE" "cg"      "I"       
# [6,] "cg20253340" "chr1" "68849" "68850" "TRUE" "cg"      "II"      

probesKeep<-as.data.frame(probesKeep,stringsAsFactors=FALSE)
rownames(probesKeep)<-probesKeep$id

str(probesKeep)
# 'data.frame':	485577 obs. of  7 variables:
#  $ id        : chr  "cg13869341" "cg14008030" "cg12045430" "cg20826792" ...
#  $ chr       : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start     : chr  "15865" "18827" "29407" "29425" ...
#  $ end       : chr  "15866" "18828" "29408" "29426" ...
#  $ isBad     : chr  "TRUE" "TRUE" "TRUE" "TRUE" ...
#  $ probeType : chr  "cg" "cg" "cg" "cg" ...
#  $ designType: chr  "I" "II" "I" "I" ...

table(probesKeep$isBad,probesKeep$designType)
  #            I     II
  # FALSE 115147 306286
  # TRUE   20354  43790

probesKeep$start<-as.integer(probesKeep$start)
probesKeep$end<-as.integer(probesKeep$end)
probesKeep$isBad<-as.logical(probesKeep$isBad)

##check default platform ids versus uppdated annotations
table( probesKeep$id %in% rownames(betaData) )
 # FALSE   TRUE 
 #    65 485512 

table( sub("\\d+","",probesKeep$id[ ! probesKeep$id %in% rownames(betaData) ]) )
# rs 
# 65 

table( rownames(betaData) %in% probesKeep$id )
#   TRUE 
# 485512 

length( intersect( rownames(betaData),rownames(probesKeep) ) )
#[1] 485512

isect<-intersect( rownames(betaData),rownames(probesKeep) )

betaData<-betaData[isect,]

probesKeep<-probesKeep[isect,]

all.equal( rownames(betaData) , rownames(probesKeep) )
#[1] TRUE

rm(isect)

##beta
dim(betaData)
#[1] 485512    669

sum(probesKeep$isBad)
#[1] 64144

rmIndex<- ! probesKeep$isBad

sum(rmIndex)
#[1] 421368
 
sum(rmIndex)+sum(!rmIndex)
#[1] 485512

##trim beta values
betaData <- betaData[rmIndex,]

dim(betaData)
#[1] 421368    669

##trim probe annotations
probesKeep<-probesKeep[rmIndex,]

dim(probesKeep)
#[1] 421368      7

all.equal( rownames(betaData) , rownames(probesKeep) )
#[1] TRUE

##Only need matrices and annots going forward
rm(rgSetFnorm,rmIndex)

ls()
#  [1] "betaData"            "data.atac.core"      "data.cn.core"       
#  [4] "data.counts.core"    "data.fpkm.core"      "data.mutations.core"
#  [7] "data.seg.core"       "data.uq.core"        "hm.manifest"        
# [10] "HOME"                "MANIFEST"            "methSet"            
# [13] "mutations.core"      "probesKeep"          "sampleMat"          
# [16] "samples.me.core"     "TUMOR_TYPE"         

rm(methSet,hm.manifest)

gc()

################################################################################
###Check if batch effects present in data set using Martin Lauss "swamp"-package

str(sampleMat)
# 'data.frame':	669 obs. of  6 variables:
#  $ Sample_ID  : chr  "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
#  $ Sample_Type: chr  "Primary Solid Tumor" "Primary Solid Tumor" "Primary Solid Tumor" "Primary Solid Tumor" ...
#  $ Slide      : chr  "9993943013" "9993943013" "9993943005" "9993943017" ...
#  $ Array      : chr  "R04C01" "R01C02" "R02C02" "R01C01" ...
#  $ Basename   : chr  "/tmp/Rtmp505mX4/arrays/9993943013/9993943013_R04C01" "/tmp/Rtmp505mX4/arrays/9993943013/9993943013_R01C02" "/tmp/Rtmp505mX4/arrays/9993943005/9993943005_R02C02" "/tmp/Rtmp505mX4/arrays/9993943017/9993943017_R01C01" ...
#  $ filenames  : chr  "/tmp/Rtmp505mX4/arrays/9993943013/9993943013_R04C01" "/tmp/Rtmp505mX4/arrays/9993943013/9993943013_R01C02" "/tmp/Rtmp505mX4/arrays/9993943005/9993943005_R02C02" "/tmp/Rtmp505mX4/arrays/9993943017/9993943017_R01C01" ...

g<-betaData

o<-as(sampleMat[,c("Slide",
  "Array"
)],"data.frame")

str(o)
# 'data.frame':	669 obs. of  3 variables:
#  $ Slide      : chr  "9993943013" "9993943013" "9993943005" "9993943017" ...
#  $ Array      : chr  "R04C01" "R01C02" "R02C02" "R01C01" ...

##all are basically categorical
o<-data.frame(lapply(o,as.factor))
rownames(o)<-paste(o$Slide,o$Array,sep="_")

res1<-prince(g,o,top=15,permute=FALSE)

pdf(file=paste0(HOME,"/me/","princePlotPostMinfiCorrection.pdf"),width=12,height=12,useDingbats=FALSE)
prince.plot(prince=res1,notecex=1,cexRow=1.5,cexCol=1.5,note=T)
dev.off()

##seems to have done a good job of adjusting for technical variation

rm(g)
rm(o)
rm(res1)

gc()

################################################################################
###Check if infinium I vs II scaling problem is fixed by normalization

ls()
#  [1] "betaData"            "data.atac.core"      "data.cn.core"       
#  [4] "data.counts.core"    "data.fpkm.core"      "data.mutations.core"
#  [7] "data.seg.core"       "data.uq.core"        "HOME"               
# [10] "MANIFEST"            "mutations.core"      "probesKeep"         
# [13] "sampleMat"           "samples.me.core"     "TUMOR_TYPE"        

str(betaData)
 # num [1:421368, 1:669] 0.642 0.683 0.206 0.816 0.768 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:421368] "cg21870274" "cg08258224" "cg16619049" "cg18147296" ...
 #  ..$ : chr [1:669] "9993943013_R04C01" "9993943013_R01C02" "9993943005_R02C02" "9993943017_R01C01" ...

str(probesKeep)
# 'data.frame':	421368 obs. of  7 variables:
#  $ id        : chr  "cg21870274" "cg08258224" "cg16619049" "cg18147296" ...
#  $ chr       : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start     : int  69591 864703 870161 877159 898803 898915 898976 902156 902999 903106 ...
#  $ end       : int  69592 864704 870162 877160 898804 898916 898977 902157 903000 903107 ...
#  $ isBad     : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
#  $ probeType : chr  "cg" "cg" "cg" "cg" ...
#  $ designType: chr  "II" "II" "I" "II" ...

all.equal(rownames(betaData),probesKeep$id)
#[1] TRUE

pdf(file=paste0(HOME,"/me/","InfiniumIvsIIbetas.pdf"),width=12,height=12,useDingbats=F)
plot(density(betaData[probesKeep$designType=="I",1]),main=colnames(betaData)[1])
lines(density(betaData[probesKeep$designType=="II",1]),col=2)
legend("topright",c("Infinium I","Infinium II"),col=1:2,bty="n",lwd=3,cex=3)
dev.off()

##calculate distance from 1 or zero
fn1<-function(x) {
	z1<- 1-x <= .5
	z2<- x < .5
	z3<-vector(length=length(x))
	z3[z1]<-(1-x)[z1]
	z3[z2]<-(x)[z2]
	median(z3)
}

bI<-apply(betaData[probesKeep$designType=="I",],2,fn1)
bII<-apply(betaData[probesKeep$designType=="II",],2,fn1)

quantile(bI)
#        0%        25%        50%        75%       100% 
#0.05115658 0.07303714 0.07896031 0.08527074 0.13987000 

quantile(bII)
#       0%       25%       50%       75%      100% 
#0.1026016 0.1602946 0.1793252 0.1925217 0.2791591 

pdf(file=paste0(HOME,"/me/","InfiniumIvsIIdistanceToZeroOne.pdf"),width=12,height=12,useDingbats=F)
plot(density(bI),col=1,xlim=c(0,.5),main=colnames(betaData)[1],sub="median distance to zero or one")
lines(density(bII),col=2)
legend("topright",c("Infinium I","Infinium II"),col=1:2,bty="n",lwd=3,cex=3)
dev.off()

rm(bI,bII,fn1)

##difference in median distance to 0-1 between infinium I and II probes

################################################################################
##Recalibrate infinium betas

##Beta adjustment script below adapted from Markus Ringn\E9r's (Oncology & Pathology, LU) script "beta_peak_normalize.R"
  ##Citation: Holm et al. 2016 (Breast Cancer Res)

##known that infinium I and II have different dynamics leading to value-truncation for one of the infinium II probes
table(probesKeep$designType)
#      I     II 
# 115122 306246 

##define probe subsets
probes_I<- probesKeep$designType == "I"
probes_II<- probesKeep$designType == "II"

table(probes_I,probes_II)
#         probes_II
# probes_I  FALSE   TRUE
#    FALSE      0 306246
#    TRUE  115122      0

##identify zero and one peak for respective probe sets and calibrate range of one to other

##derive function for detecting (local) maxima (modified from Ringner script)
#findMaxima<-function(x,nMax=2,granularity=512,bw=.02,adjust=1,from=0,to=1,kernel="epanechnikov") {
#  aa<-density(x,kernel=kernel,n=granularity,bw=bw,adjust=adjust,from=from,to=to,na.rm=T)
#  res<-vector(length=length(aa$x))
#  for(i in 2:length(res[-1])) {
#    res[i]<-aa$y[i]>aa$y[i-1] & aa$y[i]>aa$y[i+1]
#  }
#  aa$x[res][order(aa$y[res],decreasing=T)][1:nMax]
#} ##choose 2 main peaks -> probably near zero and one

##updated function for handling samples with non-canonical beta distributions
findMaxima<-function(x,granularity=512,bw=.02,adjust=1,from=0,to=1,kernel="epanechnikov",verbose=TRUE) {
  aa<-density(x,kernel=kernel,n=granularity,bw=bw,adjust=adjust,from=from,to=to,na.rm=T)
  res<-vector(length=length(aa$x))
  for(i in 2:length(res[-1])) {
    res[i]<-aa$y[i]>aa$y[i-1] & aa$y[i]>aa$y[i+1]
  }
  y<-aa$x[res][order(aa$y[res],decreasing=T)] 
	if(verbose) {
		cat("found peaks at:",y,"\n")
		cat("choosing min/max:",range(y),"\n") ###CHANGED here to pick out min-max peaks as these are typically the correct ones
	}
	range(y) 
}

##define results matrix
meth.cal<-matrix(0,ncol=ncol(betaData),nrow=nrow(betaData))
rownames(meth.cal)<-rownames(betaData)
colnames(meth.cal)<-colnames(betaData)

##define vector for flagging IDs
flaggedIDs<-NULL

pdf(file=paste0(HOME,"/me/","infiniumAdjustmentPlots.pdf"),width=8,height=8)
par(mfrow=c(2,2),font=2,font.sub=2,font.lab=2,font.axis=2,las=1)
sink(paste0(HOME,"/me/","infiniumAdjustmentStats.txt"))

for ( i in 1:ncol(betaData)) {
  cat("#####Sample ",i,":",colnames(betaData)[i],"\n")
  ##do pre-calibaration I vs II plot
  plot(density(betaData[probes_I,i],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T),main=colnames(betaData)[i],xlim=c(-.05,1.05),col="orange",lwd=3,xlab="pre-adjustment beta",cex.main=.9,ylim=c(0,8))
  lines(density(betaData[probes_II,i],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T),col="darkgreen",lwd=3)

  ##define maxima
  tempMax_I<-findMaxima(betaData[probes_I,i])
  tempMax_II<-findMaxima(betaData[probes_II,i])
  ##function not oriented so peak near zero can be in vector slot 2 -> fix
  tempMax_I<-tempMax_I[order(tempMax_I)]
  tempMax_II<-tempMax_II[order(tempMax_II)]
  ##add plot "FLAG" if peaks not well separated and close to right place..
  FLAG<-FALSE
  if( tempMax_I[1] >.2 | tempMax_I[2] <.8 | tempMax_II[1] >.2 | tempMax_II[2] <.8 ) {
    FLAG<-TRUE
    flaggedIDs<-c(flaggedIDs,colnames(betaData)[i])
  }

  ##add to plot
  abline(v=tempMax_I,col="orange",lwd=2,lty=2)
  abline(v=tempMax_II,col="darkgreen",lwd=2,lty=2)
  text(x=.5,y=c(3.5,3),labels=c("Infinium_I ","Infinium_II"),col=c("orange","darkgreen"),font=2,cex=1)
  if(FLAG)   text(x=.5,y=.5,labels=c("FLAGGED"),col=2,font=2,cex=2)

  if(FLAG) {
    cat("Something funky with this sample!!!\n")
    cat("infinium_I peaks:",tempMax_I,"\n",sep="\t")
    cat("infinium_II peaks:",tempMax_II,"\n",sep="\t")
  } else {
    cat("infinium_I peaks:",tempMax_I,"\n",sep="\t")
    cat("infinium_II peaks:",tempMax_II,"\n",sep="\t")
  }

  ##do calibration scale both to same range..
  meth.cal[probes_I,i]<-(betaData[probes_I,i]-tempMax_I[1])/(tempMax_I[2]-tempMax_I[1])
  meth.cal[probes_II,i]<-(betaData[probes_II,i]-tempMax_II[1])/(tempMax_II[2]-tempMax_II[1])

  ##cap ends so that nothing over 1 or below zero..
  meth.cal[probes_I,i][meth.cal[probes_I,i] < 0 ]<-0
  meth.cal[probes_I,i][meth.cal[probes_I,i] > 1 ]<-1
  meth.cal[probes_II,i][meth.cal[probes_II,i] < 0 ]<-0
  meth.cal[probes_II,i][meth.cal[probes_II,i] > 1 ]<-1

  ##do post-calibration I vs II plot
  plot(density(meth.cal[probes_I,i],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T),main=colnames(meth.cal)[i],xlim=c(-.05,1.05),col="orange",lwd=3,xlab="post-adjustment beta",cex.main=.9,ylim=c(0,8))
  lines(density(meth.cal[probes_II,i],kernel="epanechnikov",bw=0.02,from=0,to=1,na.rm=T),col="darkgreen",lwd=3)
  text(x=.5,y=c(3.5,3),labels=c("Infinium_I ","Infinium_II"),col=c("orange","darkgreen"),font=2,cex=1)
  abline(v=findMaxima(meth.cal[probes_I,i],verbose=F),col="orange",lwd=2,lty=2)
  abline(v=findMaxima(meth.cal[probes_II,i],verbose=F),col="darkgreen",lwd=2,lty=2)

  cat("\n")
} ; rm( list = c("i","tempMax_I","tempMax_II","FLAG") )
cat("##FlaggedIDs:\n")
cat(flaggedIDs,sep="\n")
sink()
dev.off()

##some peaks not where they should be??
flaggedIDs
#[1] "9993943017_R06C02"

writeLines(flaggedIDs,con=paste0(HOME,"/me/","infiniumAdjustmentFlaggedIDs.txt"))

if(!is.null(flaggedIDs)) {
	write(paste(sampleMat[flaggedIDs,"Sample_ID"],flaggedIDs,sep="\t"),file=paste0(HOME,"/me/","methylationSampleBlacklist.txt"),append=TRUE)
}

rm(probes_I)
rm(probes_II)
rm(findMaxima)

##write adjusted table to file
con<-gzfile(paste0(HOME,"/me/","infiniumAdjustedNormalizedData.txt.gz"),"w")
write.table(meth.cal,file=con,sep="\t",quote=F)
close(con)
rm(con)

##HERE!!

################################################################################
###Check magnitude of beta change for probes

pdf(file=paste0(HOME,"/me/","InfiniumIvsIIadjustmentMagnitude.pdf"),width=12,height=12,useDingbats=FALSE)
par(mfrow=c(2,1))

##median shift in beta
res<-NULL

for (i in 1:ncol(betaData)) {
  res<-c(res, median( abs(betaData[,i]-meth.cal[,i]) ) )
}
quantile(res)
#         0%        25%        50%        75%       100% 
# 0.02571348 0.03745477 0.04002998 0.04284042 0.05921211 

plot(res,col=c(1,2)[1+as.integer(colnames(meth.cal) %in% flaggedIDs)],ylab="MAD",
	main="Infinium I/II adjustment",sub="median absolute beta shift",pch=16)
legend("topright",c("Flagged"),col=2,bty="n",pch=16,cex=1)

##pre/post sample corrlations
res<-NULL

for (i in 1:ncol(betaData)) {
  res<-c(res, cor(betaData[,i],meth.cal[,i]) )
}
quantile(res)
#       0%       25%       50%       75%      100% 
#0.9991500 0.9994938 0.9995245 0.9995539 0.9996684 

plot(res,col=c(1,2)[1+as.integer(colnames(meth.cal) %in% flaggedIDs)],ylab="Correlation",
	main="Infinium I/II adjustment",sub="correlation pre/post adjustment",pch=16)
legend("bottomright",c("Flagged"),col=2,bty="n",pch=16,cex=1)

dev.off()

table(res<.98)
# FALSE 
#   669 

##individual probe correlations
tmp1<-apply(betaData,1,sd)!=0
tmp2<-apply(meth.cal,1,sd)!=0

table(tmp1,tmp2) ##capping at 0/1 leads to sd=0 rows
#       tmp2
# tmp1    FALSE   TRUE
#   TRUE   1834 419534

##no zero SD rows in original data due to there being some signal always

set.seed(20201004)
rowSample<-sample(which(tmp1 & tmp2),10000)

res<-NULL

for (i in rowSample) {
  res<-c(res, cor(betaData[i,],meth.cal[i,]) )
}
quantile(res)
#        0%       25%       50%       75%      100% 
# 0.1085867 0.9633940 0.9946673 0.9981274 0.9999691 

lapply(split(res,probesKeep$designType[rowSample]),quantile)
# $I
#        0%       25%       50%       75%      100% 
# 0.1085867 0.9623843 0.9976341 0.9997909 0.9999691 

# $II
#        0%       25%       50%       75%      100% 
# 0.1100834 0.9637807 0.9938648 0.9976329 0.9994696 

pdf(file=paste0(HOME,"/me/","InfiniumIvsIIadjustmentProbeCorr.pdf"),width=12,height=12,useDingbats=F)
par(mfrow=c(2,1))

boxplot(res~probesKeep$designType[rowSample],las=1,cex.sub=1.5,
	lwd=2,ylim=c(0,1),xlab="",main="10K random probes",sub="Infinium I and II pre/post adj correlation")

##effect of flagged samples?
res<-NULL

for (i in rowSample) {
  res<-c(res, cor(betaData[i,!colnames(meth.cal) %in% flaggedIDs],meth.cal[i,!colnames(meth.cal) %in% flaggedIDs]) )
}
quantile(res,na.rm=T)
#        0%        25%        50%        75%       100% 
#0.07193345 0.98950482 0.99864126 0.99959297 0.99994293 

lapply(split(res,probesKeep$designType[rowSample]),quantile,na.rm=T)
# $I
#        0%       25%       50%       75%      100% 
# 0.1085842 0.9624229 0.9976607 0.9997929 0.9999693 

# $II
#        0%       25%       50%       75%      100% 
# 0.1100664 0.9654577 0.9942198 0.9978137 0.9994689 

boxplot(res~probesKeep$designType[rowSample],las=1,cex.sub=1.5,
	lwd=2,ylim=c(0,1),xlab="",main="10K random probes minus flagged samples",sub="Infinium I and II pre/post adj correlation")

dev.off()

rm(res,rowSample,tmp1,tmp2,i)

################################################################################
###PCA on calibrated mat to see that large effects not introduced

##with bad calibaration samples
g<-meth.cal

o<-as(sampleMat[,c("Slide",
  "Array"
)],"data.frame")

str(o)
# 'data.frame':	669 obs. of  3 variables:
#  $ Slide      : chr  "9993943013" "9993943013" "9993943005" "9993943017" ...
#  $ Array      : chr  "R04C01" "R01C02" "R02C02" "R01C01" ...

##all are basically categorical
o<-data.frame(lapply(o,as.factor))
rownames(o)<-paste(o$Slide,o$Array,sep="_")

res1<-prince(g,o,top=15,permute=FALSE)

pdf(file=paste0(HOME,"/me/","princePlotPostInfiniumScaling.pdf"),width=12,height=12,useDingbats=FALSE)
prince.plot(prince=res1,notecex=1,cexRow=1.5,cexCol=1.5,note=T)
dev.off()

##no major changes to variance structure..

rm(g)
rm(o)
rm(res1)

gc()

################################################################################
##check data and move to save..

ls()
#  [1] "betaData"            "data.atac.core"      "data.cn.core"       
#  [4] "data.counts.core"    "data.fpkm.core"      "data.mutations.core"
#  [7] "data.seg.core"       "data.uq.core"        "flaggedIDs"         
# [10] "HOME"                "MANIFEST"            "meth.cal"           
# [13] "mutations.core"      "probesKeep"          "sampleMat"          
# [16] "samples.me.core"     "TUMOR_TYPE"         

rm(flaggedIDs)

all.equal(rownames(sampleMat),colnames(meth.cal))
#[1] TRUE
all.equal(rownames(sampleMat),colnames(betaData))
#[1] TRUE

all.equal(rownames(meth.cal),rownames(betaData))
#[1] TRUE

dim(meth.cal)
#[1] 421368    669

dim(betaData)
#[1] 421368    669

################################################################################
##save and exit

ls()
#  [1] "betaData"            "data.atac.core"      "data.cn.core"       
#  [4] "data.counts.core"    "data.fpkm.core"      "data.mutations.core"
#  [7] "data.seg.core"       "data.uq.core"        "HOME"               
# [10] "MANIFEST"            "meth.cal"            "mutations.core"     
# [13] "probesKeep"          "sampleMat"           "samples.me.core"    
# [16] "TUMOR_TYPE"     

##beta matrix 
con<-gzfile(paste0(HOME,"/me/","data450k_minfiNormalizedBeta.txt.gz"),"w")
write.table(betaData,file=con,sep="\t",quote=F)
close(con)
rm(con)

##adjusted beta matrix 
con<-gzfile(paste0(HOME,"/me/","data450k_minfiNormalized_adjustedBeta.txt.gz"),"w")
write.table(meth.cal,file=con,sep="\t",quote=F)
close(con)
rm(con)

##sample matrix 
write.table(sampleMat,file=paste0(HOME,"/me/","minfiPipeline_sampleMatrix.txt"),sep="\t",quote=FALSE)

##save three main objects
betaFinal<-round(meth.cal,3)

all.equal(rownames(betaFinal),probesKeep$id)
#[1] TRUE
all.equal(colnames(betaFinal),rownames(sampleMat))
#[1] TRUE

betaFinal[1:4,1:5]
#            9993943013_R04C01 9993943013_R01C02 9993943005_R02C02
# cg21870274             0.718             0.443             0.731
# cg08258224             0.766             0.749             0.975
# cg16619049             0.189             0.451             0.585
# cg18147296             0.925             0.832             0.929
#            9993943017_R01C01 9993943013_R02C01
# cg21870274             0.912             0.682
# cg08258224             0.983             0.940
# cg16619049             0.565             0.371
# cg18147296             0.868             0.900

##rounded beta matrix 
con<-gzfile(paste0(HOME,"/me/","data450k_minfiNormalized_adjustedBeta_round3dec.txt.gz"),"w")
write.table(betaFinal,file=con,sep="\t",quote=F)
close(con)
rm(con)

save(betaFinal,sampleMat,probesKeep,file=paste0(HOME,"/me/","workspace_minfiNormInfiniumAdjFinalBetas.RData"))

save(probesKeep,file=paste0(HOME,"/me/","object_450k_probesKeep.RData"))

##get md5 to check reproducibility 
	##use final rounded matrix for checking as others may differ due to default rounding options and floating point diffs or similar...

md5sum(paste0(HOME,"/me/","data450k_minfiNormalized_adjustedBeta_round3dec.txt.gz"))
# I:/data/tcgaBrca/me/data450k_minfiNormalized_adjustedBeta_round3dec.txt.gz 
#                                         "edf08072a379ae91bf7d6b8eefa51bf8" 

q("no")
###END