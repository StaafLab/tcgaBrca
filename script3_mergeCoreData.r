#####======================================================================#####
### Get TCGA/GDC data for TCGA BRCA
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
#dir.create(paste0(HOME,"/","merge"),recursive=TRUE)

################################################################################
##load required packages

################################################################################
##Define functions

##add download header
parseTcgaName<-function(x="TCGA-NN-NNNN-NNN-NNN-NNNN-NN",level="sample_type") {
	x<-as.character(x)
	if(level=="1") { return( sub("(TCGA)-..-....-...-...-....-..","\\1",x) ) }
	if(level=="2") { return( sub("(TCGA-..)-....-...-...-....-..","\\1",x) ) }
	if(level=="3") { return( sub("(TCGA-..-....)-...-...-....-..","\\1",x) ) }
	if(level=="4") { return( sub("(TCGA-..-....-...)-...-....-..","\\1",x) ) }
	if(level=="5") { return( sub("(TCGA-..-....-...-...)-....-..","\\1",x) ) }
	if(level=="6") { return( sub("(TCGA-..-....-...-...-....)-..","\\1",x) ) }
	if(level=="7") { return( sub("(TCGA-..-....-...-...-....-..)","\\1",x) ) }	
	if(level=="patient") { return( sub("(TCGA-..-....)-...-...-....-..","\\1",x) ) }	
	if(level=="sample") { return( sub("(TCGA-..-....-..).-...-....-..","\\1",x) ) }	
	if(level=="batch") { return( sub("TCGA-..-....-...-...-(....)-..","\\1",x) ) }	
	if(level=="sample_type") { return( sub("TCGA-..-....-(..).-...-....-..","\\1",x) ) }	
	if(level=="atac_patient") { return( sub("(TCGA-..-....)-...-..-....-..","\\1",x) ) }	
	if(level=="atac_sample") { return( sub("(TCGA-..-....-..).-..-....-..","\\1",x) ) }	
	if(level=="atac_type") { return( sub("TCGA-..-....-(..).-..-....-..","\\1",x) ) }	
	if(level=="atac_id") { return( sub("(TCGA-..-....-...)-..-....-..","\\1",x) ) }	
}

################################################################################
##get core data levels and merge

ls()
#[1] "HOME"       "MANIFEST"   "TUMOR_TYPE"

load(file=paste0(HOME,"/","coreData_gexCnWesAtac_unique_samples.RData"))

ls()
#  [1] "data.atac"           "data.cn"             "data.counts"        
#  [4] "data.fpkm"           "data.seg"            "data.uq"            
#  [7] "HOME"                "MANIFEST"            "mutations.consensus"
# [10] "TUMOR_TYPE"         

samples.me<-read.table(file=paste0(HOME,"/","me/","me_unique_samples.txt"),sep="\t",as.is=TRUE,col.names=c("tcga_id","array"))

str(samples.me)
# 'data.frame':   893 obs. of  2 variables:
#  $ tcga_id: chr  "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
#  $ array  : chr  "9993943013_R04C01" "9993943013_R01C02" "9993943005_R02C02" "9993943017_R01C01" ...
 
################################################################################
##get core data level overlaps

##atac
s.atac<-colnames(data.atac)

##all seq-arrays have same samples
s.gex<-colnames(data.counts)
s.gex2<-colnames(data.fpkm)
s.gex3<-colnames(data.uq)

all.equal(colnames(data.counts),colnames(data.fpkm))
#[1] TRUE
all.equal(colnames(data.counts),colnames(data.uq))
#[1] TRUE

##cn
s.cn<-colnames(data.cn)
s.cn2<-unique(data.seg[,"Sample_ID"])

all.equal(intersect(s.cn,s.cn2),s.cn)
#[1] TRUE
all.equal(intersect(s.cn2,s.cn2),s.cn2)
#[1] TRUE

##me
s.me<-samples.me$tcga_id

##wes
s.wes<-parseTcgaName(unique(mutations.consensus$Tumor_Sample_Barcode),level=4)

length(unique(s.atac))==length(s.atac)
#[1] TRUE

length(unique(s.gex))==length(s.gex)
#[1] TRUE
length(unique(s.gex2))==length(s.gex2)
#[1] TRUE
length(unique(s.gex3))==length(s.gex3)
#[1] TRUE

length(unique(s.cn))==length(s.cn)
#[1] TRUE
length(unique(s.cn2))==length(s.cn2)
#[1] TRUE

length(unique(s.me))==length(s.me)
#[1] TRUE

length(unique(s.wes))==length(s.wes)
#[1] TRUE

s.all<-intersect(intersect(intersect(s.gex,s.cn),s.me),s.wes)

table( sub(".+\\-","",s.all) )
# 01A 01B 
# 666   3 

length( unique(sub("-...$","",s.all)) )
#[1] 669

##intersect with atac
length(intersect(s.atac,s.all))
#[1] 39

length(intersect(s.atac,s.gex))
#[1] 39
length(intersect(s.atac,s.cn))
#[1] 39
length(intersect(s.atac,s.me))
#[1] 39
length(intersect(s.atac,s.wes))
#[1] 67

s.all.atac<-intersect(intersect(intersect(intersect(s.gex,s.cn),s.me),s.wes),s.atac)

table( sub(".+\\-","",s.atac) )
# 01A 01B 
#  73   1 

table( sub(".+\\-","",s.all.atac) )
# 01A 
#  39 

################################################################################
##filter down to core data sets

##
ls()
#  [1] "data.atac"           "data.cn"             "data.counts"        
#  [4] "data.fpkm"           "data.seg"            "data.uq"            
#  [7] "HOME"                "MANIFEST"            "mutations.consensus"
# [10] "parseTcgaName"       "s.all"               "s.all.atac"         
# [13] "s.atac"              "s.cn"                "s.cn2"              
# [16] "s.gex"               "s.gex2"              "s.gex3"             
# [19] "s.me"                "s.wes"               "samples.me"         
# [22] "TUMOR_TYPE"         

##cn
data.cn.core<-data.cn[,colnames(data.cn) %in% s.all]
mode(data.cn.core)<-"integer"

table(rowSums(is.na(data.cn.core)))
#     0     1     2     3     4     5     6     7     8     9    10    12    13 
# 57219  1940   339   147    58    68    55    52     2     3     4     2     4 
#    14    15    16    17    20    21    22    32    37    45    64   259   455 
#     7     1    29    12     2    42    12     1     1     1     1    41    10 
#   461   669 
#     1   392 
data.cn.core<-as.data.frame(data.cn.core,stringsAsFactors=FALSE)
data.cn.core<-data.cn.core[,s.all]

data.seg.core<-data.seg[data.seg[,"Sample_ID"] %in% s.all,]
data.seg.core<-as.data.frame(data.seg.core,stringsAsFactors=FALSE)
str(data.seg.core)
# 'data.frame':	133229 obs. of  7 variables:
#  $ Sample_ID        : chr  "TCGA-E9-A1R7-01A" "TCGA-E9-A1R7-01A" "TCGA-E9-A1R7-01A" "TCGA-E9-A1R7-01A" ...
#  $ Chromosome       : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ Start            : chr  "62920" "1756148" "3164376" "4354764" ...
#  $ End              : chr  "1741164" "3164291" "4354656" "5111634" ...
#  $ Copy_Number      : chr  "4" "3" "4" "6" ...
#  $ Major_Copy_Number: chr  "3" "2" "3" "5" ...
#  $ Minor_Copy_Number: chr  "1" "1" "1" "1" ...

data.seg.core$Start <- as.integer(data.seg.core$Start )
data.seg.core$End <- as.integer(data.seg.core$End )
data.seg.core$Copy_Number <- as.integer(data.seg.core$Copy_Number )
data.seg.core$Major_Copy_Number <- as.integer(data.seg.core$Major_Copy_Number )
data.seg.core$Minor_Copy_Number <- as.integer(data.seg.core$Minor_Copy_Number )

data.seg.core<-data.seg.core[order(data.seg.core$Sample_ID,
	as.integer( sub("X","23",sub("Y","24",sub("chr","",data.seg.core$Chromosome) ) ) ),
	data.seg.core$Start),
	]

##gex
data.counts.core<-data.counts[,colnames(data.counts) %in% s.all]
data.fpkm.core<-data.fpkm[,colnames(data.fpkm) %in% s.all]
data.uq.core<-data.uq[,colnames(data.uq) %in% s.all]

mode(data.counts.core)<-"integer"
mode(data.fpkm.core)<-"numeric"
mode(data.uq.core)<-"numeric"

data.counts.core<-data.counts.core[,s.all]
data.fpkm.core<-data.fpkm.core[,s.all]
data.uq.core<-data.uq.core[,s.all]

##me 
samples.me.core<-samples.me[samples.me$tcga_id %in% s.all,]

samples.me.core<-samples.me.core[match(s.all,samples.me.core$tcga_id),]

##wes
mutations.core<-mutations.consensus[parseTcgaName(mutations.consensus$Tumor_Sample_Barcode,level=4) %in% s.all,]
mutations.core$Tumor_Sample_Barcode<-parseTcgaName(mutations.core$Tumor_Sample_Barcode,level=4)

mutations.core<-mutations.core[order(mutations.core$Tumor_Sample_Barcode,
	as.integer( sub("X","23",sub("Y","24",sub("chr","",mutations.core$Chromosome))) ),
	as.integer( mutations.core$Start_Position),
	as.integer( mutations.core$End_Position)
),]

ls(pattern="\\.core")
# [1] "data.cn.core"     "data.counts.core" "data.fpkm.core"   "data.seg.core"   
# [5] "data.uq.core"     "mutations.core"   "samples.me.core" 

dim(data.cn.core)
#[1] 60446   669

length(intersect(rownames(data.cn.core),rownames(data.counts.core)))
#[1] 60446
length(intersect(rownames(data.cn.core),rownames(data.fpkm.core)))
#[1] 60446
length(intersect(rownames(data.fpkm.core),rownames(data.counts.core)))
#[1] 60483

dim(data.counts.core)
#[1] 60488   669

setdiff(rownames(data.fpkm.core),rownames(data.counts.core))
#character(0)

setdiff(rownames(data.counts.core),rownames(data.fpkm.core))
# [1] "__no_feature"           "__ambiguous"            "__too_low_aQual"       
# [4] "__not_aligned"          "__alignment_not_unique"

##remove 5 extra rows from "counts" -> harmonize matrices

str(data.counts.core)
 # int [1:60488, 1:669] 651 2 1877 1113 387 250 914 1646 1880 3842 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60488] "ENSG00000000003.13" "ENSG00000000005.5" "ENSG00000000419.11" "ENSG00000000457.12" ...
 #  ..$ : chr [1:669] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
str(data.fpkm.core)
 # num [1:60483, 1:669] 0 0 3.68 0 1.83 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60483] "ENSG00000242268.2" "ENSG00000270112.3" "ENSG00000167578.15" "ENSG00000273842.1" ...
 #  ..$ : chr [1:669] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
dim(data.cn.core)
#[1] 60446   669

tmp<-intersect(rownames(data.fpkm.core),rownames(data.counts.core))

data.counts.core<-data.counts.core[tmp,]
data.fpkm.core<-data.fpkm.core[tmp,]
data.uq.core<-data.uq.core[tmp,]
rm(tmp)

##make mutations matrix - filter low-impact variants
table(mutations.core$Gene %in% sub("\\..+","",rownames(data.fpkm.core)) )
# FALSE  TRUE 
#     2 63765 

mutations.core[which(!mutations.core$Gene %in% sub("\\..+","",rownames(data.fpkm.core)) ),"Gene"]
#[1] "" ""

table(mutations.core$IMPACT)
    # HIGH      LOW MODERATE MODIFIER 
    # 5470    13205    34226    10866 

table(mutations.core$Variant_Classification,mutations.core$IMPACT)
  #                         HIGH   LOW MODERATE MODIFIER
  # 3'Flank                    0     0        0      472
  # 3'UTR                      0     0        0     4485
  # 5'Flank                    0     0        0      408
  # 5'UTR                      0     0        0     1604
  # Frame_Shift_Del         1128     0        0        0
  # Frame_Shift_Ins          437     0        0        0
  # IGR                        0     0        0        2
  # In_Frame_Del               0     0      195        0
  # In_Frame_Ins               0     0       27        0
  # Intron                     0     0        0     2907
  # Missense_Mutation          0     0    34003        3
  # Nonsense_Mutation       2849     0        0        0
  # Nonstop_Mutation          47     0        0        0
  # RNA                        0     0        0      985
  # Silent                     0 12334        0        0
  # Splice_Region              0   871        1        0
  # Splice_Site              966     0        0        0
  # Translation_Start_Site    43     0        0        0

ind<-which(mutations.core$IMPACT %in% c("HIGH","MODERATE"))

mutations.tmp<-mutations.core[ind,]

##grab row in main mutations matrix mapping to position in sample-gene matrix
mutations.tmp<-lapply(
	split(
		paste(mutations.tmp$Gene,ind,sep="|"),
		mutations.tmp$Tumor_Sample_Barcode),
	function(x) {
		y<-match(sub("\\..+","",rownames(data.counts.core)),sub("\\|.+","",x) ) 
		y[!is.na(y)]<-sub(".+\\|","",x)[y[!is.na(y)]]
		y
	})

mutations.tmp<-do.call("cbind",mutations.tmp)
rownames(mutations.tmp)<-sub("\\..+","",rownames(data.counts.core))
rm(ind)

str(mutations.tmp)
 # chr [1:60483, 1:666] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : NULL
 #  ..$ : chr [1:666] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

sum(is.na(mutations.tmp))
#[1] 40243680

mutations.tmp[is.na(mutations.tmp)]<-"0"

sum(mutations.tmp=="0")
#[1] 40243680

##Three samples only have low-impact variants..
mode(mutations.tmp)<-"integer"

mutations.filter.core<-matrix(0,nrow=nrow(data.counts.core),
	ncol=length(s.all),
	dimnames=list(sub("\\..+","",rownames(data.counts.core)),s.all)
	)
mutations.filter.core[rownames(mutations.tmp),colnames(mutations.tmp)]<-mutations.tmp
mode(mutations.filter.core)<-"integer"

( sum(mutations.filter.core==0)-sum(mutations.tmp==0) ) / (3*nrow(mutations.filter.core))
#[1] 1

test<-apply(mutations.filter.core,1,function(x) unique(mutations.core$Gene[x[x!=0]]) )

table( names(test)[rowSums(mutations.filter.core)!=0] == unlist(test) )
#  TRUE 
# 13468 
rm(test)

head(rownames(mutations.filter.core)[order(rowSums(mutations.filter.core),decreasing=TRUE)])
# [1] "ENSG00000121879" "ENSG00000141510" "ENSG00000155657" "ENSG00000039068"
# [5] "ENSG00000107485" "ENSG00000181143"

##mutations about expected?
mutations.core$Hugo_Symbol[match(head(rownames(mutations.filter.core)[order(rowSums(mutations.filter.core),decreasing=TRUE)],12),
	mutations.core$Gene
	)]
 # [1] "PIK3CA" "TP53"   "TTN"    "CDH1"   "GATA3"  "MUC16"  "KMT2C"  "MAP3K1"
 # [9] "SYNE1"  "PTEN"   "RYR2"   "NCOR1" 

rm(mutations.tmp)

data.mutations.core<-mutations.filter.core
rm(mutations.filter.core)

###intersect with with atac data
intersect(colnames(data.atac),s.all)
#  [1] "TCGA-3C-AALJ-01A" "TCGA-4H-AAAK-01A" "TCGA-A2-A0SV-01A" "TCGA-A2-A0SW-01A"
#  [5] "TCGA-A2-A0SX-01A" "TCGA-A2-A0T4-01A" "TCGA-A2-A0T5-01A" "TCGA-A2-A0T6-01A"
#  [9] "TCGA-A2-A0T7-01A" "TCGA-A2-A0YC-01A" "TCGA-A2-A0YD-01A" "TCGA-A2-A0YF-01A"
# [13] "TCGA-A2-A0YG-01A" "TCGA-A2-A0YH-01A" "TCGA-A2-A0YJ-01A" "TCGA-A2-A0YK-01A"
# [17] "TCGA-A2-A0YL-01A" "TCGA-A2-A0YT-01A" "TCGA-A7-A0D9-01A" "TCGA-A7-A13E-01A"
# [21] "TCGA-A7-A13F-01A" "TCGA-AO-A03L-01A" "TCGA-AO-A0JB-01A" "TCGA-AO-A0JG-01A"
# [25] "TCGA-AO-A0JM-01A" "TCGA-AO-A124-01A" "TCGA-AR-A0TP-01A" "TCGA-AR-A0TV-01A"
# [29] "TCGA-AR-A0U0-01A" "TCGA-BH-A0B5-01A" "TCGA-BH-A0BA-01A" "TCGA-BH-A0BZ-01A"
# [33] "TCGA-BH-A0DP-01A" "TCGA-BH-A0DV-01A" "TCGA-BH-A0E0-01A" "TCGA-BH-A0HP-01A"
# [37] "TCGA-BH-A1EV-01A" "TCGA-C8-A8HR-01A" "TCGA-S3-AA0Z-01A"

data.atac.core<-data.atac[,intersect(colnames(data.atac),s.all)]

str(data.atac.core)
 # num [1:215920, 1:39] 1.031 1.667 1.839 1.375 0.769 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:215920] "chr1:17234-17733|BRCA_2" "chr1:180634-181133|BRCA_3" "chr1:181207-181706|BRCA_4" "chr1:183557-184056|BRCA_5" ...
 #  ..$ : chr [1:39] "TCGA-3C-AALJ-01A" "TCGA-4H-AAAK-01A" "TCGA-A2-A0SV-01A" "TCGA-A2-A0SW-01A" ...
 
################################################################################
##gather data

ls(pattern="\\.core")
# [1] "data.atac.core"      "data.cn.core"        "data.counts.core"   
# [4] "data.fpkm.core"      "data.mutations.core" "data.seg.core"      
# [7] "data.uq.core"        "mutations.core"      "samples.me.core"    

save(list=ls(pattern="\\.core"),file=paste0(HOME,"/","coreData_gexCnWesMeAtac_unique_matched_samples.RData"))

##Next step is to normalize methylation
	##Then write function to extract LOH etc info from "data.seg.core" for given locus

q("no")
###END