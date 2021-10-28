
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
dir.create(paste0(HOME,"/","gex"),recursive=TRUE)
dir.create(paste0(HOME,"/","cn"),recursive=TRUE)
dir.create(paste0(HOME,"/","me"),recursive=TRUE)
dir.create(paste0(HOME,"/","wes"),recursive=TRUE)
dir.create(paste0(HOME,"/","atac"),recursive=TRUE)

################################################################################
##load required packages

library(tools)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
     install.packages("BiocManager")
 }

if (!requireNamespace("minfi", quietly = TRUE)) {
     install.packages("minfi")
}

library(minfi)

################################################################################
##Define functions

##separate tab-sep lines 
doSep<-function(x) { unlist(strsplit(x,"\\t")) }

##separate pipe-sep lines 
doSep2<-function(x) { unlist(strsplit(x,"\\|")) }

##change pipe-sep to tab-sep 
doUnTab<-function(x) { gsub("\\|","\t",x) }

##paste GDC file manifest header to download info  --V1
#doUnSep<-function(x,y) { paste0(x,"|",y) }

##remove blank spaces from strings
doUnSpace<-function(x) { gsub("(\\[|\\])","",gsub("( |\\t|\\n)","_",x)) }

##function for file id sanity check  --V1
#doUnPaste<-function(x) { sub(".+\\|(.+)\\t.+\\..+","\\1",x) }

##function for file id sanity check
doUnPaste<-function(x) { sub("\\t.+","",x) }

##function for file id sanity check
getIdat<-function(x) { sub(".+\\t(\\d{10}_R\\d+C\\d+_.+\\idat)\\t.+","\\1",x) }

##translate tumor type to TCGA tumor type
idToTcga<-function(x) {list("brca"="TCGA-BRCA",
	"luad"="TCGA-LUAD",
	"lusc"="TCGA-LUSC"
	)[[x]] 
}

##translate tumor type to TCGA tumor type
sampleDecode<-function(x) { 
	l<-list("01"="Primary Solid Tumor",
	"02"="Recurrent Solid Tumor",
	"03"="Primary Blood Derived Cancer - Peripheral Blood",
	"04"="Recurrent Blood Derived Cancer - Bone Marrow",
	"05"="Additional - New Primary",
	"06"="Metastatic",
	"07"="Additional Metastatic",
	"08"="Human Tumor Original Cells",
	"09"="Primary Blood Derived Cancer - Bone Marrow",
	"10"="Blood Derived Normal",
	"11"="Solid Tissue Normal",
	"12"="Buccal Cell Normal",
	"13"="EBV Immortalized Normal",
	"14"="Bone Marrow Normal",
	"15"="sample type 15",
	"16"="sample type 16",
	"20"="Control Analyte",
	"40"="Recurrent Blood Derived Cancer - Peripheral Blood",
	"50"="Cell Lines",
	"60"="Primary Xenograft Tissue",
	"61"="Cell Line Derived Xenograft Tissue",
	"99"="sample type 99")[x]
	l<-unlist(lapply(l,function(y) { if(is.null(y)) { return(NA) } else { return(y) } }),use.names=FALSE)
	return(l)
}

##add GDC download header
addDownloadHeader<-function(x,type=1) {
	if(type==1) {	
		c(paste(c("id",
			"filename",
			"md5",
			"size",
			"state"),
		collapse="\t"),
		x)
	} else {
		c(paste(c("id",
			"filename",
			"md5",
			"size"),
		collapse="\t"),
		x)
	}
}

##add download header
addSampleHeader<-function(x) {
		c(paste(c("File_ID",
			"File_Name",
			"Data_Category",
			"Data_Type",
			"Project_ID",
			"Case_ID",
			"Sample_ID",
			"Sample_Type"),
		collapse="\t"),
		x)
}

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

# c("doSep",
# "doSep2",
# "doUnTab",
# "doUnSpace",
# "doUnPaste",
# "getIdat",
# "idToTcga",
# "sampleDecode",
# "addDownloadHeader",
# "addSampleHeader",
# "parseTcgaName"
# )

################################################################################
##get core expression

m<-read.table(paste0(HOME,"/","processedManifests/","coreSetManifest_gex.txt"),
	header=TRUE,sep="\t",as.is=TRUE
)

files<-addDownloadHeader(doUnTab(unique(m$manifest)))

#get tempdir for file download
temp<-tempdir()

if( Sys.info()["sysname"] == "Windows" ) {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client.exe",MANIFEST),temp)
	writeLines(files,paste0(temp,"/manifest_tmp.txt"))
	setwd(temp)

	} else {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
	writeLines(files,paste0(temp,"/manifest_tmp.txt"))
	setwd(temp)

	#do syscommand to do download
	system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )

	}

# #transfer to tempdir
# file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
# writeLines(files,paste0(temp,"/manifest_tmp.txt"))
# setwd(temp)

# #do syscommand to do download
# system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )

cmd<-paste0(
  paste0(getwd(),"/",list.files(".",pattern="gdc-client")),
  paste0(" download -m  ",
  paste0(getwd(),"/",list.files(".",pattern="manifest_tmp"))
  ))
system(cmd)
rm(cmd)

##get file paths and connect file to sample
paths<-unlist(lapply(m$manifest,function(x) {
	paste(doSep2(x)[c(1,2)],collapse="/")
}))

##counts
paths.counts<-paths[grep("htseq\\.counts\\.gz",paths)]
samples.counts<-m$Sample_ID[grep("htseq\\.counts\\.gz",paths)]

data.counts<-lapply(paths.counts,readLines)

rn<-sub("\\t.+","",data.counts[[1]])
data.counts<-do.call("cbind",lapply(data.counts,function(x) {
	sub(".+\\t","",x)
}))
rownames(data.counts)<-rn ; rm(rn)
colnames(data.counts)<-sub(".+/","",paths.counts)

#remove duplicates - choose one for each sample
keep.counts<-unlist(lapply(split(colnames(data.counts),samples.counts),function(x) x[1]))
data.counts<-data.counts[,keep.counts]
colnames(data.counts)<-names(keep.counts)

##write out kept counts-files
writeLines(paste(names(keep.counts),keep.counts,sep="\t",collapse="\n"),
	con=paste0(HOME,"/","gex/","gex_counts_unique_samples.txt")
	)

rm(paths.counts,samples.counts,keep.counts)


##fpkm
paths.fpkm<-paths[grep("FPKM\\.txt\\.gz",paths)]
samples.fpkm<-m$Sample_ID[grep("FPKM\\.txt\\.gz",paths)]

data.fpkm<-lapply(paths.fpkm,readLines)

rn<-sub("\\t.+","",data.fpkm[[1]])
data.fpkm<-do.call("cbind",lapply(data.fpkm,function(x) {
	sub(".+\\t","",x)
}))
rownames(data.fpkm)<-rn ; rm(rn)
colnames(data.fpkm)<-sub(".+/","",paths.fpkm)

#remove duplicates - choose one
keep.fpkm<-unlist(lapply(split(colnames(data.fpkm),samples.fpkm),function(x) x[1]))
data.fpkm<-data.fpkm[,keep.fpkm]
colnames(data.fpkm)<-names(keep.fpkm)

##write out kept counts-files
writeLines(paste(names(keep.fpkm),keep.fpkm,sep="\t",collapse="\n"),
	con=paste0(HOME,"/","gex/","gex_fpkm_unique_samples.txt")
	)

rm(paths.fpkm,samples.fpkm,keep.fpkm)

##fpkm UQ
paths.uq<-paths[grep("FPKM-UQ\\.txt\\.gz",paths)]
samples.uq<-m$Sample_ID[grep("FPKM-UQ\\.txt\\.gz",paths)]

data.uq<-lapply(paths.uq,readLines)

rn<-sub("\\t.+","",data.uq[[1]])
data.uq<-do.call("cbind",lapply(data.uq,function(x) {
	sub(".+\\t","",x)
}))
rownames(data.uq)<-rn ; rm(rn)
colnames(data.uq)<-sub(".+/","",paths.uq)

#remove duplicates - choose one
keep.uq<-unlist(lapply(split(colnames(data.uq),samples.uq),function(x) x[1]))
data.uq<-data.uq[,keep.uq]
colnames(data.uq)<-names(keep.uq)

##write out kept counts-files
writeLines(paste(names(keep.uq),keep.uq,sep="\t",collapse="\n"),
	con=paste0(HOME,"/","gex/","gex_fpkm-uq_unique_samples.txt")
	)

rm(paths.uq,samples.uq,keep.uq)

#clean
unlink(paste0(getwd(),"/",list.files(".",pattern="gdc-client")))
unlink(paste0(getwd(),"/",list.files(".",pattern="manifest_tmp")))
unlink(sub("/.+","",paths),recursive=TRUE)
rm(m,paths)

ls()
#  [1] "addDownloadHeader" "addSampleHeader"   "data.counts"      
#  [4] "data.fpkm"         "data.uq"           "doSep"            
#  [7] "doSep2"            "doUnPaste"         "doUnSpace"        
# [10] "doUnTab"           "files"             "getIdat"          
# [13] "HOME"              "idToTcga"          "MANIFEST"         
# [16] "parseTcgaName"     "sampleDecode"      "temp"             
# [19] "TUMOR_TYPE"  

str(data.counts)
 # chr [1:60488, 1:1217] "651" "2" "1877" "1113" "387" "250" "914" "1646" ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60488] "ENSG00000000003.13" "ENSG00000000005.5" "ENSG00000000419.11" "ENSG00000000457.12" ...
 #  ..$ : chr [1:1217] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

str(data.fpkm)
 # chr [1:60483, 1:1217] "0.0" "0.0" "3.67636339739" "0.0" "1.82613852501" ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60483] "ENSG00000242268.2" "ENSG00000270112.3" "ENSG00000167578.15" "ENSG00000273842.1" ...
 #  ..$ : chr [1:1217] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

str(data.uq)
 # chr [1:60483, 1:1217] "0.0" "0.0" "76027.9793838" "0.0" "37764.9342908" ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60483] "ENSG00000242268.2" "ENSG00000270112.3" "ENSG00000167578.15" "ENSG00000273842.1" ...
 #  ..$ : chr [1:1217] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

all.equal(colnames(data.counts),colnames(data.fpkm))
#[1] TRUE

all.equal(colnames(data.counts),colnames(data.uq))
#[1] TRUE

#Tumors vs normals?
table( sub(".+-(..).","\\1",colnames(data.counts)) )
#   01   06   11 
# 1097    7  113 

length( sub("-...$","",colnames(data.counts)) )
#[1] 1217

length( unique( sub("-...$","",colnames(data.counts)) ) )
# [1] 1092

#save as RData
save(data.counts,file=paste0(HOME,"/","gex/","gex_counts_unique_samples.RData")	)

save(data.fpkm,file=paste0(HOME,"/","gex/","gex_fpkm_unique_samples.RData")	)

save(data.uq,file=paste0(HOME,"/","gex/","gex_fpkm-uq_unique_samples.RData")	)

################################################################################
##get core copynumbers

m<-read.table(paste0(HOME,"/","processedManifests/","coreSetManifest_ascat.txt"),
	header=TRUE,sep="\t",as.is=TRUE
)

files<-addDownloadHeader(doUnTab(unique(m$manifest)))

#get tempdir for file download
temp<-tempdir()

if( Sys.info()["sysname"] == "Windows" ) {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client.exe",MANIFEST),temp)
	writeLines(files,paste0(temp,"/manifest_tmp.txt"))
	setwd(temp)

	} else {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
	writeLines(files,paste0(temp,"/manifest_tmp.txt"))
	setwd(temp)

	#do syscommand to do download
	system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )

	}

# #transfer to tempdir
# file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
# writeLines(files,paste0(temp,"/manifest_tmp.txt"))
# setwd(temp)

# #do syscommand to do download
# system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )

cmd<-paste0(
  paste0(getwd(),"/",list.files(".",pattern="gdc-client")),
  paste0(" download -m  ",
  paste0(getwd(),"/",list.files(".",pattern="manifest_tmp"))
  ))
system(cmd)
rm(cmd)

##get file paths and connect file to sample
paths<-unlist(lapply(m$manifest,function(x) {
	paste(doSep2(x)[c(1,2)],collapse="/")
}))

table( sub(".+\\.(gene.+)|.+\\.(ascat.+)","\\1\\2",paths) )
# ascat2.allelic_specific.seg.txt      gene_level_copy_number.tsv 
#                            1084                            1084 

##gene-level
paths.cn<-paths[grep(".+\\.(gene.+)",paths)]
samples.cn<-m$Sample_ID[grep(".+\\.(gene.+)",paths)]

head(samples.cn)
# [1] "TCGA-AO-A03L-01A, TCGA-AO-A03L-10A" "TCGA-E2-A2P5-10A, TCGA-E2-A2P5-01A"
# [3] "TCGA-C8-A273-10A, TCGA-C8-A273-01A" "TCGA-AC-A6IX-10A, TCGA-AC-A6IX-01A"
# [5] "TCGA-D8-A146-01A, TCGA-D8-A146-10A" "TCGA-A8-A09R-01A, TCGA-A8-A09R-10A"

##each sample is tumor-normal pair
table( gsub("TCGA-..-....-(..).,.+TCGA-..-....-(..).","\\1_\\2",samples.cn) )
# 01_10 01_11 06_10 06_11 10_01 10_06 11_01 
#   487    45     2     3   504     2    41 

table( sub(".+-","",sub(".*(TCGA-..-....-0..).*","\\1",samples.cn)) )
#  01A  01B  06A 
# 1058   19    7 

##separate tumor id from paired normal id
samples.cn<-sub(".*(TCGA-..-....-0..).*","\\1",samples.cn)

length(samples.cn)-length(unique(samples.cn))
#[1] 2

##read data
data.cn<-lapply(paths.cn,readLines)

table( unlist(lapply(data.cn,length)) )
# 60447 
#  1084 

x<-sub("\\t.+","",data.cn[[1]])
table( unlist(lapply(data.cn,function(y) { y<-sub("\\t.+","",y) ; all(y == x) })) )
# TRUE 
# 1084 
rm(x) 

#cut out data
rn<-sub("\\t.+","",data.cn[[1]][-1])
data.cn<-do.call("cbind",lapply(data.cn,function(x) {
	sub("ENSG.+\\t.+\\tchr.+\\t\\d+\\t\\d+(\\t.*\\t.*\\t.*)","\\1\\2\\3",x[-1])
}))
rownames(data.cn)<-rn ; rm(rn)
colnames(data.cn)<-sub(".+/","",paths.cn)

#remove duplicates - choose one for each sample
keep.cn<-unlist(lapply(split(colnames(data.cn),samples.cn),function(x) x[1]))
data.cn<-data.cn[,keep.cn]
colnames(data.cn)<-names(keep.cn)

str(data.cn)
 # chr [1:60446, 1:1082] "\t\t\t" "\t\t\t" "\t\t\t" "\t\t\t" "\t\t\t" ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60446] "ENSG00000223972.5" "ENSG00000227232.5" "ENSG00000278267.1" "ENSG00000243485.3" ...
 #  ..$ : chr [1:1082] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

##for each sample gene-copy is three values = cn \t min cn \t max cn
x<-(apply(data.cn,2,function(q) { sub("\\t(.*)\\t(.*)\\t(.*)","\\1",q) }))
y<-(apply(data.cn,2,function(q) { sub("\\t(.*)\\t(.*)\\t(.*)","\\2",q) }))
z<-(apply(data.cn,2,function(q) { sub("\\t(.*)\\t(.*)\\t(.*)","\\3",q) }))

sum(x == y) / length(x)
#[1] 0.9993181
sum(x == z) / length(x)
# [1] 0.9993435
sum(y == z) / length(x)
#[1] 0.9986872

##only save actual cn.
all.equal(rownames(data.cn),rownames(x))
#[1] TRUE
all.equal(colnames(data.cn),colnames(x))
#[1] TRUE

data.cn<-x
rm(x,y,z)

##write out kept counts-files
writeLines(paste(names(keep.cn),keep.cn,sep="\t",collapse="\n"),
	con=paste0(HOME,"/","cn/","cn_geneCopy_unique_samples.txt")
	)

rm(paths.cn,samples.cn,keep.cn)

##allele-level
paths.seg<-paths[grep(".+\\.(ascat.+)",paths)]
samples.seg<-m$Sample_ID[grep(".+\\.(ascat.+)",paths)]

head(samples.seg)
# [1] "TCGA-A8-A07F-01A, TCGA-A8-A07F-10A" "TCGA-E9-A1R7-10A, TCGA-E9-A1R7-01A"
# [3] "TCGA-B6-A0IC-01A, TCGA-B6-A0IC-10A" "TCGA-D8-A1JE-01A, TCGA-D8-A1JE-10A"
# [5] "TCGA-A7-A26E-10A, TCGA-A7-A26E-01A" "TCGA-B6-A0I9-01A, TCGA-B6-A0I9-10A"

##each sample is tumor-normal pair
table( gsub("TCGA-..-....-(..).,.+TCGA-..-....-(..).","\\1_\\2",samples.seg) )
# 01_10 01_11 06_10 06_11 10_01 10_06 11_01 
#   487    45     2     3   504     2    41 

table( sub(".+-","",sub(".*(TCGA-..-....-0..).*","\\1",samples.seg)) )
#  01A  01B  06A 
# 1058   19    7 

##separate tumor id from paired normal id
samples.seg<-sub(".*(TCGA-..-....-0..).*","\\1",samples.seg)

length(samples.seg)-length(unique(samples.seg))
#[1] 2

head(readLines(paths.seg[1]),2)
# [1] "GDC_Aliquot\tChromosome\tStart\tEnd\tCopy_Number\tMajor_Copy_Number\tMinor_Copy_Number"
# [2] "1609caae-f341-4227-ac2a-a76a6c253e95\tchr1\t62920\t61395882\t2\t1\t1"   

##read data
data.seg<-lapply(paths.seg,readLines)

quantile( unlist(lapply(data.seg,length)) )
 #   0%   25%   50%   75%  100% 
 # 49.0 117.0 174.5 266.0 709.0 

##trim data - only keep essential data

#cut out data
cn<-unlist(strsplit(data.seg[[1]][1],"\\t"))
x<-do.call("rbind",lapply(data.seg,function(x) {
	do.call("rbind",lapply(x[-1],function(y) unlist(strsplit(y,"\\t")) ))
}))
colnames(x)<-cn ; rm(cn)

##remap samples
y<-sub(".+/TCGA-.+\\.(.+)\\.ascat2\\..+","\\1",paths.seg)
z<-samples.seg[match(x[,1],y)]

#remove duplicates - choose one for each sample - change identifiers
keep.seg<-unlist(lapply(split(x[,1],z),function(q) unique(q)[1]))

table(x[,1] %in% keep.seg)
 # FALSE   TRUE 
 #   292 223227 

x<-x[x[,1] %in% keep.seg,]
z<-samples.seg[match(x[,1],y)]

x[,1]<-z
rm(y,z)

data.seg<-x

str(data.seg)
 # chr [1:223227, 1:7] "TCGA-A8-A07F-01A" "TCGA-A8-A07F-01A" ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : NULL
 #  ..$ : chr [1:7] "GDC_Aliquot" "Chromosome" "Start" "End" ...

colnames(data.seg)[1]<-"Sample_ID"

##segment LOH-estimates per gene can be derived in later scripts..

##write out kept seg-files
writeLines(paste(names(keep.seg),keep.seg,sep="\t",collapse="\n"),
	con=paste0(HOME,"/","cn/","seg_geneCopy_unique_samples.txt")
	)

rm(paths.seg,samples.seg,keep.seg)

#clean
unlink(paste0(getwd(),"/",list.files(".",pattern="gdc-client")))
unlink(paste0(getwd(),"/",list.files(".",pattern="manifest_tmp")))
unlink(sub("/.+","",paths),recursive=TRUE)
rm(m,paths,x)

#save as RData
save(data.seg,file=paste0(HOME,"/","cn/","cn_seg_unique_samples.RData")	)

save(data.cn,file=paste0(HOME,"/","cn/","cn_gene_unique_samples.RData")	)

################################################################################
##get core methylation

m<-read.table(paste0(HOME,"/","processedManifests/","coreSetManifest_me.txt"),
	header=TRUE,sep="\t",as.is=TRUE
)

files<-addDownloadHeader(doUnTab(unique(m$manifest)))

length(files)
#[1] 1791

#get tempdir for file download
temp<-tempdir()

if( Sys.info()["sysname"] == "Windows" ) {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client.exe",MANIFEST),temp)
	writeLines(files,paste0(temp,"/manifest_tmp.txt"))
	setwd(temp)

	} else {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
	writeLines(files,paste0(temp,"/manifest_tmp.txt"))
	setwd(temp)

	#do syscommand to do download
	system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )
}

# #transfer to tempdir
# file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
# writeLines(files,paste0(temp,"/manifest_tmp.txt"))
# setwd(temp)

# #do syscommand to do download
# system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )

cmd<-paste0(
  paste0(getwd(),"/",list.files(".",pattern="gdc-client")),
  paste0(" download -m  ",
  paste0(getwd(),"/",list.files(".",pattern="manifest_tmp"))
  ))
system(cmd)
rm(cmd)

##get file paths and connect file to sample
paths<-unlist(lapply(m$manifest,function(x) {
	paste(doSep2(x)[c(1,2)],collapse="/")
}))

table( sub(".+/.+_([RG]..)\\.idat","\\1",paths) )
# Grn Red 
# 895 895 

length(unique(m$Sample_ID))
#[1] 893

length(unique(sub("_[RG]..\\.idat","",m$File_Name)))
#[1] 895

length(unique(sub("_R..C.._[RG]..\\.idat","",m$File_Name)))
#[1] 98

##read in as minfi-object for later normalization -> needs a sample-sheet
dir.create("arrays")
for (i in unique(sub("_R..C.._[RG]..\\.idat","",m$File_Name))) {
	dir.create(paste0("arrays/",i))
} ; rm(i)

for (i in paths) {
	file.copy( i , paste0("arrays/",
		sub(".+/(.+)_R..C.._[RG]..\\.idat","\\1/",i),
		sub(".+/(.+_R..C.._[RG]..\\.idat)","\\1",i)
		) 
	)
} ; rm(i)

length(unlist(lapply(dir("arrays",full.names=TRUE),function(x) list.files(x,full.names=TRUE))))
#[1] 1790

head(unlist(lapply(dir("arrays",full.names=TRUE),function(x) list.files(x,full.names=TRUE))))
# [1] "arrays/3999997079/3999997079_R01C01_Grn.idat"
# [2] "arrays/3999997079/3999997079_R01C01_Red.idat"
# [3] "arrays/3999997079/3999997079_R01C02_Grn.idat"
# [4] "arrays/3999997079/3999997079_R01C02_Red.idat"
# [5] "arrays/3999997079/3999997079_R02C01_Grn.idat"
# [6] "arrays/3999997079/3999997079_R02C01_Red.idat"

bname<-sub(".+/","",unlist(lapply(dir("arrays",full.names=TRUE),function(x) list.files(x,full.names=TRUE))))
bname2<-unlist(lapply(dir("arrays",full.names=TRUE),function(x) list.files(x,full.names=TRUE)))

bname<-bname[grep("Grn\\.idat",bname)]
bname2<-bname2[grep("Grn\\.idat",bname2)]
bname2<-sub("_Grn\\.idat","",bname2)

length(bname)
#[1] 895
length(bname2)
#[1] 895

##only grab annotation for green (same for red) - make into minfi sample sheet
sampleSheet<-as.data.frame(cbind(Sample_ID=m$Sample_ID[match(bname,m$File_Name)],
	Sample_Type=m$Sample_Type[match(bname,m$File_Name)],
	Slide=sub("_.+","",bname),
	Array=sub(".+_(.+)_[GR].+","\\1",bname),
	Basename=paste0(getwd(),"/",bname2)
	),stringsAsFactors=FALSE)

str(sampleSheet)
# 'data.frame':	895 obs. of  5 variables:
#  $ Sample_ID  : chr  "TCGA-PL-A8LX-01A" "TCGA-Z7-A8R6-01A" "TCGA-PL-A8LY-01A" "TCGA-S3-AA17-01A" ...
#  $ Sample_Type: chr  "Primary Solid Tumor" "Primary Solid Tumor" "Primary Solid Tumor" "Primary Solid Tumor" ...
#  $ Slide      : chr  "3999997079" "3999997079" "3999997079" "3999997079" ...
#  $ Array      : chr  "R01C01" "R01C02" "R02C01" "R03C01" ...
#  $ Basename   : chr  "/tmp/Rtmp505mX4/arrays/3999997079/3999997079_R01C01" "/tmp/Rtmp505mX4/arrays/3999997079/3999997079_R01C02" "/tmp/Rtmp505mX4/arrays/3999997079/3999997079_R02C01" "/tmp/Rtmp505mX4/arrays/3999997079/3999997079_R03C01" ...

rm(bname,bname2)

##trim to get rid of duplicate samples
keep.me<-unlist(lapply(split(1:nrow(sampleSheet),sampleSheet$Sample_ID),function(x) x[1]))

length(keep.me)
#[1] 893

sampleSheet<-sampleSheet[keep.me,]

methSet <- read.metharray.exp(targets = sampleSheet)

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

##write out kept me-files
writeLines(paste(sampleSheet$Sample_ID,paste(sampleSheet$Slide,sampleSheet$Array,sep="_"),sep="\t",collapse="\n"),
	con=paste0(HOME,"/","me/","me_unique_samples.txt")
	)

save(methSet,file=paste0(HOME,"/","me/","me_minfi_RGset_unique_samples.RData")	)

rm(keep.me,methSet,sampleSheet)

#clean
unlink(paste0(getwd(),"/",list.files(".",pattern="gdc-client")))
unlink(paste0(getwd(),"/",list.files(".",pattern="manifest_tmp")))
unlink(sub("/.+","",paths),recursive=TRUE)
unlink("arrays",recursive=TRUE)
rm(m,paths,files)

gc()

################################################################################
##get wes

##download-manifest repeated for each sample and wes-pipeline but in reality downloading 4 unique files is enough
m<-read.table(paste0(HOME,"/","processedManifests/","coreSetManifest_wes.txt"),
	header=TRUE,sep="\t",as.is=TRUE
)

files<-addDownloadHeader(doUnTab(unique(m$manifest)))

files
# [1] "id\tfilename\tmd5\tsize\tstate"                                                                                                                                                 
# [2] "b8ca5856-9819-459c-87c5-94e91aca4032\tTCGA.BRCA.muse.b8ca5856-9819-459c-87c5-94e91aca4032.DR-10.0.somatic.maf.gz\tf814c909403da10f3d2200d2a4b45f42\t23039649\treleased"         
# [3] "995c0111-d90b-4140-bee7-3845436c3b42\tTCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz\t4a00a1c3611c2ddcb2eed684c451eea7\t31058116\treleased"       
# [4] "7dd592e3-5950-4438-96d5-3c718aca3f13\tTCGA.BRCA.somaticsniper.7dd592e3-5950-4438-96d5-3c718aca3f13.DR-10.0.somatic.maf.gz\t21d64325c5780be0b1fadcc55fe5e3d4\t16910889\treleased"
# [5] "6c93f518-1956-4435-9806-37185266d248\tTCGA.BRCA.varscan.6c93f518-1956-4435-9806-37185266d248.DR-10.0.somatic.maf.gz\t23f303806e40f6f30fc52b5a9c6eb327\t24161125\treleased"      

#get tempdir for file download
temp<-tempdir()

if( Sys.info()["sysname"] == "Windows" ) {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client.exe",MANIFEST),temp)
	writeLines(files,paste0(temp,"/manifest_tmp.txt"))
	setwd(temp)

	} else {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
	writeLines(files,paste0(temp,"/manifest_tmp.txt"))
	setwd(temp)

	#do syscommand to do download
	system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )
}

# #transfer to tempdir
# file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
# writeLines(files,paste0(temp,"/manifest_tmp.txt"))
# setwd(temp)

# #do syscommand to do download
# system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )

cmd<-paste0(
  paste0(getwd(),"/",list.files(".",pattern="gdc-client")),
  paste0(" download -m  ",
  paste0(getwd(),"/",list.files(".",pattern="manifest_tmp"))
  ))
system(cmd)
rm(cmd)

##get file paths and connect file to sample
paths<-unlist(lapply(unique(m$manifest),function(x) {
	paste(doSep2(x)[c(1,2)],collapse="/")
}))

data.wes<-lapply(paths,function(x) {
	readLines(x)
})

str(data.wes)
# List of 4
#  $ : chr [1:90975] "#version gdc-1.0.0" "#filedate 20170929" "#annotation.spec gdc-1.0.1-public" "#n.analyzed.samples 986" ...
#  $ : chr [1:120994] "#version gdc-1.0.0" "#filedate 20170930" "#annotation.spec gdc-1.0.1-public" "#n.analyzed.samples 986" ...
#  $ : chr [1:64898] "#version gdc-1.0.0" "#filedate 20170929" "#annotation.spec gdc-1.0.1-public" "#n.analyzed.samples 986" ...
#  $ : chr [1:93618] "#version gdc-1.0.0" "#filedate 20170929" "#annotation.spec gdc-1.0.1-public" "#n.analyzed.samples 986" ...

##grep 4 data types
data.muse<-data.wes[[ grep("muse",paths) ]]

data.mutect<-data.wes[[ grep("mutect",paths) ]]

data.somsniper<-data.wes[[ grep("somaticsniper",paths) ]]

data.varscan<-data.wes[[ grep("varscan",paths) ]]

##separate on tab
data.muse<-lapply(data.muse,doSep)

data.mutect<-lapply(data.mutect,doSep)

data.somsniper<-lapply(data.somsniper,doSep)

data.varscan<-lapply(data.varscan,doSep)

##check 
table( unlist(lapply(data.muse,length)) )
    # 1   120 
    # 5 90970 

table( unlist(lapply(data.mutect,length)) )
     # 1    120 
     # 5 120989 

table( unlist(lapply(data.somsniper,length)) )
    # 1   120 
    # 5 64893 

table( unlist(lapply(data.varscan,length)) )
    # 1   120 
    # 5 93613 

##remove extra rows
data.muse<-data.muse[ unlist(lapply(data.muse,length)) == 120 ]

data.mutect<-data.mutect[ unlist(lapply(data.mutect,length)) == 120 ]

data.somsniper<-data.somsniper[ unlist(lapply(data.somsniper,length)) == 120 ]

data.varscan<-data.varscan[ unlist(lapply(data.varscan,length)) == 120 ]

##merge into matrix
data.muse<-do.call("rbind",data.muse)
cols<-data.muse[1,]
colnames(data.muse)<-cols
data.muse<-data.muse[-1,]
data.muse<-as.data.frame(data.muse,stringsAsFactors=FALSE)

data.mutect<-do.call("rbind",data.mutect)
cols<-data.mutect[1,]
colnames(data.mutect)<-cols
data.mutect<-data.mutect[-1,]
data.mutect<-as.data.frame(data.mutect,stringsAsFactors=FALSE)

data.somsniper<-do.call("rbind",data.somsniper)
cols<-data.somsniper[1,]
colnames(data.somsniper)<-cols
data.somsniper<-data.somsniper[-1,]
data.somsniper<-as.data.frame(data.somsniper,stringsAsFactors=FALSE)

data.varscan<-do.call("rbind",data.varscan)
cols<-data.varscan[1,]
colnames(data.varscan)<-cols
data.varscan<-data.varscan[-1,]
data.varscan<-as.data.frame(data.varscan,stringsAsFactors=FALSE)

##check variants
all.equal(colnames(data.muse),colnames(data.mutect))
#[1] TRUE

all.equal(colnames(data.muse),colnames(data.somsniper))
#[1] TRUE

all.equal(colnames(data.muse),colnames(data.varscan))
#[1] TRUE

table(data.muse$Variant_Type)
#   SNP 
# 90969 
table(data.mutect$Variant_Type)
  #  DEL    INS    SNP 
  # 5209   3739 112040 
table(data.somsniper$Variant_Type)
#   SNP 
# 64892 
table(data.varscan$Variant_Type)
 #  DEL   INS   SNP 
 # 5243  1590 86779 

table( sub("\\(.+","",data.muse$SIFT) )
#                                           deleterious 
#                      43868                      24305 
# deleterious_low_confidence                  tolerated 
#                       2092                      18921 
#   tolerated_low_confidence 
#                       1783 

table( sub("\\(.+","",data.muse$PolyPhen) )
          #                    benign possibly_damaging probably_damaging 
          #   40331             21658              9821             17672 
          # unknown 
          #    1487 

table(sub("\\(.+","",data.muse$SIFT),sub("\\(.+","",data.muse$PolyPhen))
  #                                  benign possibly_damaging probably_damaging
  #                            40331   1174               638              1119
  # deleterious                    0   4732              5569             13682
  # deleterious_low_confidence     0    731               530               659
  # tolerated                      0  13628              2921              2125
  # tolerated_low_confidence       0   1393               163                87
                            
  #                            unknown
  #                                606
  # deleterious                    322
  # deleterious_low_confidence     172
  # tolerated                      247
  # tolerated_low_confidence       140

table(data.muse$One_Consequence,data.muse$IMPACT)                                    
  #                                      HIGH   LOW MODERATE MODIFIER
  # 3_prime_UTR_variant                    0     0        0     6012
  # 5_prime_UTR_variant                    0     0        0     2198
  # coding_sequence_variant                0     0        0        2
  # downstream_gene_variant                0     0        0      609
  # incomplete_terminal_codon_variant      0     2        0        0
  # intergenic_variant                     0     0        0        3
  # intron_variant                         0     0        0     4182
  # mature_miRNA_variant                   0     0        0       16
  # missense_variant                       0     0    50806        0
  # non_coding_transcript_exon_variant     0     0        0     1461
  # regulatory_region_variant              0     0        0        2
  # splice_acceptor_variant              799     0        0        0
  # splice_donor_variant                 528     0        0        0
  # splice_region_variant                  0  1193        0        0
  # start_lost                            62     0        0        0
  # stop_gained                         4480     0        0        0
  # stop_lost                             67     0        0        0
  # stop_retained_variant                  0    59        0        0
  # synonymous_variant                     0 17932        0        0
  # upstream_gene_variant                  0     0        0      556
 
table(sub("\\(.+","",data.muse$SIFT),data.muse$IMPACT)                                    
  #                            HIGH   LOW MODERATE MODIFIER
  #                             5876 19186     3765    15041
  # deleterious                   24     0    24281        0
  # deleterious_low_confidence    26     0     2066        0
  # tolerated                      6     0    18915        0
  # tolerated_low_confidence       4     0     1779        0

table(sub("\\(.+","",data.muse$PolyPhen),data.muse$IMPACT)                                    
  #                    HIGH   LOW MODERATE MODIFIER
  #                    5874 19186      230    15041
  # benign               25     0    21633        0
  # possibly_damaging    13     0     9808        0
  # probably_damaging    19     0    17653        0
  # unknown               5     0     1482        0

##do consensus of methods and keep moderate and high impact variants.
muse.combn<-apply(data.muse,1,function(x) {
	paste(x[c("Tumor_Sample_Barcode",
		"Hugo_Symbol",
		"Entrez_Gene_Id",
		"Gene",
		"Chromosome",
		"Start_Position",
		"End_Position",
		"Strand",
		"Reference_Allele",
		"Tumor_Seq_Allele1",
		"Tumor_Seq_Allele2",
		"IMPACT",
		"Variant_Classification",
		"Variant_Type")],collapse="|")
})

mutect.combn<-apply(data.mutect,1,function(x) {
	paste(x[c("Tumor_Sample_Barcode",
		"Hugo_Symbol",
		"Entrez_Gene_Id",
		"Gene",
		"Chromosome",
		"Start_Position",
		"End_Position",
		"Strand",
		"Reference_Allele",
		"Tumor_Seq_Allele1",
		"Tumor_Seq_Allele2",
		"IMPACT",		
		"Variant_Classification",
		"Variant_Type")],collapse="|")
})

somsniper.combn<-apply(data.somsniper,1,function(x) {
	paste(x[c("Tumor_Sample_Barcode",
		"Hugo_Symbol",
		"Entrez_Gene_Id",
		"Gene",
		"Chromosome",
		"Start_Position",
		"End_Position",
		"Strand",
		"Reference_Allele",
		"Tumor_Seq_Allele1",
		"Tumor_Seq_Allele2",
		"IMPACT",		
		"Variant_Classification",
		"Variant_Type")],collapse="|")
})

varscan.combn<-apply(data.varscan,1,function(x) {
	paste(x[c("Tumor_Sample_Barcode",
		"Hugo_Symbol",
		"Entrez_Gene_Id",
		"Gene",
		"Chromosome",
		"Start_Position",
		"End_Position",
		"Strand",
		"Reference_Allele",
		"Tumor_Seq_Allele1",
		"Tumor_Seq_Allele2",
		"IMPACT",		
		"Variant_Classification",
		"Variant_Type")],collapse="|")
})

all<-unique(c(muse.combn,mutect.combn,somsniper.combn,varscan.combn))

head(all)
# [1] "TCGA-BH-A0HO-01A-11W-A050-09|FOXD3|27022|ENSG00000187140|chr1|63323535|63323535|+|G|G|A|LOW|Silent|SNP"                   
# [2] "TCGA-BH-A0HO-01A-11W-A050-09|KPRP|448834|ENSG00000203786|chr1|152760163|152760163|+|G|G|A|MODERATE|Missense_Mutation|SNP" 
# [3] "TCGA-BH-A0HO-01A-11W-A050-09|SFT2D2|375035|ENSG00000213064|chr1|168246650|168246650|+|T|T|G|MODIFIER|3'UTR|SNP"           
# [4] "TCGA-BH-A0HO-01A-11W-A050-09|KIAA1407|57577|ENSG00000163617|chr3|113965097|113965097|+|A|A|C|MODIFIER|3'UTR|SNP"          
# [5] "TCGA-BH-A0HO-01A-11W-A050-09|FNDC3B|64778|ENSG00000075420|chr3|172133474|172133474|+|A|A|G|MODERATE|Missense_Mutation|SNP"
# [6] "TCGA-BH-A0HO-01A-11W-A050-09|RAI14|26064|ENSG00000039560|chr5|34812204|34812204|+|A|A|C|MODERATE|Missense_Mutation|SNP"   

length(unique(all))
#[1] 134066

all.combn<-cbind(muse=as.integer(!is.na(match(all,muse.combn))),
mutect=as.integer(!is.na(match(all,mutect.combn))),
somsniper=as.integer(!is.na(match(all,somsniper.combn))),
varscan=as.integer(!is.na(match(all,varscan.combn))))

table(rowSums(all.combn))
#     1     2     3     4 
# 34178 20468 22333 57087 

all.map<-cbind(muse=match(all,muse.combn),
mutect=match(all,mutect.combn),
somsniper=match(all,somsniper.combn),
varscan=match(all,varscan.combn))

table(sub(".+\\|","",all))
  #  DEL    INS    SNP 
  # 7356   4295 122415 

table(sub(".+\\|","",all[rowSums(all.combn)>1]))
 #  DEL   INS   SNP 
 # 3096  1034 95758 

##remove mutations with less than 2 hits
all<-all[rowSums(all.combn)>1]

all.combn<-cbind(muse=as.integer(!is.na(match(all,muse.combn))),
mutect=as.integer(!is.na(match(all,mutect.combn))),
somsniper=as.integer(!is.na(match(all,somsniper.combn))),
varscan=as.integer(!is.na(match(all,varscan.combn))))

table(rowSums(all.combn))
#     2     3     4 
# 20468 22333 57087 

all.map<-cbind(muse=match(all,muse.combn),
mutect=match(all,mutect.combn),
somsniper=match(all,somsniper.combn),
varscan=match(all,varscan.combn))

table(sub(".+\\|","",all))
 #  DEL   INS   SNP 
 # 3096  1034 95758 

mutations.consensus<-do.call("rbind",lapply(all,doSep2))
colnames(mutations.consensus)<-c("Tumor_Sample_Barcode",
		"Hugo_Symbol",
		"Entrez_Gene_Id",
		"Gene",
		"Chromosome",
		"Start_Position",
		"End_Position",
		"Strand",
		"Reference_Allele",
		"Tumor_Seq_Allele1",
		"Tumor_Seq_Allele2",
		"IMPACT",		
		"Variant_Classification",
		"Variant_Type")
mutations.consensus<-as.data.frame(mutations.consensus,stringsAsFactors=FALSE)

table(mutations.consensus$Variant_Classification,mutations.consensus$IMPACT)
  #                         HIGH   LOW MODERATE MODIFIER
  # 3'Flank                    0     0        0      646
  # 3'UTR                      0     0        0     6326
  # 5'Flank                    0     0        0      592
  # 5'UTR                      0     0        0     2314
  # Frame_Shift_Del         1986     0        0        0
  # Frame_Shift_Ins          768     0        0        0
  # IGR                        0     0        0        5
  # In_Frame_Del               0     0      334        0
  # In_Frame_Ins               0     0       42        0
  # Intron                     0     0        0     4463
  # Missense_Mutation          0     0    53798        5
  # Nonsense_Mutation       4701     0        0        0
  # Nonstop_Mutation          72     0        0        0
  # RNA                        0     0        0     1571
  # Silent                     0 19397        0        0
  # Splice_Region              0  1284        2        0
  # Splice_Site             1521     0        0        0
  # Translation_Start_Site    61     0        0        0

dim(all.map)
#[1] 99888     4
dim(mutations.consensus)
#[1] 99888    14

##map to original row of SNV-files
mutations.consensus$muse_row<-all.map[,"muse"]
mutations.consensus$mutect_row<-all.map[,"mutect"]
mutations.consensus$somsniper_row<-all.map[,"somsniper"]
mutations.consensus$varscan_row<-all.map[,"varscan"]

unique(mutations.consensus$Tumor_Sample_Barcode)

##write out kept me-files
writeLines(paste(parseTcgaName(unique(mutations.consensus$Tumor_Sample_Barcode),level=4),"Muse, Mutect, SomaticSniper, Varscan",sep="\t",collapse="\n"),
	con=paste0(HOME,"/","wes/","wes_unique_samples.txt")
	)

save(mutations.consensus,file=paste0(HOME,"/","wes/","wes_mutations_consensus_unique_samples.RData")	)

save(data.muse,file=paste0(HOME,"/","wes/","wes_mutations_muse_unique_samples.RData")	)
save(data.mutect,file=paste0(HOME,"/","wes/","wes_mutations_mutect_unique_samples.RData")	)
save(data.somsniper,file=paste0(HOME,"/","wes/","wes_mutations_somaticsniper_unique_samples.RData")	)
save(data.varscan,file=paste0(HOME,"/","wes/","wes_mutations_varscan_unique_samples.RData")	)

#clean
rm(all,all.combn,all.map,cols,data.wes,muse.combn,mutect.combn,somsniper.combn,varscan.combn)

unlink(paste0(getwd(),"/",list.files(".",pattern="gdc-client")))
unlink(paste0(getwd(),"/",list.files(".",pattern="manifest_tmp")))
unlink(sub("/.+","",paths),recursive=TRUE)
unlink("arrays",recursive=TRUE)
rm(m,paths,files,temp)

gc()

################################################################################
##get atac - manifest for "extended set"

m<-read.table(paste0(HOME,"/","processedManifests/","extendedSetManifest_atac.txt"),
	header=TRUE,sep="\t",as.is=TRUE
)

key<-read.table(paste0(HOME,"/","processedManifests/","stanfordToTcga_sampleMapping.txt"),
	header=FALSE,sep="\t",as.is=TRUE,col.names=c("stanford_id","tcga_id")
)

str(key)
# 'data.frame':	141 obs. of  2 variables:
#  $ stanford_id: chr  "BRCA-000CFD9F-ADDF-4304-9E60-6041549E189C-X017-S06-L011-B1-T1-P040" "BRCA-000CFD9F-ADDF-4304-9E60-6041549E189C-X017-S06-L012-B1-T2-P046" "BRCA-01112370-4F6F-4A20-9BE0-7975C3465268-X017-S04-L007-B1-T1-P042" "BRCA-01112370-4F6F-4A20-9BE0-7975C3465268-X017-S04-L008-B1-T2-P044" ...
#  $ tcga_id    : chr  "TCGA-A7-A13F-01A" "TCGA-A7-A13F-01A" "TCGA-A7-A13E-01A" "TCGA-A7-A13E-01A" ...

key$stanford_id<-gsub("-","_",key$stanford_id)
key<-key[order(key$tcga_id),]

##only one file for log2 normalized ATAC counts
files<-addDownloadHeader(doUnTab(unique(m$manifest)))

files
# [1] "id\tfilename\tmd5\tsize\tstate"                                                                                                                      
# [2] "38b8f311-f3a4-4746-9829-b8e3edb9c157\tTCGA-ATAC_Cancer_Type-specific_Count_Matrices_log2norm_counts.zip\tc31893739baa05a95dc57f6ab38150b4\t632769396"

#get tempdir for file download
temp<-tempdir()

if( Sys.info()["sysname"] == "Windows" ) {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client.exe",MANIFEST),temp)
	writeLines(files,paste0(temp,"/manifest_tmp.txt"))
	setwd(temp)

	} else {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
	writeLines(files,paste0(temp,"/manifest_tmp.txt"))
	setwd(temp)

	#do syscommand to do download
	system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )
}

# #transfer to tempdir
# file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
# writeLines(files,paste0(temp,"/manifest_tmp.txt"))
# setwd(temp)

##do syscommand to do download
#system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )

cmd<-paste0(
  paste0(getwd(),"/",list.files(".",pattern="gdc-client")),
  paste0(" download -m  ",
  paste0(getwd(),"/",list.files(".",pattern="manifest_tmp"))
  ))
system(cmd)
rm(cmd)

##get file paths and connect file to sample
paths<-unlist(lapply(unique(m$manifest),function(x) {
	paste(doSep2(x)[c(1,2)],collapse="/")
}))

data.atac<-do.call("rbind",lapply(paths,function(x) {
	unzip(zipfile=x,list=TRUE)
}))

data.atac<-data.atac[grep(TUMOR_TYPE,data.atac[,1],ignore.case=TRUE),1]

unzip(zipfile=paths,files=data.atac)

data.atac<-readLines(con=data.atac)

data.atac<-lapply(data.atac,doSep)

table(unlist(lapply(data.atac,length)))
#    146 
# 215921 

data.atac<-do.call("rbind",data.atac)

colnames(data.atac)<-data.atac[1,]
data.atac<-data.atac[-1,]

##separate annotations from data
anno.atac<-data.atac[,1:5]
data.atac<-data.atac[,-1:-5]
mode(data.atac)<-"numeric"

all(colnames(data.atac) %in% key$stanford_id)
#[1] TRUE
all(key$stanford_id %in% colnames(data.atac))
#[1] TRUE

data.atac<-data.atac[,key$stanford_id]

data.atac<-t(apply(data.atac,1,function(x) {
	unlist(lapply(split(x,key$tcga_id),mean))
}))

dim(data.atac)
#[1] 215920     74

##add rownames - drop "score" column
anno.atac<-apply(anno.atac,1,function(x) { 
	paste(paste(x[1],paste(x[2:3],collapse="-"),sep=":"),x[4],sep="|")
})

rownames(data.atac)<-anno.atac
rm(anno.atac)

##500bp consensus format..
quantile(as.integer(sub(".+:.+\\-(.+)\\|.+","\\1",rownames(data.atac)))-
	as.integer(sub(".+:(.+)\\-.+\\|.+","\\1",rownames(data.atac)))
 )
 #  0%  25%  50%  75% 100% 
 # 499  499  499  499  499 

save(data.atac,file=paste0(HOME,"/","atac/","atac_avgLog2normPeaks_unique_samples.RData"))

unlink(paste0(getwd(),"/",list.files(".",pattern="gdc-client")))
unlink(paste0(getwd(),"/",list.files(".",pattern="manifest_tmp")))
unlink(sub("/.+","",paths),recursive=TRUE)
unlink("arrays",recursive=TRUE)
rm(m,key,paths,files,temp)

rm(list=c("doSep",
"doSep2",
"doUnTab",
"doUnSpace",
"doUnPaste",
"getIdat",
"idToTcga",
"sampleDecode",
"addDownloadHeader",
"addSampleHeader",
"parseTcgaName"
))

ls()
#  [1] "data.atac"           "data.cn"             "data.counts"        
#  [4] "data.fpkm"           "data.muse"           "data.mutect"        
#  [7] "data.seg"            "data.somsniper"      "data.uq"            
# [10] "data.varscan"        "HOME"                "MANIFEST"           
# [13] "mutations.consensus" "TUMOR_TYPE"         

##these not needed
rm(data.muse,data.mutect,data.somsniper,data.varscan)

ls(pattern="data|mutations")
# [1] "data.atac"           "data.cn"             "data.counts"        
# [4] "data.fpkm"           "data.seg"            "data.uq"            
# [7] "mutations.consensus"

save(list=ls(pattern="data|mutations"),file=paste0(HOME,"/","coreData_gexCnWesAtac_unique_samples.RData"))

rm(HOME,MANIFEST,TUMOR_TYPE)

q("no")
##END