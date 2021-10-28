#####======================================================================#####
### Get TCGA/GDC data for TCGA BRCA data
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncoloy and Pathology

################################################################################
##Set home directory

##downloads will require substantial amount of space, e.g. 1-2Tb for TCGA-BRCA

##set/create own home directory below:

##work
dir.create("~/hdd1/tcgaBrca")
HOME<-"~/hdd1/tcgaBrca"
MANIFEST<-"~/Documents/tcgaBrca/manifest"
##home
dir.create("I:/data/tcgaBrca")
HOME<-"I:/data/tcgaBrca"
MANIFEST<-"F:/gitProjects/tcgaBrca/manifest"

##tumor type
TUMOR_TYPE<-"brca"

list.files(MANIFEST,full.names=T)
#[1] "/home/med-mai/Documents/tcgaBrca/manifest/atac"  
#[2] "/home/med-mai/Documents/tcgaBrca/manifest/brca"  
#[5] "/home/med-mai/Documents/tcgaBrca/manifest/pancan"

##manifests by data type for GDC BRCA data files
list.files(paste0(MANIFEST,"/brca"),full.names=T)
#[1] "/home/med-mai/Documents/tcgaBrca/manifest/brca/image"      
#[2] "/home/med-mai/Documents/tcgaBrca/manifest/brca/me450"      
#[3] "/home/med-mai/Documents/tcgaBrca/manifest/brca/rnaseq"     
#[4] "/home/med-mai/Documents/tcgaBrca/manifest/brca/rnaseq_hg38"
#[5] "/home/med-mai/Documents/tcgaBrca/manifest/brca/snp6"       
#[6] "/home/med-mai/Documents/tcgaBrca/manifest/brca/wes"        

##manifest for TCGA ATAC-seq data files (Corces et al. 2018, Science)
list.files(paste0(MANIFEST,"/atac"),full.names=T)
#[1] "/home/med-mai/Documents/tcgaBrca/manifest/atac/ATACseq-AWG_Open_GDC-Manifest.txt"

##Tumor purity estimates for all TCGA tumors 
list.files(paste0(MANIFEST,"/pancan"),full.names=T)
#[1] "/home/med-mai/Documents/tcgaBrca/manifest/pancan/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
#[2] "/home/med-mai/Documents/tcgaBrca/manifest/pancan/TCGA-CDR-SupplementalTableS1.xlsx"          

##output directory for derived files 
dir.create(paste0(HOME,"/processedManifests"))

##add files for clinical, putity etc
dir.create(paste0(HOME,"/sampleAnno"))

file.copy(paste0(MANIFEST,"/pancan/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"),paste0(HOME,"/sampleAnno/"))

file.copy(paste0(MANIFEST,"/pancan/TCGA-CDR-SupplementalTableS1.xlsx"),paste0(HOME,"/sampleAnno/"))

file.copy(paste0(MANIFEST,"/pancan/clinical_PANCAN_patient_with_followup.tsv"),paste0(HOME,"/sampleAnno/"))

file.copy(paste0(MANIFEST,"/pancan/merged_sample_quality_annotations.tsv"),paste0(HOME,"/sampleAnno/"))

################################################################################
##load required packages

library(tools)

library(gridExtra)

################################################################################
##Define functions

##separate tab-sep lines 
doSep<-function(x) { unlist(strsplit(x,"\\t")) }

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

################################################################################
##Load data manifests and parse

##create big list of all sample-file combos with GDC-downloader manifest

##1. RNAseq
TEMP_PATH<-paste(MANIFEST,TUMOR_TYPE,"rnaseq_hg38",sep="/")

FILE<-list.files(TEMP_PATH,pattern="gdc_sample_sheet",full.names=T)
FILE2<-list.files(TEMP_PATH,pattern="manifest",full.names=T)

tmp<-readLines(FILE)
tmp2<-readLines(FILE2)

##F1
tmp.h<-tmp[1] ; tmp<-tmp[-1]
tmp<-do.call("rbind",lapply(tmp,doSep))
colnames(tmp)<- doUnSpace(doSep(tmp.h)) ; rm(tmp.h)
tmp<-as.data.frame(tmp,stringsAsFactors=F)

##F2
tmp2<-tmp2[-1]

#tmp2<-unlist(lapply(tmp2,doUnSep,x=tmp2.h))  ##v1

all.equal(unlist(lapply(tmp2,doUnPaste)),tmp$File_ID)
#[1] TRUE

##save manifest
tmp.gex<-tmp
tmp.gex$manifest<-gsub("\\t","|",tmp2) #change tab-sep to pipe for easier splitting

##2.1 SNP6 - ascat
TEMP_PATH<-paste(MANIFEST,TUMOR_TYPE,"snp6/ascat2",sep="/")

FILE<-list.files(TEMP_PATH,pattern="gdc_sample_sheet",full.names=T)
FILE2<-list.files(TEMP_PATH,pattern="manifest",full.names=T)

tmp<-readLines(FILE)
tmp2<-readLines(FILE2)

##F1
tmp.h<-tmp[1] ; tmp<-tmp[-1]
tmp<-do.call("rbind",lapply(tmp,doSep))
colnames(tmp)<- doUnSpace(doSep(tmp.h)) ; rm(tmp.h)
tmp<-as.data.frame(tmp,stringsAsFactors=F)

##F2
tmp2<-tmp2[-1]

#tmp2<-unlist(lapply(tmp2,doUnSep,x=tmp2.h))  ##v1

all.equal(unlist(lapply(tmp2,doUnPaste)),tmp$File_ID)
#[1] TRUE

##save manifest
tmp.ascat<-tmp
tmp.ascat$manifest<-gsub("\\t","|",tmp2)

##2.2 SNP6 - dnacopy
TEMP_PATH<-paste(MANIFEST,TUMOR_TYPE,"snp6/dnacopy",sep="/")

FILE<-list.files(TEMP_PATH,pattern="gdc_sample_sheet",full.names=T)
FILE2<-list.files(TEMP_PATH,pattern="manifest",full.names=T)

tmp<-readLines(FILE)
tmp2<-readLines(FILE2)

##F1
tmp.h<-tmp[1] ; tmp<-tmp[-1]
tmp<-do.call("rbind",lapply(tmp,doSep))
colnames(tmp)<- doUnSpace(doSep(tmp.h)) ; rm(tmp.h)
tmp<-as.data.frame(tmp,stringsAsFactors=F)

##F2
tmp2<-tmp2[-1]

#tmp2<-unlist(lapply(tmp2,doUnSep,x=tmp2.h))  ##v1

all.equal(unlist(lapply(tmp2,doUnPaste)),tmp$File_ID)
#[1] TRUE

##save manifest
tmp.dnacopy<-tmp
tmp.dnacopy$manifest<-gsub("\\t","|",tmp2)

##3.1 Image - dx
TEMP_PATH<-paste(MANIFEST,TUMOR_TYPE,"image/dx",sep="/")

FILE<-list.files(TEMP_PATH,pattern="gdc_sample_sheet",full.names=T)
FILE2<-list.files(TEMP_PATH,pattern="manifest",full.names=T)

tmp<-readLines(FILE)
tmp2<-readLines(FILE2)

##F1
tmp.h<-tmp[1] ; tmp<-tmp[-1]
tmp<-do.call("rbind",lapply(tmp,doSep))
colnames(tmp)<- doUnSpace(doSep(tmp.h)) ; rm(tmp.h)
tmp<-as.data.frame(tmp,stringsAsFactors=F)

##F2
tmp2<-tmp2[-1]

#tmp2<-unlist(lapply(tmp2,doUnSep,x=tmp2.h))  ##v1

all.equal(unlist(lapply(tmp2,doUnPaste)),tmp$File_ID)
#[1] TRUE

##save manifest
tmp.idx<-tmp
tmp.idx$manifest<-gsub("\\t","|",tmp2)

##3.2 Image - tissue
TEMP_PATH<-paste(MANIFEST,TUMOR_TYPE,"image/tissue",sep="/")

FILE<-list.files(TEMP_PATH,pattern="gdc_sample_sheet",full.names=T)
FILE2<-list.files(TEMP_PATH,pattern="manifest",full.names=T)

tmp<-readLines(FILE)
tmp2<-readLines(FILE2)

##F1
tmp.h<-tmp[1] ; tmp<-tmp[-1]
tmp<-do.call("rbind",lapply(tmp,doSep))
colnames(tmp)<- doUnSpace(doSep(tmp.h)) ; rm(tmp.h)
tmp<-as.data.frame(tmp,stringsAsFactors=F)

##F2
tmp2<-tmp2[-1]

#tmp2<-unlist(lapply(tmp2,doUnSep,x=tmp2.h))  ##v1

all.equal(unlist(lapply(tmp2,doUnPaste)),tmp$File_ID)
#[1] TRUE

##save manifest
tmp.its<-tmp
tmp.its$manifest<-gsub("\\t","|",tmp2)

##4 me450
TEMP_PATH<-paste(MANIFEST,TUMOR_TYPE,"me450",sep="/")

FILE<-list.files(TEMP_PATH,pattern="sdrf",full.names=T)
FILE2<-list.files(TEMP_PATH,pattern="manifest",full.names=T)

tmp<-readLines(FILE)
tmp2<-readLines(FILE2)

##F1 - sdrf with sample to file mapping
tmp.h<-tmp[1] ; tmp<-tmp[-1]
tmp<-do.call("rbind",lapply(tmp,doSep))
colnames(tmp)<- doUnSpace(doSep(tmp.h)) ; rm(tmp.h)
tmp<-as.data.frame(tmp,stringsAsFactors=F)

##F2
tmp2<-tmp2[-1]
#tmp2<-unlist(lapply(tmp2,doUnSep,x=tmp2.h))

i1<-intersect( getIdat(tmp2),tmp$Array_Data_File )

tmp<-tmp[match(i1,tmp$Array_Data_File),]
tmp2<-tmp2[match(i1,getIdat(tmp2))]
rm(i1)

##file is legacy and does not have same structure as uniform GDC-files
	##build format similar to other data types
tmpN<-data.frame(
	matrix("",
			nrow=length(tmp2),
			ncol=length(doSep(addSampleHeader(NULL)))
		),stringsAsFactors =FALSE
)
colnames(tmpN)<-doSep(addSampleHeader(NULL))

##generic columns
tmpN$Data_Category<-"Methylation Profiling"
tmpN$Data_Type<-"Illumina 450K methylation"
tmpN$Project_ID<-idToTcga(TUMOR_TYPE)
tmpN$Data_Category<-"DNA methylation"

##specific columns
tmpN$File_ID<-doUnPaste(tmp2) ##unlike other data this only has fileID in manifest
tmpN$File_Name<-getIdat(tmp2)
tmpN$Case_ID<-parseTcgaName(tmp$Comment_TCGA_Barcode,level="patient")
tmpN$Sample_ID<-parseTcgaName(tmp$Comment_TCGA_Barcode,level=4)
tmpN$Sample_Type<-sampleDecode(parseTcgaName(tmp$Comment_TCGA_Barcode,level="sample_type"))

all.equal(getIdat(tmp2),tmpN$File_Name)
#[1] TRUE
all.equal(getIdat(tmp2),tmp$Array_Data_File)
#[1] TRUE

##save manifest
tmp.me<-tmpN
tmp.me$manifest<-gsub("\\t","|",tmp2)

rm(FILE,FILE2,TEMP_PATH)
rm(tmp,tmp2,tmpN)

##5 WES -> Temp-get file to grab sample ID's

#get manifest
TEMP_PATH<-paste(MANIFEST,TUMOR_TYPE,"wes",sep="/")
FILE<-list.files(TEMP_PATH,pattern="gdc_manifest_wes",full.names=T)
FILE2<-readLines(FILE)

#get tempdir for file download
temp<-tempdir()

if( Sys.info()["sysname"] == "Windows" ) {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client.exe",MANIFEST),temp)
	file.copy(FILE,temp)
	setwd(temp)

	} else {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
	file.copy(FILE,temp)
	setwd(temp)

	#do syscommand to do download
	system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )

	}

cmd<-paste0(
  paste0(getwd(),"/",list.files(".",pattern="gdc-client")),
  paste0(" download -m  ",
  paste0(getwd(),"/",list.files(".",pattern="manifest"))
  ))
system(cmd)

#clean
unlink(paste0(getwd(),"/",list.files(".",pattern="gdc-client")))
unlink(paste0(getwd(),"/",list.files(".",pattern="manifest")))
rm(cmd)

#get file links
tmp<-paste0(temp,"/",unlist(lapply(list.files(),function(x) list.files(x,pattern="maf.gz",full.names=TRUE))))

#back and clean
setwd("~/")
rm(temp)

tmp.files<-lapply(tmp,readLines)

# tmp.samples<-lapply(tmp.files,function(x) {
# 	setdiff(sub(".+(TCGA-..-....).+|.+","\\1",x),c(""))
# 	})

tmp.samples<-lapply(tmp.files,function(x) {
	setdiff(sub(".+(TCGA-..-....-...-...-....-..)\\t(TCGA-..-....-...-...-....-..).+|.+","\\1|\\2",x),c("","|"))
	})

tmp.samples2<-lapply(1:length(tmp.samples),function(x) {
	y<-unique(tmp.samples[[x]])
	paste(y,
		sub(".+TCGA\\.....\\.(\\w+)\\..+","\\1",tmp[x]),
		sub(".+(TCGA\\..+)","\\1",tmp[x]),
		sep="|")
})

tmp.samples2<-unlist(tmp.samples2)
tmp.samples2<-tmp.samples2[order(tmp.samples2)]

tmpN<-data.frame(
	matrix("",
			nrow=length(tmp.samples2),
			ncol=length(doSep(addSampleHeader(NULL)))
		),stringsAsFactors =FALSE
)
colnames(tmpN)<-doSep(addSampleHeader(NULL))

##generic columns
tmpN$Data_Category<-"Whole exome sequencing"
tmpN$Data_Type<-"Somatic mutation"
tmpN$Project_ID<-idToTcga(TUMOR_TYPE)

##specific columns
tmpN$File_Name<-sub(".+\\|","",tmp.samples2)
tmpN$File_ID<-sub("\\t.+","",FILE2[unlist(lapply(tmpN$File_Name,function(x) grep(x,FILE2)))])

tmpN$Data_Type<-paste(tmpN$Data_Type,sub("TCGA\\.....\\.(\\w+)\\..+","\\1",tmpN$File_Name),sep=" - ")

tmpN$Case_ID<-parseTcgaName(sub("\\|.+","",tmp.samples2),level="patient")
tmpN$Sample_ID<-parseTcgaName(sub("\\|.+","",tmp.samples2),level=4)
tmpN$Sample_Type<-sampleDecode(parseTcgaName(sub("\\|.+","",tmp.samples2),level="sample_type"))

##save manifest
tmp.wes<-tmpN
tmp.wes$manifest<-gsub("\\t","|",FILE2[unlist(lapply(tmpN$File_Name,function(x) grep(x,FILE2)))])

rm(FILE,FILE2,TEMP_PATH)
rm(tmp,tmp.files,tmp.samples,tmp.samples2,tmpN)

##6 ATAC-seq -> Temp-get id-mappings

#get manifest
TEMP_PATH<-paste(MANIFEST,"atac",sep="/")
FILE<-list.files(TEMP_PATH,pattern="Manifest",full.names=T)
#get id-mapping
FILE2<-readLines(FILE)
FILE2<-c(FILE2[1],FILE2[grep("TCGA_identifier_mapping.txt",FILE2)])
#get normalized counts zip-file - only manifest for later download..
FILE3<-readLines(FILE)
FILE3<-FILE3[grep("TCGA-ATAC_Cancer_Type-specific_Count_Matrices_log2norm_counts.zip",FILE3)]

#get tempdir for file download
temp<-tempdir()

if( Sys.info()["sysname"] == "Windows" ) {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client.exe",MANIFEST),temp)
	writeLines(FILE2,paste0(temp,"/manifest_atac.txt"))
	setwd(temp)

	} else {

	#transfer to tempdir
	file.copy(sub("manifest","gdcTransferTool/gdc-client",MANIFEST),temp)
	writeLines(FILE2,paste0(temp,"/manifest_atac.txt"))
	setwd(temp)

	#do syscommand to do download
	system( paste0('chmod 777 ',paste0(getwd(),"/",list.files(".",pattern="gdc-client"))) )

	}

cmd<-paste0(
  paste0(getwd(),"/",list.files(".",pattern="gdc-client")),
  paste0(" download -m  ",
  paste0(getwd(),"/",list.files(".",pattern="manifest_atac"))
  ))
system(cmd)

#clean
unlink(paste0(getwd(),"/",list.files(".",pattern="gdc-client")))
unlink(paste0(getwd(),"/",list.files(".",pattern="manifest_atac")))
rm(cmd)

#get file links
tmp<-paste0(temp,"/",unlist(lapply(list.files(),function(x) list.files(x,pattern="TCGA_identifier_mapping.txt",full.names=TRUE))))

tmp2<-readLines(tmp)

table(sub("(\\w\\w\\w\\w)-.+|.+","\\1",tmp2))
#      ACCx BLCA BRCA CESC CHOL COAD ESCA GBMx HNSC KIRC KIRP LGGx LIHC LUAD LUSC 
#    1   18   20  141    7   10   81   33   18   16   32   68   24   34   44   32 
# MESO PCPG PRAD SKCM STAD TGCT THCA UCEC 
#   13   18   52   24   41   18   28   24 

tmp2<-tmp2[grep(TUMOR_TYPE,tmp2,ignore.case=TRUE)]

tmpN<-data.frame(
	matrix("",
			nrow=length(tmp2),
			ncol=length(doSep(addSampleHeader(NULL)))
		),stringsAsFactors =FALSE
)
colnames(tmpN)<-doSep(addSampleHeader(NULL))

##generic columns
tmpN$Data_Category<-"ATAC sequencing"
tmpN$Data_Type<-"log2 normalized peak counts"
tmpN$Project_ID<-idToTcga(TUMOR_TYPE)

##specific columns
#atac data not separated by sample - one matrix for all
tmpN$File_Name<-sub("(.+)\\t(.+)\\t(.+)\\t(.+)","\\2",FILE3)
tmpN$File_ID<-sub("(.+)\\t(.+)\\t(.+)\\t(.+)","\\1",FILE3)

#has Stanfor and TCGA id's
tmpN$Case_ID<-parseTcgaName(sub(".+\\t","",tmp2),level="atac_patient")
tmpN$Sample_ID<-parseTcgaName(sub(".+\\t","",tmp2),level="atac_sample")
tmpN$Sample_Type<-sampleDecode(parseTcgaName(sub(".+\\t","",tmp2),level="atac_type"))

##save manifest
tmp.atac<-tmpN
tmp.atac$manifest<-gsub("\\t","|",FILE3)

##save sample mapping Stanford-TCGA - put in processedManifests
writeLines( paste(sub("\\t.+","",tmp2),parseTcgaName(sub(".+\\t","",tmp2),level="atac_id"),sep="\t",collapse="\n"),
	con=paste0(HOME,"/","processedManifests/","stanfordToTcga_sampleMapping.txt")
)

rm(FILE,FILE2,FILE3,TEMP_PATH)
rm(tmp,tmp2,tmpN,temp)

################################################################################
##Gather tumor and manifest data and get sample statistics

obj<-ls(pattern="tmp\\.")

obj
# [1] "tmp.ascat"   "tmp.atac"    "tmp.dnacopy" "tmp.gex"     "tmp.idx"    
# [6] "tmp.its"     "tmp.me"      "tmp.wes"   

##Get unique sample ID's and map data levels to these
samples.all<-unlist(lapply(obj,function(x) {
	y<-get(x)
	return( unique(sub("(TCGA-..-....).+","\\1",y$Case_ID)) )
}))

length(samples.all)
#[1] 7263
length(unique(samples.all))
#[1] 1098

samples.all<-unique(samples.all)

samples.overlap<-do.call("cbind",lapply(obj,function(x) {
	y<-get(x)
	samples.all %in% unique(sub("(TCGA-..-....).+","\\1",y$Case_ID))
}))

colnames(samples.overlap)<-sub("tmp\\.","",obj)
rownames(samples.overlap)<-samples.all

crosstab<-matrix(0,
			nrow=ncol(samples.overlap),
			ncol=ncol(samples.overlap)
		)
colnames(crosstab)<-colnames(samples.overlap)
rownames(crosstab)<-colnames(samples.overlap)

crosstab[lower.tri(crosstab)]<-apply(combn(ncol(samples.overlap),2),2,function(x) sum(samples.overlap[,x[1]] & samples.overlap[,x[2]]) )
diag(crosstab)<-colSums(samples.overlap)
crosstab<-as.data.frame(crosstab,stringsAsFactors=FALSE)

#output image with crosstab stats
png(paste0(HOME,"/","processedManifests/","data_summary.png"),width = 300*6, height = 300*3, units = "px",res=300)
par(mar=c(0.1,0.1,0.1,0.1))
grid.table(crosstab)
dev.off()

rm(crosstab,samples.all,obj)

################################################################################
##Gather final manifests

##write out matrices 
ls()
#  [1] "addDownloadHeader" "addSampleHeader"   "doSep"            
#  [4] "doUnPaste"         "doUnSpace"         "getIdat"          
#  [7] "HOME"              "idToTcga"          "MANIFEST"         
# [10] "parseTcgaName"     "sampleDecode"      "samples.overlap"  
# [13] "tmp.ascat"         "tmp.atac"          "tmp.dnacopy"      
# [16] "tmp.gex"           "tmp.idx"           "tmp.its"          
# [19] "tmp.me"            "tmp.wes"           "TUMOR_TYPE"       

ls(pattern="tmp\\.")
# [1] "tmp.ascat"   "tmp.atac"    "tmp.dnacopy" "tmp.gex"     "tmp.idx"    
# [6] "tmp.its"     "tmp.me"      "tmp.wes"    

##core data manifest
	##gex + me450 + ascat copy + wes

sel <- rowSums( samples.overlap[,c("ascat","gex","me","wes")] ) >= 1

table( sel )
# sel
# TRUE 
# 1098 

samples.core<-rownames(samples.overlap[sel,])

##write to file

#ascat
write.table(tmp.ascat[sub("(TCGA-..-....).+","\\1",tmp.ascat$Case_ID) %in% samples.core,],
	file=paste0(HOME,"/","processedManifests/","coreSetManifest_ascat.txt"),
		sep="\t",quote=FALSE,row.names=FALSE
		)

#gex
write.table(tmp.gex[tmp.gex$Case_ID %in% samples.core,],
	file=paste0(HOME,"/","processedManifests/","coreSetManifest_gex.txt"),
		sep="\t",quote=FALSE,row.names=FALSE
		)

#me
write.table(tmp.me[tmp.me$Case_ID %in% samples.core,],
	file=paste0(HOME,"/","processedManifests/","coreSetManifest_me.txt"),
		sep="\t",quote=FALSE,row.names=FALSE
		)

#wes
write.table(tmp.wes[tmp.wes$Case_ID %in% samples.core,],
	file=paste0(HOME,"/","processedManifests/","coreSetManifest_wes.txt"),
		sep="\t",quote=FALSE,row.names=FALSE
		)

##extended data manifest
	##gex + me450 + ascat copy + wes + atac
sel <- rowSums( samples.overlap[,c("ascat","atac","gex","me","wes")] ) == 5

table( sel )
# sel
# FALSE  TRUE 
#  1059    39 

samples.ext<-rownames(samples.overlap[sel,])

rm(sel)

#atac
write.table(tmp.atac[sub("(TCGA-..-....).+","\\1",tmp.atac$Case_ID) %in% samples.ext,],
	file=paste0(HOME,"/","processedManifests/","extendedSetManifest_atac.txt"),
		sep="\t",quote=FALSE,row.names=FALSE
		)

save.image(paste0(HOME,"/","processedManifests/","manifestObject.RData") )

q("no")
###END