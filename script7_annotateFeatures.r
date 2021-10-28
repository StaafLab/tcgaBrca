#####======================================================================#####
### Create feature annotation object(s) for 450K platform
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
if( !file.exists( paste0(HOME,"/","annotateFeatures") )   ) {
	dir.create(paste0(HOME,"/","annotateFeatures"),recursive=TRUE)
}

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

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db") }

library("org.Hs.eg.db")

writeLines(apply(metadata(org.Hs.eg.db),1,paste,collapse=": "),con=paste0(HOME,"/annotateFeatures/org.Hs.eg.db_versionInfo.txt"))

if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene") }

library("TxDb.Hsapiens.UCSC.hg38.knownGene")

writeLines(apply(metadata(TxDb.Hsapiens.UCSC.hg38.knownGene),1,paste,collapse=": "),con=paste0(HOME,"/annotateFeatures/TxDb.Hsapiens.UCSC.hg38.knownGene_versionInfo.txt"))

##md5
if(!requireNamespace("tools", quietly = TRUE)) {
  install.packages("tools") }

library(tools)

if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38") }

library("BSgenome.Hsapiens.UCSC.hg38")

sink(file=paste0(HOME,"/annotateFeatures/BSgenome.Hsapiens.UCSC.hg38_versionInfo.txt"))
BSgenome.Hsapiens.UCSC.hg38
sink()

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings") }

library("Biostrings")

################################################################################
##Get core tumor set from methylation data

load(file=paste0(HOME,"/","workspace_blacklistFiltered_atacCnGexMeWes_withSampleAnnotations.RData"))

ls()
#  [1] "betaOrig"   "dataAtac"   "dataCn"     "dataMut"    "dataSeg"   
#  [6] "gexCounts"  "gexFpkm"    "gexUq"      "HOME"       "MANIFEST"  
# [11] "sampleAnno" "TUMOR_TYPE"

#object - "betaAdj" - purity adjusted betas
load(file=paste0(HOME,"/me/data450k_",nrow(betaOrig),"x",ncol(betaOrig),"_minfiNormalized_ringnerAdjusted_purityAdjusted_tumorBetaValues.RData"))

str(betaAdj)
 # num [1:421368, 1:630] 0.742 0.756 0.118 0.967 0.904 0.7 0.966 1 0.985 1 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:421368] "cg21870274" "cg08258224" "cg16619049" "cg18147296" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

#object - "betaNorm" - inferred normal methylation
load(file=paste0(HOME,"/me/data450k_",nrow(betaOrig),"x",ncol(betaOrig),"_minfiNormalized_ringnerAdjusted_purityAdjusted_normalBetaValues.RData"))

str(betaNorm)
 # num [1:421368, 1:630] 0.584 0.824 0.59 0.686 0.662 0.313 0.861 0.997 0.992 0.996 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:421368] "cg21870274" "cg08258224" "cg16619049" "cg18147296" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...
 
##get probesKeep
load(file=paste0(HOME,"/me/object_450k_probesKeep.RData"))

################################################################################
###Get transcript/gene data for TCGA GEX 

str(gexFpkm)
 # num [1:60483, 1:630] 0 0 3.68 0 1.83 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:60483] "ENSG00000242268.2" "ENSG00000270112.3" "ENSG00000167578.15" "ENSG00000273842.1" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

identical(rownames(gexFpkm),rownames(gexCounts))
#[1] TRUE

identical(rownames(gexFpkm),rownames(gexUq))
#[1] TRUE

##https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files

##Use GDC files - GDC.h38 GENCODE v22 GTF (used in RNA-Seq alignment and by HTSeq)
# tmp<-tempfile()
# download.file("https://api.gdc.cancer.gov/data/25aa497c-e615-4cb7-8751-71f744f9691f",tmp)

#as.vector(md5sum(tmp))

# genes.gtf<-readLines(tmp)

##Use GDC files - GDC.h38 GENCODE TSV
tmp<-tempfile()
download.file("https://api.gdc.cancer.gov/data/b011ee3e-14d8-4a97-aed4-e0b10f6bbe82",tmp)

as.vector(md5sum(tmp))
#[1] "3c1b05ce5c978370ed5e1c2f074f44b0"

genes.csv<-read.table(tmp,header=T,sep="\t",as.is=TRUE)

str(genes.csv)
# 'data.frame':   60483 obs. of  12 variables:
#  $ gene_id    : chr  "ENSG00000223972.5" "ENSG00000238009.5" "ENSG00000230415.1" "ENSG00000236335.1" ...
#  $ gene_name  : chr  "DDX11L1" "RP11-34P13.7" "RP5-902P8.10" "RP4-591L5.1" ...
#  $ seqname    : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start      : int  11869 89295 1275223 30409560 32752910 32785646 32818018 32951574 32996608 33033863 ...
#  $ end        : int  14409 133723 1280420 30411638 32753901 32786116 32897826 32956349 33097230 33034810 ...
#  $ strand     : chr  "+" "-" "+" "-" ...
#  $ gene_type  : chr  "transcribed_unprocessed_pseudogene" "lincRNA" "lincRNA" "lincRNA" ...
#  $ gene_status: chr  "KNOWN" "NOVEL" "NOVEL" "NOVEL" ...
#  $ havana_gene: chr  "OTTHUMG00000000961.2" "OTTHUMG00000001096.2" "OTTHUMG00000002234.2" "OTTHUMG00000003682.1" ...
#  $ full_length: int  2541 44429 5198 2079 992 471 79809 4776 100623 948 ...
#  $ exon_length: int  1735 3726 513 507 992 471 8685 3095 4364 948 ...
#  $ exon_num   : int  9 17 5 3 1 1 4 2 123 1 ...

length(intersect(rownames(gexFpkm),genes.csv$gene_id))
#[1] 60483

all(genes.csv$gene_id %in% rownames(gexFpkm) )
#[1] TRUE
all(rownames(gexFpkm) %in% genes.csv$gene_id )
#[1] TRUE

length(unique(genes.csv$gene_id))
#[1] 60483

length(unique(sub("\\..*","",genes.csv$gene_id)))
#[1] 60483
head(unique(sub("\\..*","",genes.csv$gene_id)))
#[1] "ENSG00000223972" "ENSG00000238009" "ENSG00000230415" "ENSG00000236335" "ENSG00000213842" "ENSG00000227337"

##every ensg is unique
rownames(genes.csv)<-sub("\\..*","",genes.csv$gene_id)
genes.csv$gene_id2<-sub("\\..*","",genes.csv$gene_id)

length(unique(genes.csv$gene_name))
#[1] 58387

table(genes.csv$gene_type)
 #          3prime_overlapping_ncrna                          antisense                          IG_C_gene                    IG_C_pseudogene 
 #                                29                               5565                                 14                                  9 
 #                         IG_D_gene                          IG_J_gene                    IG_J_pseudogene                          IG_V_gene 
 #                                37                                 18                                  3                                142 
 #                   IG_V_pseudogene                            lincRNA                       macro_lncRNA                              miRNA 
 #                               180                               7656                                  1                               4093 
 #                          misc_RNA                            Mt_rRNA                            Mt_tRNA                         non_coding 
 #                              2298                                  2                                 22                                  3 
 #            polymorphic_pseudogene               processed_pseudogene               processed_transcript                     protein_coding 
 #                                59                              10304                                484                              19814 
 #                        pseudogene                           ribozyme                               rRNA                             scaRNA 
 #                                36                                  8                                544                                 49 
 #                    sense_intronic                  sense_overlapping                             snoRNA                              snRNA 
 #                               920                                197                                961                               1896 
 #                              sRNA                                TEC                          TR_C_gene                          TR_D_gene 
 #                                20                               1045                                  5                                  3 
 #                         TR_J_gene                    TR_J_pseudogene                          TR_V_gene                    TR_V_pseudogene 
 #                                73                                  4                                106                                 30 
 #  transcribed_processed_pseudogene     transcribed_unitary_pseudogene transcribed_unprocessed_pseudogene    translated_processed_pseudogene 
 #                               443                                  1                                663                                  1 
 # translated_unprocessed_pseudogene                 unitary_pseudogene             unprocessed_pseudogene                           vaultRNA 
 #                                 1                                169                               2574                                  1 

table(genes.csv$gene_status)
   # KNOWN    NOVEL PUTATIVE 
   # 41346    19008      129 

table(genes.csv$gene_type,genes.csv$gene_status)
  #                                    KNOWN NOVEL PUTATIVE
  # 3prime_overlapping_ncrna               1    28        0
  # antisense                            244  5321        0
  # IG_C_gene                             14     0        0
  # IG_C_pseudogene                        9     0        0
  # IG_D_gene                             37     0        0
  # IG_J_gene                             18     0        0
  # IG_J_pseudogene                        3     0        0
  # IG_V_gene                            137     5        0
  # IG_V_pseudogene                      180     0        0
  # lincRNA                              413  7231       12
  # macro_lncRNA                           1     0        0
  # miRNA                               1707  2386        0
  # misc_RNA                             948  1350        0
  # Mt_rRNA                                2     0        0
  # Mt_tRNA                               22     0        0
  # non_coding                             3     0        0
  # polymorphic_pseudogene                59     0        0
  # processed_pseudogene               10302     2        0
  # processed_transcript                  66   415        3
  # protein_coding                     19308   392      114
  # pseudogene                             6    30        0
  # ribozyme                               3     5        0
  # rRNA                                 516    28        0
  # scaRNA                                18    31        0
  # sense_intronic                        25   895        0
  # sense_overlapping                      7   190        0
  # snoRNA                               382   579        0
  # snRNA                               1798    98        0
  # sRNA                                   0    20        0
  # TEC                                 1045     0        0
  # TR_C_gene                              5     0        0
  # TR_D_gene                              3     0        0
  # TR_J_gene                             72     1        0
  # TR_J_pseudogene                        4     0        0
  # TR_V_gene                            105     1        0
  # TR_V_pseudogene                       30     0        0
  # transcribed_processed_pseudogene     443     0        0
  # transcribed_unitary_pseudogene         1     0        0
  # transcribed_unprocessed_pseudogene   663     0        0
  # translated_processed_pseudogene        1     0        0
  # translated_unprocessed_pseudogene      1     0        0
  # unitary_pseudogene                   169     0        0
  # unprocessed_pseudogene              2574     0        0
  # vaultRNA                               1     0        0

##use all - can select subsets for analysis if needed
  ##add ENTREZ as well for those with available

##With entrez and using fpkm-rownames
db<-org.Hs.eg.db

id<-sub("\\..*","",rownames(gexFpkm))

head(id)
# [1] "ENSG00000242268" "ENSG00000270112" "ENSG00000167578" "ENSG00000273842"
# [5] "ENSG00000078237" "ENSG00000146083"

keytypes(org.Hs.eg.db)
#  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
#  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"        
# [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
# [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"     
# [26] "UNIPROT"    

#sym<-mapIds(db, keys=id, column=c("SYMBOL"), keytype="ENSEMBL",multiVals="first")
sym<-mapIds(db, keys=id, column=c("ENTREZID"), keytype="ENSEMBL",multiVals="first")

str(sym)
 # Named chr [1:60483] "100507661" NA "53916" NA "57103" "22838" NA "55567" ...
 # - attr(*, "names")= chr [1:60483] "ENSG00000242268" "ENSG00000270112" "ENSG00000167578" "ENSG00000273842" ...

length(sym)
#[1] 60483

length(unique(sym))
#[1] 25555

table(is.na(sym))
# FALSE  TRUE 
# 25629 34854 

length(unique(names(sym)))
#[1] 60483

sym[is.na(sym)]<-""

table(sym=="")
# FALSE  TRUE 
# 25629 34854 

all(names(sym) %in% genes.csv$gene_id2)
#[1] TRUE

all(genes.csv$gene_id2 %in% names(sym))
#[1] TRUE

length(intersect(genes.csv$gene_id2 , names(sym))) == length(sym)
#[1] TRUE

sym<-sym[genes.csv$gene_id2]

all(names(sym) == genes.csv$gene_id2)
#[1] TRUE

genes.csv$orgHs_entrez<-as.vector(sym)

##Add orgHs - symbol
sym<-mapIds(db, keys=id, column=c("SYMBOL"), keytype="ENSEMBL",multiVals="first")

str(sym)
 # Named chr [1:60483] "LINC02082" NA "RAB4B" NA "TIGAR" "RNF44" NA "DNAH3" "RPL23A" "EHD4-AS1" NA NA NA "ARL8B" NA "CALB2" "MFSD3" NA "LINC00636" "PIGV" ...
 # - attr(*, "names")= chr [1:60483] "ENSG00000242268" "ENSG00000270112" "ENSG00000167578" "ENSG00000273842" ...

length(sym)
#[1] 60483

length(unique(sym))
#[1] 25555

table(is.na(sym))
# FALSE  TRUE 
# 25629 34854 

length(unique(names(sym)))
#[1] 60483

sym[is.na(sym)]<-""

table(sym=="")
# FALSE  TRUE 
# 25629 34854 

all(names(sym) %in% genes.csv$gene_id2)
#[1] TRUE

all(genes.csv$gene_id2 %in% names(sym))
#[1] TRUE

length(intersect(genes.csv$gene_id2 , names(sym))) == length(sym)
#[1] TRUE

sym<-sym[genes.csv$gene_id2]

all(names(sym) == genes.csv$gene_id2)
#[1] TRUE

genes.csv$orgHs_symbol<-as.vector(sym)

table(genes.csv$gene_name==genes.csv$orgHs_symbol)
# FALSE  TRUE 
# 38377 22106 

table(genes.csv$orgHs_entrez!="",genes.csv$orgHs_symbol!="")
  #       FALSE  TRUE
  # FALSE 34854     0
  # TRUE      0 25629

table(genes.csv$orgHs_entrez!="",genes.csv$gene_status)
  #       KNOWN NOVEL PUTATIVE
  # FALSE 19313 15440      101
  # TRUE  22033  3568       28

table(genes.csv$gene_type,genes.csv$orgHs_entrez!="")
  #                                    FALSE  TRUE
  # 3prime_overlapping_ncrna              28     1
  # antisense                           4137  1428
  # IG_C_gene                             14     0
  # IG_C_pseudogene                        9     0
  # IG_D_gene                             37     0
  # IG_J_gene                             18     0
  # IG_J_pseudogene                        3     0
  # IG_V_gene                            142     0
  # IG_V_pseudogene                      180     0
  # lincRNA                             5500  2156
  # macro_lncRNA                           1     0
  # miRNA                               2675  1418
  # misc_RNA                            2286    12
  # Mt_rRNA                                2     0
  # Mt_tRNA                               21     1
  # non_coding                             3     0
  # polymorphic_pseudogene                14    45
  # processed_pseudogene               10203   101
  # processed_transcript                 257   227
  # protein_coding                       706 19108
  # pseudogene                            36     0
  # ribozyme                               7     1
  # rRNA                                 524    20
  # scaRNA                                31    18
  # sense_intronic                       842    78
  # sense_overlapping                    172    25
  # snoRNA                               540   421
  # snRNA                               1853    43
  # sRNA                                  19     1
  # TEC                                 1022    23
  # TR_C_gene                              4     1
  # TR_D_gene                              3     0
  # TR_J_gene                             73     0
  # TR_J_pseudogene                        4     0
  # TR_V_gene                            106     0
  # TR_V_pseudogene                       30     0
  # transcribed_processed_pseudogene     350    93
  # transcribed_unitary_pseudogene         0     1
  # transcribed_unprocessed_pseudogene   362   301
  # translated_processed_pseudogene        1     0
  # translated_unprocessed_pseudogene      0     1
  # unitary_pseudogene                   118    51
  # unprocessed_pseudogene              2520    54
  # vaultRNA                               1     0

str(genes.csv)
# 'data.frame':   60483 obs. of  15 variables:
#  $ gene_id     : chr  "ENSG00000223972.5" "ENSG00000238009.5" "ENSG00000230415.1" "ENSG00000236335.1" ...
#  $ gene_name   : chr  "DDX11L1" "RP11-34P13.7" "RP5-902P8.10" "RP4-591L5.1" ...
#  $ seqname     : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start       : int  11869 89295 1275223 30409560 32752910 32785646 32818018 32951574 32996608 33033863 ...
#  $ end         : int  14409 133723 1280420 30411638 32753901 32786116 32897826 32956349 33097230 33034810 ...
#  $ strand      : chr  "+" "-" "+" "-" ...
#  $ gene_type   : chr  "transcribed_unprocessed_pseudogene" "lincRNA" "lincRNA" "lincRNA" ...
#  $ gene_status : chr  "KNOWN" "NOVEL" "NOVEL" "NOVEL" ...
#  $ havana_gene : chr  "OTTHUMG00000000961.2" "OTTHUMG00000001096.2" "OTTHUMG00000002234.2" "OTTHUMG00000003682.1" ...
#  $ full_length : int  2541 44429 5198 2079 992 471 79809 4776 100623 948 ...
#  $ exon_length : int  1735 3726 513 507 992 471 8685 3095 4364 948 ...
#  $ exon_num    : int  9 17 5 3 1 1 4 2 123 1 ...
#  $ gene_id2    : chr  "ENSG00000223972" "ENSG00000238009" "ENSG00000230415" "ENSG00000236335" ...
#  $ orgHs_entrez: chr  "100287102" "100996442" "" "" ...
#  $ orgHs_symbol: chr  "DDX11L1" "LOC100996442" "" "" ...

##check chromosomes
table(genes.csv$seqname)
 # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrM  chrX  chrY 
 # 5397  2306  3381  3047  1383  2289  2247  2597  3111  1206  2997  4150  1436   880  1385  3163  2633  2993  3001  2980  2444  2350    37  2476   594 

##are coordinates BED?
all(genes.csv$end-genes.csv$start == genes.csv$full_length)
#[1] FALSE

##are coordinates actual?
all(1+(genes.csv$end-genes.csv$start) == genes.csv$full_length)
#[1] TRUE

##Make genes GRanges
geneCoords<-GRanges(seqnames=genes.csv$seqname,
  ranges=IRanges(start=genes.csv$start,
    end=genes.csv$end),
  strand=genes.csv$strand
)
names(geneCoords)<-genes.csv$gene_id

geneCoords$gene_id<-genes.csv$gene_id
geneCoords$gene_id2<-genes.csv$gene_id2
geneCoords$gene_name<-genes.csv$gene_name
geneCoords$orgHs_entrez<-genes.csv$orgHs_entrez
geneCoords$orgHs_symbol<-genes.csv$orgHs_symbol
geneCoords$gene_type<-genes.csv$gene_type
geneCoords$gene_status<-genes.csv$gene_status

geneCoords
# GRanges object with 60483 ranges and 7 metadata columns:
#                      seqnames            ranges strand |            gene_id        gene_id2    gene_name orgHs_entrez orgHs_symbol
#                         <Rle>         <IRanges>  <Rle> |        <character>     <character>  <character>  <character>  <character>
#    ENSG00000223972.5     chr1       11869-14409      + |  ENSG00000223972.5 ENSG00000223972      DDX11L1    100287102      DDX11L1
#    ENSG00000238009.5     chr1      89295-133723      - |  ENSG00000238009.5 ENSG00000238009 RP11-34P13.7    100996442 LOC100996442
#    ENSG00000230415.1     chr1   1275223-1280420      + |  ENSG00000230415.1 ENSG00000230415 RP5-902P8.10                          
#    ENSG00000236335.1     chr1 30409560-30411638      - |  ENSG00000236335.1 ENSG00000236335  RP4-591L5.1                          
#    ENSG00000213842.2     chr3 32752910-32753901      + |  ENSG00000213842.2 ENSG00000213842      SUGT1P2                          
#                  ...      ...               ...    ... .                ...             ...          ...          ...          ...
#    ENSG00000226955.1     chr3 32592972-32593312      - |  ENSG00000226955.1 ENSG00000226955   AC104306.4                          
#    ENSG00000231552.1     chr3 32620903-32621902      - |  ENSG00000231552.1 ENSG00000231552      IGBP1P3                          
#    ENSG00000237676.1     chr3 32635222-32635554      - |  ENSG00000237676.1 ENSG00000237676      RPL30P4                          
#   ENSG00000182973.17     chr3 32685145-32773875      + | ENSG00000182973.17 ENSG00000182973       CNOT10        25904       CNOT10
#    ENSG00000251224.1     chr3 32730635-32737454      - |  ENSG00000251224.1 ENSG00000251224   CNOT10-AS1                          
#                                               gene_type gene_status
#                                             <character> <character>
#    ENSG00000223972.5 transcribed_unprocessed_pseudogene       KNOWN
#    ENSG00000238009.5                            lincRNA       NOVEL
#    ENSG00000230415.1                            lincRNA       NOVEL
#    ENSG00000236335.1                            lincRNA       NOVEL
#    ENSG00000213842.2               processed_pseudogene       KNOWN
#                  ...                                ...         ...
#    ENSG00000226955.1               processed_pseudogene       KNOWN
#    ENSG00000231552.1               processed_pseudogene       KNOWN
#    ENSG00000237676.1               processed_pseudogene       KNOWN
#   ENSG00000182973.17                     protein_coding       KNOWN
#    ENSG00000251224.1                          antisense       NOVEL
#   -------
#   seqinfo: 25 sequences from an unspecified genome; no seqlengths

save(geneCoords,genes.csv,file=paste0(HOME,"/annotateFeatures/object_gdc_GENCODE22_gene_annotations.RData"))

################################################################################
###Make GRanges for TCGA ATAC

str(dataAtac)
 # num [1:215920, 1:39] 1.031 1.667 1.839 1.375 0.769 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:215920] "chr1:17234-17733|BRCA_2" "chr1:180634-181133|BRCA_3" "chr1:181207-181706|BRCA_4" "chr1:183557-184056|BRCA_5" ...
 #  ..$ : chr [1:39] "TCGA-3C-AALJ-01A" "TCGA-4H-AAAK-01A" "TCGA-A2-A0SV-01A" "TCGA-A2-A0SW-01A" ...

table(sub(":.+","",rownames(dataAtac)))
#  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX 
# 21640  9733 11414 10140  4142  6186  6812  7042 13211  4038  9379 15782  7036  2942  4259 12401  7463  9337 15170 10545 15026  8240  3982 

quantile(as.integer(sub("chr.+:(.+)-.+","\\1",rownames(dataAtac))))
    #    0%       25%       50%       75%      100% 
    # 10123  33926919  65673646 114911183 248945392 
quantile(as.integer(sub("chr.+:(.+)-(.+)\\|.+","\\2",rownames(dataAtac))))
    #    0%       25%       50%       75%      100% 
    # 10622  33927418  65674145 114911682 248945891 

##500 bp peaks
quantile(as.integer(sub("chr.+:(.+)-(.+)\\|.+","\\2",rownames(dataAtac))) - 
  as.integer(sub("chr.+:(.+)-.+","\\1",rownames(dataAtac)))
  )
 #  0%  25%  50%  75% 100% 
 # 499  499  499  499  499 

atacCoords<-data.frame(seqnames=sub(":.+","",rownames(dataAtac)),
  start=as.integer(sub("chr.+:(.+)-.+","\\1",rownames(dataAtac))),
  end=as.integer(sub("chr.+:(.+)-(.+)\\|.+","\\2",rownames(dataAtac))),
  stringsAsFactors=FALSE
)
atacCoords$name<-sub(".+\\|","",rownames(dataAtac))
names(atacCoords)<-rownames(dataAtac)

str(atacCoords)
# 'data.frame':   215920 obs. of  4 variables:
#  $ seqnames: chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start   : int  17234 180634 181207 183557 184247 191249 267757 778474 779599 827304 ...
#  $ end     : int  17733 181133 181706 184056 184746 191748 268256 778973 780098 827803 ...
#  $ name    : chr  "BRCA_2" "BRCA_3" "BRCA_4" "BRCA_5" ...

atacRanges<-GRanges(seqnames=atacCoords$seqnames,
  ranges=IRanges(start=atacCoords$start,
    end=atacCoords$end)
)
names(atacRanges)<-paste(as.character(seqnames(atacRanges)),paste(start(atacRanges),end(atacRanges),sep="-"),sep=":")
#names(atacRanges)<-rownames(dataAtac)

quantile(width(atacRanges))
 #  0%  25%  50%  75% 100% 
 # 500  500  500  500  500 

##no overlaps
length(reduce(atacRanges))
#[1] 215920
length(atacRanges)
#[1] 215920

################################################################################
###Make GRanges for Illumina 450K

identical(rownames(betaAdj),rownames(probesKeep))
#[1] TRUE

head(probesKeep)
#                    id  chr  start    end isBad probeType designType
# cg21870274 cg21870274 chr1  69591  69592 FALSE        cg         II
# cg08258224 cg08258224 chr1 864703 864704 FALSE        cg         II
# cg16619049 cg16619049 chr1 870161 870162 FALSE        cg          I
# cg18147296 cg18147296 chr1 877159 877160 FALSE        cg         II
# cg13938959 cg13938959 chr1 898803 898804 FALSE        cg         II
# cg12445832 cg12445832 chr1 898915 898916 FALSE        cg         II

table(probesKeep$chr)
#  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX  chrY 
# 41098 21206 25633 21624 10540 13204 13119 18627 24356  5280 21984 30550  9354  3378  7222 22315 17646 21279 31170 25024 18194  8605  9718   242 

pAnno<-GRanges(seqnames=probesKeep$chr,
  ranges=IRanges(start=probesKeep$start,
    end=probesKeep$end)
)
names(pAnno)<-rownames(probesKeep)

##give CpGs single coordinate - use 3' in CG - otherwise more complicated overlap calculations
start(pAnno)<-end(pAnno)

pAnno
# GRanges object with 421368 ranges and 0 metadata columns:
#              seqnames    ranges strand
#                 <Rle> <IRanges>  <Rle>
#   cg21870274     chr1     69592      *
#   cg08258224     chr1    864704      *
#   cg16619049     chr1    870162      *
#   cg18147296     chr1    877160      *
#   cg13938959     chr1    898804      *
#          ...      ...       ...    ...
#   cg11021362     chrY  21594791      *
#   cg25914522     chrY  22306397      *
#   cg21106100     chrY  26409390      *
#   cg08265308     chrY  26409404      *
#   cg14273923     chrY  26409766      *
#   -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

################################################################################
##create simple annotation matrix 

identical(rownames(betaAdj),rownames(probesKeep))
#[1] TRUE
identical(rownames(betaAdj),names(pAnno))
#[1] TRUE

annoObj<-data.frame(illuminaID=rownames(probesKeep),stringsAsFactors=FALSE)
rownames(annoObj)<-rownames(probesKeep)

annoObj$chr<-""
annoObj[rownames(annoObj),"chr"]<-probesKeep$chr

annoObj$start<-0
annoObj[rownames(annoObj),"start"]<-probesKeep$start

annoObj$end<-0
annoObj[rownames(annoObj),"end"]<-probesKeep$end

################################################################################
##add "knowngene" mask - rest are for genes with gex

db<-org.Hs.eg.db

p0<-genes(TxDb.Hsapiens.UCSC.hg38.knownGene)               
p0<-p0[as.character(seqnames(p0)) %in% paste0("chr",c(1:22,"X"))]
seqlevels(p0, pruning.mode="coarse") <- paste("chr",c(1:22,"X"),sep="")

#sym<-select(db, keys=p0$gene_id, columns="SYMBOL", keytype="ENTREZID",multiVals="first")
sym<-as.character( mapIds(db, keys=p0$gene_id, column=c("ENTREZID"), keytype="ENTREZID",multiVals="first") )

sym2<-as.character( mapIds(db, keys=p0$gene_id, column=c("SYMBOL"), keytype="ENTREZID",multiVals="first") )

sym2[is.na(sym2)]<-""

sym<-as.data.frame(cbind(ENTREZID=sym,SYMBOL=sym2),stringsAsFactors=FALSE)

str(sym)
# 'data.frame':   25668 obs. of  2 variables:
#  $ ENTREZID: chr  "1" "10" "100" "1000" ...
#  $ SYMBOL  : chr  "A1BG" "NAT2" "ADA" "CDH2" ...

all(sym$ENTREZID==names(p0))
#[1] TRUE

p0$symbol<-sym$SYMBOL

##Gb in genic regions
sum(width(reduce(p0,ignore.strand=TRUE)))/1e9
#[1] 1.515932

p1<-p0
p1<-p1[as.character(seqnames(p1)) %in% paste0("chr",c(1:22,"X"))]
seqlevels(p1, pruning.mode="coarse") <- paste("chr",c(1:22,"X"),sep="")

p2<-promoters(p1,upstream=5000,downstream=0)                                    ###Need to take orientation into account!!
end(p2)[as.vector(strand(p2))=="+"]<-end(p2)[as.vector(strand(p2))=="+"]-500
start(p2)[as.vector(strand(p2))=="-"]<-start(p2)[as.vector(strand(p2))=="-"]+500

p3<-promoters(p1,upstream=0,downstream=5000)
end(p3)[as.vector(strand(p3))=="-"]<-end(p3)[as.vector(strand(p3))=="-"]-500
start(p3)[as.vector(strand(p3))=="+"]<-start(p3)[as.vector(strand(p3))=="+"]+500

p4<-promoters(p1,upstream=500,downstream=500)

##any knownGene overlap?
ol1<-findOverlaps(p1,pAnno)
ol1<-unlist(lapply(split(p1$symbol[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasUCSCknownGeneOverlap<-0
annoObj[names(ol1),"hasUCSCknownGeneOverlap"]<-as.integer(ol1!="")

annoObj$nameUCSCknownGeneOverlap<-""
annoObj[names(ol1),"nameUCSCknownGeneOverlap"]<-ol1

annoObj$numberUCSCknownGeneOverlap<-0
annoObj[names(ol1),"numberUCSCknownGeneOverlap"]<-unlist(lapply(ol1,function(x) length(unlist(strsplit(x,"_")))))

##promoter overlap
ol1<-findOverlaps(p4,pAnno)
ol1<-unlist(lapply(split(names(p4)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasUCSCknownGenePromOverlap<-0
annoObj[names(ol1),"hasUCSCknownGenePromOverlap"]<-as.integer(ol1!="")

##up5000 gene overlap
ol1<-findOverlaps(p2,pAnno)
ol1<-unlist(lapply(split(names(p2)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasUCSCknownGeneUp5kbOverlap<-0
annoObj[names(ol1),"hasUCSCknownGeneUp5kbOverlap"]<-as.integer(ol1!="")

##dn5000 gene overlap
ol1<-findOverlaps(p3,pAnno)
ol1<-unlist(lapply(split(names(p3)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasUCSCknownGeneDn5kbOverlap<-0
annoObj[names(ol1),"hasUCSCknownGeneDn5kbOverlap"]<-as.integer(ol1!="")

##any knownTx promoter overlap
p0<-transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)               
p0<-p0[as.character(seqnames(p0)) %in% paste0("chr",c(1:22,"X"))]
seqlevels(p0, pruning.mode="coarse") <- paste("chr",c(1:22,"X"),sep="")

p1<-p0
p1<-p1[as.character(seqnames(p1)) %in% paste0("chr",c(1:22,"X"))]
seqlevels(p1, pruning.mode="coarse") <- paste("chr",c(1:22,"X"),sep="")

p2<-promoters(p1,upstream=5000,downstream=0)                                    ###Need to take orientation into account!!
end(p2)[as.vector(strand(p2))=="+"]<-end(p2)[as.vector(strand(p2))=="+"]-500
start(p2)[as.vector(strand(p2))=="-"]<-start(p2)[as.vector(strand(p2))=="-"]+500
names(p2)<-1:length(p2)

p3<-promoters(p1,upstream=0,downstream=5000)
end(p3)[as.vector(strand(p3))=="-"]<-end(p3)[as.vector(strand(p3))=="-"]-500
start(p3)[as.vector(strand(p3))=="+"]<-start(p3)[as.vector(strand(p3))=="+"]+500
names(p3)<-1:length(p3)

p4<-promoters(p1,upstream=500,downstream=500)
names(p4)<-1:length(p4)

ol1<-findOverlaps(p4,pAnno)
ol1<-unlist(lapply(split(names(p4)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasUCSCknownTxPromOverlap<-0
annoObj[names(ol1),"hasUCSCknownTxPromOverlap"]<-as.integer(ol1!="")

ol1<-findOverlaps(p2,pAnno)
ol1<-unlist(lapply(split(names(p2)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasUCSCknownTxUp5kbOverlap<-0
annoObj[names(ol1),"hasUCSCknownTxUp5kbOverlap"]<-as.integer(ol1!="")

##dn5000 gene overlap
ol1<-findOverlaps(p3,pAnno)
ol1<-unlist(lapply(split(names(p3)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasUCSCknownTxDn5kbOverlap<-0
annoObj[names(ol1),"hasUCSCknownTxDn5kbOverlap"]<-as.integer(ol1!="")

rm(sym,p0,p1,p2,p3,p4)

##may have some contradictory annotations b/c multiple transcripts may overlap for same or different genes

##add annotation "distal"
annoObj$ucscKnownGeneIsDistal<-0
annoObj[ annoObj$hasUCSCknownGeneUp5kbOverlap==0 &
  annoObj$hasUCSCknownGeneDn5kbOverlap==0 &
  annoObj$hasUCSCknownGenePromOverlap==0 
  ,"ucscKnownGeneIsDistal"]<-1

annoObj$ucscKnownGeneIsGenic<-0
annoObj[ annoObj$hasUCSCknownGeneOverlap!=0 
  ,"ucscKnownGeneIsGenic"]<-1

annoObj$ucscKnownGeneIsDistalNonGenic<-0
annoObj[ annoObj$hasUCSCknownGeneUp5kbOverlap==0 &
  annoObj$hasUCSCknownGeneDn5kbOverlap==0 &
  annoObj$hasUCSCknownGenePromOverlap==0 &
  annoObj$hasUCSCknownGeneOverlap==0 
  ,"ucscKnownGeneIsDistalNonGenic"]<-1

table(annoObj$ucscKnownGeneIsDistal,annoObj$ucscKnownGeneIsDistalNonGenic)
  #        0      1
  # 0 206341      0
  # 1 138020  77007

annoObj$ucscKnownGeneIsGeneBody<-0
annoObj[ annoObj$hasUCSCknownGenePromOverlap==0 &
  annoObj$hasUCSCknownGeneOverlap==1 
  ,"ucscKnownGeneIsGeneBody"]<-1

table(annoObj$ucscKnownGeneIsDistal,annoObj$ucscKnownGeneIsGeneBody)
  #        0      1
  # 0 132249  74092
  # 1  77007 138020
table(annoObj$ucscKnownGeneIsDistalNonGenic,annoObj$ucscKnownGeneIsGeneBody)
  #        0      1
  # 0 132249 212112
  # 1  77007      0

##annotations do not strictly take into account that catogories are fulfulled for same gene -> but ok compromise still

################################################################################
##add more detailed annotations for genes with gex

##do overlaps
p1<-geneCoords

##use ensg+proper symbol as gene identiifer
names(p1)<-paste(geneCoords$gene_id2,geneCoords$gene_name,sep="|")

p1<-p1[as.character(seqnames(p1)) %in% paste0("chr",c(1:22,"X"))]
seqlevels(p1, pruning.mode="coarse") <- paste("chr",c(1:22,"X"),sep="")

p2<-promoters(p1,upstream=5000,downstream=0)                                    ###Need to take orientation into account!!
end(p2)[as.vector(strand(p2))=="+"]<-end(p2)[as.vector(strand(p2))=="+"]-500
start(p2)[as.vector(strand(p2))=="-"]<-start(p2)[as.vector(strand(p2))=="-"]+500

p3<-promoters(p1,upstream=0,downstream=5000)
end(p3)[as.vector(strand(p3))=="-"]<-end(p3)[as.vector(strand(p3))=="-"]-500
start(p3)[as.vector(strand(p3))=="+"]<-start(p3)[as.vector(strand(p3))=="+"]+500

p4<-promoters(p1,upstream=500,downstream=500)

ol1<-findOverlaps(p1,pAnno)
ol1<-unlist(lapply(split(names(p1)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasGeneOverlap<-0
annoObj[names(ol1),"hasGeneOverlap"]<-as.integer(ol1!="")

annoObj$nameGeneOverlap<-""
annoObj[names(ol1),"nameGeneOverlap"]<-ol1

annoObj$numberGeneOverlap<-0
annoObj[names(ol1),"numberGeneOverlap"]<-unlist(lapply(ol1,function(x) length(unlist(strsplit(x,"_")))))

##up5000 gene overlap
ol1<-findOverlaps(p2,pAnno)
ol1<-unlist(lapply(split(names(p2)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasUp5kbOverlap<-0
annoObj[names(ol1),"hasUp5kbOverlap"]<-as.integer(ol1!="")

annoObj$nameUp5kbOverlap<-""
annoObj[names(ol1),"nameUp5kbOverlap"]<-ol1

annoObj$numberUp5kbOverlap<-0
annoObj[names(ol1),"numberUp5kbOverlap"]<-unlist(lapply(ol1,function(x) length(unlist(strsplit(x,"_")))))

##dn5000 gene overlap
ol1<-findOverlaps(p3,pAnno)
ol1<-unlist(lapply(split(names(p3)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasDn5kbOverlap<-0
annoObj[names(ol1),"hasDn5kbOverlap"]<-as.integer(ol1!="")

annoObj$nameDn5kbOverlap<-""
annoObj[names(ol1),"nameDn5kbOverlap"]<-ol1

annoObj$numberDn5kbOverlap<-0
annoObj[names(ol1),"numberDn5kbOverlap"]<-unlist(lapply(ol1,function(x) length(unlist(strsplit(x,"_")))))

##prom gene overlap
ol1<-findOverlaps(p4,pAnno)
ol1<-unlist(lapply(split(names(p4)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="_")))

annoObj$hasPromOverlap<-0
annoObj[names(ol1),"hasPromOverlap"]<-as.integer(ol1!="")

annoObj$namePromOverlap<-""
annoObj[names(ol1),"namePromOverlap"]<-ol1

annoObj$numberPromOverlap<-0
annoObj[names(ol1),"numberPromOverlap"]<-unlist(lapply(ol1,function(x) length(unlist(strsplit(x,"_")))))

##atac overlap
ol1<-findOverlaps(atacRanges,pAnno)
ol1<-unlist(lapply(split(names(atacRanges)[queryHits(ol1)],names(pAnno)[subjectHits(ol1)]),function(x) paste(unique(x),collapse="|")))

annoObj$hasAtacOverlap<-0
annoObj[names(ol1),"hasAtacOverlap"]<-as.integer(ol1!="")

annoObj$nameAtacOverlap<-""
annoObj[names(ol1),"nameAtacOverlap"]<-ol1

annoObj$numberAtacOverlap<-0
annoObj[names(ol1),"numberAtacOverlap"]<-unlist(lapply(ol1,function(x) length(unlist(strsplit(x,"\\|")))))

table(annoObj$hasAtacOverlap,annoObj$hasGeneOverlap)
  #        0      1
  # 0  75244 200010
  # 1  42614 103500

table(annoObj$hasAtacOverlap,annoObj$hasUp5kbOverlap | annoObj$hasDn5kbOverlap)
  #    FALSE   TRUE
  # 0 170853 104401
  # 1  81828  64286

table(annoObj$hasUp5kbOverlap , annoObj$hasDn5kbOverlap)
  #        0      1
  # 0 252681  73308
  # 1  58080  37299 ##some overlap 2 proms or more

table(annoObj$hasDn5kbOverlap,annoObj$hasPromOverlap)
  #        0      1
  # 0 212469  98292
  # 1  86190  24417

table(annoObj$hasUp5kbOverlap,annoObj$hasPromOverlap)
  #        0      1
  # 0 225043 100946
  # 1  73616  21763

table(annoObj$hasAtacOverlap,annoObj$hasPromOverlap)
  #        0      1
  # 0 226106  49148
  # 1  72553  73561

##add annotation "distal"
annoObj$isDistal<-0
annoObj[ annoObj$hasUp5kbOverlap==0 &
  annoObj$hasDn5kbOverlap==0 &
  annoObj$hasPromOverlap==0 
  ,"isDistal"]<-1

annoObj$isGenic<-0
annoObj[ annoObj$hasGeneOverlap!=0 
  ,"isGenic"]<-1

annoObj$isDistalNonGenic<-0
annoObj[ annoObj$hasUp5kbOverlap==0 &
  annoObj$hasDn5kbOverlap==0 &
  annoObj$hasPromOverlap==0 &
  annoObj$hasGeneOverlap==0 
  ,"isDistalNonGenic"]<-1

annoObj$isGeneBody<-0
annoObj[ annoObj$hasPromOverlap==0 &
  annoObj$hasGeneOverlap==1 
  #& annoObj$hasUp5kbOverlap==0 #if this hashed then can be body even if up5k for other gene
  ,"isGeneBody"]<-1

str(annoObj)
# 'data.frame':   421368 obs. of  35 variables:
#  $ illuminaID                   : chr  "cg21870274" "cg08258224" "cg16619049" "cg18147296" ...
#  $ chr                          : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start                        : num  69591 864703 870161 877159 898803 ...
#  $ end                          : num  69592 864704 870162 877160 898804 ...
#  $ hasUCSCknownGeneOverlap      : num  1 0 1 1 0 0 0 0 0 0 ...
#  $ nameUCSCknownGeneOverlap     : chr  "OR4F5" "" "FAM41C" "FAM41C" ...
#  $ numberUCSCknownGeneOverlap   : num  1 0 1 1 0 0 0 0 0 0 ...
#  $ hasUCSCknownGeneUp5kbOverlap : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hasUCSCknownGeneDn5kbOverlap : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ hasUCSCknownTxPromOverlap    : num  0 0 1 1 0 0 0 0 0 0 ...
#  $ hasUCSCknownTxUp5kbOverlap   : num  0 0 1 0 0 0 0 1 1 1 ...
#  $ hasUCSCknownTxDn5kbOverlap   : num  1 1 0 1 0 0 0 0 0 0 ...
#  $ ucscKnownGeneIsDistal        : num  0 1 0 0 1 1 1 1 1 1 ...
#  $ ucscKnownGeneIsGenic         : num  1 0 1 1 0 0 0 0 0 0 ...
#  $ ucscKnownGeneIsDistalNonGenic: num  0 1 0 0 1 1 1 1 1 1 ...
#  $ ucscKnownGeneIsGeneBody      : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ hasGeneOverlap               : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ nameGeneOverlap              : chr  "ENSG00000186092|OR4F5" "" "ENSG00000230368|FAM41C" "" ...
#  $ numberGeneOverlap            : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ hasUp5kbOverlap              : num  0 0 1 0 0 0 0 1 1 1 ...
#  $ nameUp5kbOverlap             : chr  "" "" "ENSG00000234711|TUBB8P11" "" ...
#  $ numberUp5kbOverlap           : num  0 0 1 0 0 0 0 1 1 1 ...
#  $ hasDn5kbOverlap              : num  1 0 0 1 0 0 0 0 0 0 ...
#  $ nameDn5kbOverlap             : chr  "ENSG00000186092|OR4F5" "" "" "ENSG00000234711|TUBB8P11" ...
#  $ numberDn5kbOverlap           : num  1 0 0 1 0 0 0 0 0 0 ...
#  $ hasPromOverlap               : num  0 0 0 1 0 0 0 0 0 0 ...
#  $ namePromOverlap              : chr  "" "" "" "ENSG00000230368|FAM41C" ...
#  $ numberPromOverlap            : num  0 0 0 1 0 0 0 0 0 0 ...
#  $ hasAtacOverlap               : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ nameAtacOverlap              : chr  "" "" "chr1:869670-870169" "" ...
#  $ numberAtacOverlap            : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ isDistal                     : num  0 1 0 0 1 1 1 0 0 0 ...
#  $ isGenic                      : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ isDistalNonGenic             : num  0 1 0 0 1 1 1 0 0 0 ...
#  $ isGeneBody                   : num  1 0 1 0 0 0 0 0 0 0 ...

table(annoObj$chr,annoObj$hasAtacOverlap)
  #           0     1
  # chr1  26067 15031
  # chr10 14843  6363
  # chr11 16482  9151
  # chr12 13967  7657
  # chr13  7819  2721
  # chr14  8868  4336
  # chr15  8619  4500
  # chr16 12708  5919
  # chr17 14455  9901
  # chr18  3378  1902
  # chr19 12321  9663
  # chr2  20496 10054
  # chr20  5465  3889
  # chr21  2158  1220
  # chr22  4383  2839
  # chr3  14211  8104
  # chr4  12459  5187
  # chr5  14651  6628
  # chr6  19948 11222
  # chr7  17743  7281
  # chr8  11745  6449
  # chr9   5890  2715
  # chrX   6336  3382
  # chrY    242     0

##save interim object...
save(annoObj,file=paste0(HOME,"/annotateFeatures/object_noChromFilter_",nrow(annoObj),"x",ncol(annoObj),"_annotations.RData"))

rm(p1,p2,p3,p4,ol1,atacRanges,pAnno,sdAdj)

################################################################################
################################################################################
##add more annotations
################################################################################
################################################################################

################################################################################
################################################################################
##Annotate CpGs

pAnno<-GRanges(seqnames=probesKeep$chr,
  ranges=IRanges(start=probesKeep$start,
    end=probesKeep$end)
)
names(pAnno)<-rownames(probesKeep)

##give CpGs single coordinate - use 3' in CG - otherwise more complicated overlap calculations
start(pAnno)<-end(pAnno)

pAnno
# GRanges object with 421368 ranges and 0 metadata columns:
#              seqnames    ranges strand
#                 <Rle> <IRanges>  <Rle>
#   cg21870274     chr1     69592      *
#   cg08258224     chr1    864704      *
#   cg16619049     chr1    870162      *
#   cg18147296     chr1    877160      *
#   cg13938959     chr1    898804      *
#          ...      ...       ...    ...
#   cg11021362     chrY  21594791      *
#   cg25914522     chrY  22306397      *
#   cg21106100     chrY  26409390      *
#   cg08265308     chrY  26409404      *
#   cg14273923     chrY  26409766      *
#   -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

##Background genome
Hsapiens
#Human genome:
# organism: Homo sapiens (Human)
# provider: UCSC
# provider version: hg38
# release date: Dec. 2013
# release name: Genome Reference Consortium GRCh38
# 455 sequences:
#   chr1                    chr2                    chr3                   
#   chr4                    chr5                    chr6                   
#   chr7                    chr8                    chr9                   
#   chr10                   chr11                   chr12                  
#   chr13                   chr14                   chr15                  
#   ...                     ...                     ...                    
#   chrUn_KI270744v1        chrUn_KI270745v1        chrUn_KI270746v1       
#   chrUn_KI270747v1        chrUn_KI270748v1        chrUn_KI270749v1       
#   chrUn_KI270750v1        chrUn_KI270751v1        chrUn_KI270752v1       
#   chrUn_KI270753v1        chrUn_KI270754v1        chrUn_KI270755v1       
#   chrUn_KI270756v1        chrUn_KI270757v1                               
# (use 'seqnames()' to see all the sequence names, use the '$' or '[[' operator
# to access a given sequence)

##Gene models
table(geneCoords$gene_type)
 #          3prime_overlapping_ncrna                          antisense                          IG_C_gene                    IG_C_pseudogene 
 #                                29                               5565                                 14                                  9 
 #                         IG_D_gene                          IG_J_gene                    IG_J_pseudogene                          IG_V_gene 
 #                                37                                 18                                  3                                142 
 #                   IG_V_pseudogene                            lincRNA                       macro_lncRNA                              miRNA 
 #                               180                               7656                                  1                               4093 
 #                          misc_RNA                            Mt_rRNA                            Mt_tRNA                         non_coding 
 #                              2298                                  2                                 22                                  3 
 #            polymorphic_pseudogene               processed_pseudogene               processed_transcript                     protein_coding 
 #                                59                              10304                                484                              19814 
 #                        pseudogene                           ribozyme                               rRNA                             scaRNA 
 #                                36                                  8                                544                                 49 
 #                    sense_intronic                  sense_overlapping                             snoRNA                              snRNA 
 #                               920                                197                                961                               1896 
 #                              sRNA                                TEC                          TR_C_gene                          TR_D_gene 
 #                                20                               1045                                  5                                  3 
 #                         TR_J_gene                    TR_J_pseudogene                          TR_V_gene                    TR_V_pseudogene 
 #                                73                                  4                                106                                 30 
 #  transcribed_processed_pseudogene     transcribed_unitary_pseudogene transcribed_unprocessed_pseudogene    translated_processed_pseudogene 
 #                               443                                  1                                663                                  1 
 # translated_unprocessed_pseudogene                 unitary_pseudogene             unprocessed_pseudogene                           vaultRNA 
 #                                 1                                169                               2574                                  1 

table(as.character(seqnames(geneCoords)))
 # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrM  chrX  chrY 
 # 5397  2306  3381  3047  1383  2289  2247  2597  3111  1206  2997  4150  1436   880  1385  3163  2633  2993  3001  2980  2444  2350    37  2476   594 

##choose all promoters
tssCord<-promoters(geneCoords,upstream=0,downstream=1)
seqlevels(tssCord, pruning.mode="coarse") <- paste("chr",c(1:22,"X"),sep="")

length(tssCord)
#[1] 59852

##seems pretty unique
length(reduce(tssCord))
#[1] 59743

##CpGs
probeAnno
# GRanges object with 421368 ranges and 0 metadata columns:
#              seqnames    ranges strand
#                 <Rle> <IRanges>  <Rle>
#   cg21870274     chr1     69592      *
#   cg08258224     chr1    864704      *
#   cg16619049     chr1    870162      *
#   cg18147296     chr1    877160      *
#   cg13938959     chr1    898804      *
#          ...      ...       ...    ...
#   cg11021362     chrY  21594791      *
#   cg25914522     chrY  22306397      *
#   cg21106100     chrY  26409390      *
#   cg08265308     chrY  26409404      *
#   cg14273923     chrY  26409766      *
#   -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

table(seqnames(probeAnno))
#  chr1 chr22  chr2  chr3  chr4  chr5  chr6  chr7  chr8  chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21  chrX  chrY 
# 41098  7222 30550 22315 17646 21279 31170 25024 18194  8605 21206 25633 21624 10540 13204 13119 18627 24356  5280 21984  9354  3378  9718   242 

all.equal(rownames(annoObj),rownames(betaAdj))
#[1] TRUE
all.equal(rownames(annoObj),names(pAnno))
#[1] TRUE

length(pAnno)
#[1] 421368

################################################################################
################################################################################
##derive CpG/OE for all positions

##Saxonov et al. 2006 PNAS 
  #CpG/OE = obs(CpG) / ( GC content/2 )^2
  #3kb centered on TSS
    #CpG hi = 0.35

##saxonov windows 
table(width(pAnno))
#      1 
# 421368 

##use 3´coordinate 
seqCpgSaxonov<-getSeq(Hsapiens,names=as.character(seqnames(pAnno)),start=start(pAnno)-1500,end=start(pAnno)+1500)

unique(width(seqCpgSaxonov))
#[1] 3001

seqPromSaxonov<-getSeq(Hsapiens,names=as.character(seqnames(tssCord)),start=start(tssCord)-1500,end=end(tssCord)+1500)

unique(width(seqPromSaxonov))
#[1] 3001

##Saxonov et al. CpG O/E calculation
cpgOE<-function(x1,x2,seqL=3001) as.vector((x1/seqL) / ( (x2/seqL)/2 )^2)

##promoters
cgCount<-vcountPattern("CG",seqPromSaxonov)
gcCont<-letterFrequency(seqPromSaxonov,letters=c("GC"))
promOEsaxonov<-cpgOE(cgCount,gcCont)

range(promOEsaxonov)
#[1] 0.02705541 1.73724799

length(promOEsaxonov)
#[1] 59852

promClassSaxonov<-c("LCG","HCG")[ as.integer(promOEsaxonov >= .35)+1 ]

table(promClassSaxonov)
# promClassSaxonov
#   HCG   LCG 
# 24328 35524 

##all Cpgs
cgCount<-vcountPattern("CG",seqCpgSaxonov)
gcCont<-letterFrequency(seqCpgSaxonov,letters=c("GC"))
cpgOEsaxonov<-cpgOE(cgCount,gcCont)

length(cpgOEsaxonov)
#[1] 421368

cpgClassSaxonov<-c("LCG","HCG")[ as.integer(cpgOEsaxonov >= .35)+1 ]

table(cpgClassSaxonov)
# cpgClassSaxonov
#    HCG    LCG 
# 259394 161974 

ls()[grep("saxonov|Saxonov",ls())]
# [1] "cpgClassSaxonov"  "cpgOEsaxonov"     "promClassSaxonov" "promOEsaxonov"   
# [5] "seqCpgSaxonov"    "seqPromSaxonov"  

rm(seqCpgSaxonov,seqPromSaxonov)

names(cpgClassSaxonov)<-names(pAnno)
names(cpgOEsaxonov)<-names(pAnno)

names(promClassSaxonov)<-names(tssCord)
names(promOEsaxonov)<-names(tssCord)

objectsSaxonov<-list(cpgClassSaxonov=cpgClassSaxonov,
  cpgOEsaxonov=cpgOEsaxonov,
  promClassSaxonov=promClassSaxonov,
  promOEsaxonov=promOEsaxonov)

##save objects
save(objectsSaxonov,file=paste0(HOME,"/annotateFeatures/object_saxonovCpGStats.RData"))

##define some palettes
pal<-c("#551A8B","#FF0000","#2E8B57",brewer.pal(12,"Set3"))
pal2<-brewer.pal(5,"Set1")
##gexclass
pal3<-brewer.pal(5,"Dark2")
##up/down regulation - colors
pal4<-c(0,"#2E8B57","#FF4500")

##plot promoter CpG/OE
pdf(file=paste0(HOME,"/annotateFeatures/hist_cpgOeSaxonovForPromoters.pdf"),width=12,height=12,useDingbats=FALSE)
par(fig=c(0,1,0,1),mar=c(5.1,4.1,4.1,2.1),font=2,font.sub=2,font.lab=2,font.axis=2,las=1)
hist(promOEsaxonov,breaks=seq(0,max(promOEsaxonov),length.out=100),freq=F,axes=F,xlab="CpG O/E",main="GDC GENCODE v22 promoters",xlim=c(0,1),las=1)
#lines(density(promOEsaxonov),col=1,lwd=3.5)
findMaxima<-function(x,nMax=2,granularity=512,bw=.02,adjust=1,from=0,to=1,kernel="gauss",outPutX=T) {
  aa<-density(x,kernel=kernel,n=granularity,bw=bw,adjust=adjust,from=from,to=to,na.rm=T)
  res<-vector(length=length(aa$x))
  for(i in 2:length(res[-1])) {
    res[i]<-aa$y[i]>aa$y[i-1] & aa$y[i]>aa$y[i+1]
  }  #$x for coord, $y for density
  if(outPutX) {
  aa$x[res][order(aa$y[res],decreasing=T)][1:nMax]
  } else {
  aa$y[res][order(aa$y[res],decreasing=T)][1:nMax]
  }
} ##2 main peaks probably near zero and one
d<-findMaxima(promOEsaxonov,outPutX=F)
d<-d[order(d)]
d2<-density(promOEsaxonov[promClassSaxonov=="HCG"],bw=.02)
d2$y<-(d2$y/max(d2$y))*d[1]
lines(d2,col=pal2[3],lwd=3)
abline(v=findMaxima(promOEsaxonov),lwd=3,lty=2)
d3<-density(promOEsaxonov[promClassSaxonov=="LCG"],bw=.02)
d3$y<-(d3$y/max(d3$y))*d[2]
lines(d3,col=pal2[4],lwd=3)
legend(x=.65,y=max(density(promOEsaxonov)$y)+.1,legend=c("GENCODE HCP promoters","GENCODE LCG promoters"),col=pal2[3:4],pch=16,bty="n",cex=1.5)
axis(1,lwd=2)
axis(2,lwd=2)
dev.off()

rm(d,d2,d3,findMaxima)

##plot CpG CpG/OE
pdf(file=paste0(HOME,"/annotateFeatures/hist_cpgOeSaxonovForCpGs.pdf"),width=12,height=12,useDingbats=FALSE)
par(fig=c(0,1,0,1),mar=c(5.1,4.1,4.1,2.1),font=2,font.sub=2,font.lab=2,font.axis=2,las=1)
hist(cpgOEsaxonov,breaks=seq(0,max(cpgOEsaxonov),length.out=100),freq=F,axes=F,xlab="CpG O/E",main="CpGs",xlim=c(0,1),las=1)
lines(density(cpgOEsaxonov),lwd=2)
#lines(density(cpgOEsaxonov),col=1,lwd=3.5)
findMaxima<-function(x,nMax=2,granularity=512,bw=.02,adjust=1,from=0,to=1,kernel="gauss",outPutX=T) {
  aa<-density(x,kernel=kernel,n=granularity,bw=bw,adjust=adjust,from=from,to=to,na.rm=T)
  res<-vector(length=length(aa$x))
  for(i in 2:length(res[-1])) {
    res[i]<-aa$y[i]>aa$y[i-1] & aa$y[i]>aa$y[i+1]
  }  #$x for coord, $y for density
  if(outPutX) {
  aa$x[res][order(aa$y[res],decreasing=T)][1:nMax]
  } else {
  aa$y[res][order(aa$y[res],decreasing=T)][1:nMax]
  }
} ##2 main peaks probably near zero and one
d<-findMaxima(promOEsaxonov,outPutX=T)
#(d<-findMaxima(cpgOEsaxonov,outPutX=F,nMax=3)  ) ##need 3 peaks b/c peak[2] is intermed
#d2<-density(cpgOEsaxonov,bw=.02)
#lines(d2,col=1,lwd=3)
#abline(v=d,lwd=3,lty=2)
#d3<-density(cpgOEsaxonov[CpGClassSaxonov=="LCG"],bw=.02)
#d3$y<-(d3$y/max(d3$y))*d[1]
#lines(d3,col=pal2[4],lwd=3)
abline(v=d,lwd=3,lty=2,col=pal2[3:4])
legend(x=.65,y=max(density(cpgOEsaxonov)$y)+.1,legend=c("GENCODE HCP promoters","GENCODE LCG promoters"),col=pal2[3:4],pch=16,bty="n",cex=1.5)
axis(1,lwd=2)
axis(2,lwd=2)
dev.off()

rm(d,findMaxima)

rm(cpgOE,pal,pal2,pal3,pal4)

################################################################################
################################################################################
###Weber classification

##Weber class (Weber 2007)
  ## 500bp sliding window analyses w/ 5bp offset
  ## Promoters from -1300 to +200
  ## High O/E>.75 & GC>.55
  ## Low  O/E<.48
  ## else intermediate

##regions need padding for window-based calculations
  ##Weber used -200 and +200 pad for windows

##Weber windows CpGs
  ##these are not really promoters but could be treated as such..
seqCpgWeber<-getSeq(Hsapiens,names=as.character(seqnames(pAnno)),start=start(pAnno)-749-200,end=start(pAnno)+750+200)
                                                                                            
unique(width(seqCpgWeber))
#[1] 1900
   
length(seqCpgWeber)
#[1] 421368

tssCordWeber<-promoters(tssCord,upstream=1300,downstream=200)

seqPromWeber<-getSeq(Hsapiens,names=as.character(seqnames(tssCordWeber)),start=start(tssCordWeber)-200,end=end(tssCordWeber)+200)

unique(width(seqPromWeber))
#[1] 1900

################################################################################
################################################################################
###Weber for all CpGs
  
##define function for splitting sequence
  ##281 windows per sequence
weberString<-function(seq,windW=500,offset=5,width=1900) {
  a<-unlist(strsplit(as.character(seq),""),use.names=F)
  myseq<-seq(1,(width-windW)+1,by=offset)
  mywind<-lapply(myseq,function(x) paste(a[x:(x+windW-1)],collapse="") )
  return( unlist(mywind,use.names=F) )
}

length(seqCpgWeber)
#[1] 421368

system.time( replicate(10,lapply(seqCpgWeber[1],weberString)) ) /10
  #    user  system elapsed 
  # 0.035   0.000   0.035 

0.035*length(seqCpgWeber)/3600
#[1] 4.096633

##check if some N's present
test<-letterFrequency(seqCpgWeber,letters="ATCG")
table(test)
# test
#   1617   1890   1899   1900 
#      1      2      8 421357 

test<-letterFrequency(seqCpgWeber,letters="ATCGN")
table(test)
# test
#   1900 
# 421368 

##Not super many
rm(test)

##running "weberString" on entire object would use enormous amounts of memory -> run line-by-line and write to file
ff<-gzfile(paste0(HOME,"/annotateFeatures/cpgClassWeber.txt.gz"),"w")
gg<-gzfile(paste0(HOME,"/annotateFeatures/cpgOEWeber.txt.gz"),"w")
hh<-gzfile(paste0(HOME,"/annotateFeatures/cpgGCWeber.txt.gz"),"w")

for (i in 1:length(seqCpgWeber) ) {                #length(seqCpgWeber)
  # if(i %% 10000 == 0) { cat(".")}
  # if(i %% 100000 == 0) { cat("\n")}

  weberStringCpg<-lapply(seqCpgWeber[i],weberString)

  weberStringCpg<-do.call("c",weberStringCpg)

  weberStringCpg<-DNAStringSet(weberStringCpg)

  cgCount<-vcountPattern("CG",weberStringCpg)
  gcCont<-letterFrequency(weberStringCpg,letters=c("GC"),OR=0)

  weberOEcpg<- ( cgCount *500 ) / ( gcCont[,1] * gcCont[,2] )
  weberGCcpg<- ( gcCont[,1] + gcCont[,2] )/500
  rm(cgCount,gcCont)

  weberOEcpg<-round(weberOEcpg,4)
  weberGCcpg<-round(weberGCcpg,4)

  ##Weber Class
  #ICP is default density
  weberClassCpg<-"ICP" 
    if( all( weberOEcpg < .48,na.rm=T) ) {
      weberClassCpg<-"LCP"
    }
    if( any( weberOEcpg > .75 & weberGCcpg > .55,na.rm=T) ) {
      weberClassCpg<-"HCP"
    }
  #write class
  writeLines(weberClassCpg,con=ff)
  #write class
  writeLines(paste(weberOEcpg,collapse="\t"),con=gg)
  #write class
  writeLines(paste(weberGCcpg,collapse="\t"),con=hh)
  if(i%%5000==0) { cat(".") }
  if(i%%50000==0) { cat("\n") }

}

rm(i,weberStringCpg,weberOEcpg,weberGCcpg,weberClassCpg)

close(ff)
close(gg)
close(hh)

rm(ff,gg,hh)

weberClassCpg<-readLines(paste0(HOME,"/annotateFeatures/cpgClassWeber.txt.gz"))

table(weberClassCpg)
# weberClassCpg
#    HCP    ICP    LCP 
# 173901 130372 117095 

##calculate average CpG O/E and GC for weber
cpgOEweber<-vcountPattern("CG",seqCpgWeber)
cpgGCweber<-letterFrequency(seqCpgWeber,letters=c("GC"),OR=0)

cpgOEweber<- (cpgOEweber*1900) / ( cpgGCweber[,1]*cpgGCweber[,2] )
cpgGCweber<-  (cpgGCweber[,1]+cpgGCweber[,2])/1900

##do plots
pdf(file=paste0(HOME,"/annotateFeatures/hist_cpgOeWeberForCpgs.pdf"),width=12,height=12,useDingbats=FALSE)
par(fig=c(0,.5,0,.9),mar=c(5.1,4.1,4.1,2.1),font=2,font.sub=2,font.lab=2,font.axis=2,las=1)
hist(cpgOEweber,breaks=seq(0,max(cpgOEweber),length.out=71),freq=T,axes=F,ylab="",xlab="CpG O/E",main="",xlim=c(0,1.2),ylim=c(0,15e3),las=1)
axis(1,lwd=2)
axis(2,lwd=2)
par(fig=c(.5,1,.4,.9),mar=c(2.1,4.1,4.1,2.1),new=T)
hist(cpgOEweber[weberClassCpg=="LCP"],ylab="",breaks=seq(0,max(cpgOEweber),length.out=71),freq=T,axes=F,xlab="",main="",xlim=c(0,1.2),ylim=c(0,15e3),las=1,col="firebrick1")
hist(cpgOEweber[weberClassCpg=="HCP"],add=T,col="dodgerblue",freq=T,breaks=seq(0,max(cpgOEweber),length.out=71))
axis(2,lwd=2)
par(fig=c(.5,1,0,.5),mar=c(5.1,4.1,4.1,2.1),new=T)
hist(cpgOEweber[weberClassCpg=="ICP"],ylab="",breaks=seq(0,max(cpgOEweber),length.out=71),freq=T,axes=F,xlab="",main="",xlim=c(0,1.2),ylim=c(0,15e3),las=1,col="darkgreen")
axis(1,lwd=2)
axis(2,lwd=2)
dev.off()

################################################################################
################################################################################
###Weber for Promoters

##define function for splitting sequence
weberString<-function(seq,windW=500,offset=5,width=1900) { 
  a<-unlist(strsplit(as.character(seq),""),use.names=F)
  myseq<-seq(1,(width-windW)+1,by=offset)
  mywind<-lapply(myseq,function(x) paste(a[x:(x+windW-1)],collapse="") )
  return( unlist(mywind,use.names=F) )
}

##split each sequence into 500 subseqs
weberStringProm<-lapply(seqPromWeber,weberString)
 
unique(unlist(lapply(weberStringProm,width)))
#[1] 500 
unique(unlist(lapply(weberStringProm,length)))
#[1] 281

##linearize into one long string with 281 substrings from each main string
weberStringProm<-do.call("c",weberStringProm)

length(weberStringProm)
#[1] 16818412

weberStringProm<-DNAStringSet(weberStringProm) 

cgCount<-vcountPattern("CG",weberStringProm)
gcCont<-letterFrequency(weberStringProm,letters=c("GC"),OR=0)

weberOEprom<- ( cgCount * 500 ) / ( gcCont[,1] * gcCont[,2] )
weberGCprom<- ( gcCont[1,] + gcCont[,2] )/500
rm(cgCount,gcCont)

##split and get Weber class  length(seqEcrWeber)
splitOE<-split(weberOEprom,rep(1:length(seqPromWeber),each=281))
splitGC<-split(weberGCprom,rep(1:length(seqPromWeber),each=281))
rm(weberOEprom,weberGCprom)

##Weber Class
weberClassProm<-vector(length=length(splitOE))

for (i in 1:length(splitOE)) {
  weberClassProm[i]<-"ICP"
  if( all( splitOE[[i]]<.48,na.rm=T ) ) {
    weberClassProm[i]<-"LCP"
  }
  if( any( splitOE[[i]]>.75 & splitGC[[i]]>.55 , na.rm=T ) ) {
    weberClassProm[i]<-"HCP"
  }
} ; rm(i)

promOEweber<-vcountPattern("CG",seqPromWeber)
promGCweber<-letterFrequency(seqPromWeber,letters=c("GC"),OR=0)

promOEweber<- (promOEweber*1900) / ( promGCweber[,1]*promGCweber[,2] )
promGCweber<-  (promGCweber[,1]+promGCweber[,2])/1900

table(weberClassProm)
# weberClassProm
#   HCP   ICP   LCP 
# 18022 14039 27791 

pdf(file=paste0(HOME,"/annotateFeatures/hist_cpgOeWeberForPromoters.pdf"),width=12,height=12,useDingbats=FALSE)
par(fig=c(0,.5,0,.9),mar=c(5.1,4.1,4.1,2.1),font=2,font.sub=2,font.lab=2,font.axis=2,las=1)
hist(promOEweber,breaks=seq(0,max(promOEweber),length.out=101),ylim=c(0,5750),freq=T,axes=F,ylab="",xlab="CpG O/E",main="",xlim=c(0,1.5),las=1)
axis(1,lwd=2)
axis(2,lwd=2)
par(fig=c(.5,1,.4,.9),mar=c(2.1,4.1,4.1,2.1),new=T)
hist(promOEweber[weberClassProm=="LCP"],ylab="",breaks=seq(0,max(promOEweber),length.out=101),freq=T,axes=F,xlab="",main="",xlim=c(0,1.5),ylim=c(0,5750),las=1,col="firebrick1")
hist(promOEweber[weberClassProm=="HCP"],add=T,col="dodgerblue",freq=T,breaks=seq(0,max(promOEweber),length.out=101))
axis(2,lwd=2)
par(fig=c(.5,1,0,.5),mar=c(5.1,4.1,4.1,2.1),new=T)
hist(promOEweber[weberClassProm=="ICP"],ylab="",breaks=seq(0,max(promOEweber),length.out=101),freq=T,axes=F,xlab="",main="",ylim=c(0,5750),xlim=c(0,1.5),las=1,col="darkgreen")
axis(1,lwd=2)
axis(2,lwd=2)
legend(x=.5,y=2.75,legend=c("HCP promoter","ICP promoter","LCP promoter"),col=c("dodgerblue","darkgreen","firebrick1"),pch=16,bty="n",cex=1.5)
dev.off()

rm(splitOE,splitGC)

ls()[grep("weber|Weber",ls())]
 # [1] "cpgGCweber"      "cpgOEweber"      "promGCweber"     "promOEweber"    
 # [5] "seqCpgWeber"     "seqPromWeber"    "tssCordWeber"    "weberClassCpg"  
 # [9] "weberClassProm"  "weberString"     "weberStringProm"

##remove unneeded objects
rm(seqCpgWeber,seqPromWeber,tssCordWeber,weberString,weberStringProm)

gc()

names(weberClassCpg)<-names(pAnno)
names(cpgOEweber)<-names(pAnno)
names(cpgGCweber)<-names(pAnno)

names(weberClassProm)<-names(tssCord)
names(promOEweber)<-names(tssCord)
names(promGCweber)<-names(tssCord)

objectsWeber<-list(cpgClassWeber=weberClassCpg,
  cpgOEweber=cpgOEweber,
  cpgGCweber=cpgGCweber,
  promClassWeber=weberClassProm,
  promOEweber=promOEweber,
  promGCweber=promGCweber)

save(objectsWeber,file=paste0(HOME,"/annotateFeatures/object_weberCpGStats.RData"))

################################################################################
################################################################################
###Cleanup

ls()[grep("saxonov|Saxonov",ls())]
#[1] "cpgClassSaxonov"  "cpgOEsaxonov"     "objectsSaxonov"   "promClassSaxonov"
#[5] "promOEsaxonov"  

ls()[grep("weber|Weber",ls())]
#[1] "cpgGCweber"     "cpgOEweber"     "objectsWeber"   "promGCweber"   
#[5] "promOEweber"    "weberClassCpg"  "weberClassProm"

length(cpgClassSaxonov)
#[1] 421368
length(promClassSaxonov)
#[1] 59852

length(weberClassCpg)
#[1] 421368
length(weberClassProm)
#[1] 59852

length(tssCord)
#[1] 59852
length(pAnno)
#[1] 421368

head(names(pAnno))
# [1] "cg21870274" "cg08258224" "cg16619049" "cg18147296" "cg13938959" "cg12445832"

head(names(tssCord))
#[1] "ENSG00000223972.5" "ENSG00000238009.5" "ENSG00000230415.1" "ENSG00000236335.1" "ENSG00000213842.2" "ENSG00000227337.1"

##remove unneeded objects
rm(cpgClassSaxonov,cpgOEsaxonov,promClassSaxonov,promOEsaxonov)

rm(cpgGCweber,cpgOEweber,promGCweber,promOEweber,weberClassCpg,weberClassProm)

rm(pAnno,tssCord)

gc()

##Add annotation for Weber and Saxonov class as well as CpG O/E
length(intersect(annoObj$illuminaID,names(objectsWeber$cpgClassWeber)))
#[1] 421368

all.equal(annoObj$illuminaID,names(objectsWeber$cpgClassWeber))
#[1] TRUE

all(intersect(annoObj$illuminaID,names(objectsWeber$cpgClassWeber))==annoObj$illuminaID)
#[1] TRUE

ol<-intersect(annoObj$illuminaID,names(objectsWeber$cpgClassWeber))

identical(ol,annoObj$illuminaID)
#[1] TRUE

annoObj$weberOE<-objectsWeber$cpgOEweber[ol]
annoObj$weberClass<-objectsWeber$cpgClassWeber[ol]

ol<-intersect(annoObj$illuminaID,names(objectsSaxonov$cpgClassSaxonov))

identical(ol,annoObj$illuminaID)
#[1] TRUE

annoObj$saxonovOE<-objectsSaxonov$cpgOEsaxonov[ol]
annoObj$saxonovClass<-objectsSaxonov$cpgClassSaxonov[ol]

table(annoObj$weberClass,annoObj$saxonovClass)
  #        HCG    LCG
  # HCP 172695   1206
  # ICP  81434  48938
  # LCP   5265 111830

cor(annoObj$weberOE,annoObj$saxonovOE)
#[1] 0.9639102

table(annoObj$weberClass,annoObj$hasAtacOverlap)
  #          0      1
  # HCP  74163  99738
  # ICP 102551  27821
  # LCP  98540  18555

table(annoObj$weberClass,annoObj$hasPromOverlap)
  #          0      1
  # HCP  81581  92320
  # ICP 111843  18529
  # LCP 105235  11860

rm(ol)

gc()

save(annoObj,file=paste0(HOME,"/annotateFeatures/object_noChromFilter_",nrow(annoObj),"x",ncol(annoObj),"_expandedAnnotations.RData"))

################################################################################
################################################################################
###Calculate simple CpG/bp in 300 bp window

pAnno<-GRanges(seqnames=probesKeep$chr,
  ranges=IRanges(start=probesKeep$start,
    end=probesKeep$end)
)
names(pAnno)<-rownames(probesKeep)

identical(names(pAnno),rownames(annoObj))
#[1] TRUE

##give CpGs single coordinate - use 3' in CG - otherwise more complicated overlap calculations
start(pAnno)<-end(pAnno)

##base on calculation in Elife modeling methylation and CpG density..
  ##can replicate model empirically?
    ##https://elifesciences.org/articles/04094
      ##CpG content calculated in 300 bp window

##modeled methylation as exp function y = 100 / ( 1 + e^-b(x-c) ) with parameters b=-0.337 and c=6.917

##use 3´coordinate - because CpG is 2 bases do 3'-coordinate minus 150 and plus 149 bp
seqCpg<-getSeq(Hsapiens,names=as.character(seqnames(pAnno)),start=start(pAnno)-150,end=start(pAnno)+149)

unique(width(seqCpg))
#[1] 300

##Simple CpG count
cgCount<-vcountPattern("CG",seqCpg)

quantile(cgCount)
  # 0%  25%  50%  75% 100% 
  #  0    5   12   20   69 

length(cgCount)
#[1] 421368

annoObj$cgCount300Bp<- cgCount

annoObj$cgPer100Bp<- round( (cgCount/300) *100 ,4 )

cor(annoObj$cgPer100Bp,annoObj$saxonovOE)
#[1] 0.6999895

##plot subset of CpGs vs cg/100bp and compare to empirical model
mod<-function(x) 100 / ( 1+ exp(1)^(.337*(x-6.917) ) )

set.seed(12345)
pa<-sample(1:nrow(annoObj),5000)

tiff(,file=paste0(HOME,"/annotateFeatures/plot_cpgPerBpVsAvgMethylation.tiff"),width=8*300,height=8*300,res=300,unit="px",compression="lzw")
par(mar=c(5.1,4.1,4.1,2.1),font=2,font.sub=2,font.lab=2,font.axis=2,las=1)
plot(seq(0,100,by=1),mod(seq(0,100,by=1)),xlim=c(2,16),ylim=c(0,100),pch=16,axes=F,xlab="CpG per 100 bp",ylab="mean methylation",las=1)
lines(x=seq(0,100,by=1),y=mod(seq(0,100,by=1)),lwd=3,lty=2,col=1)
points(x=annoObj$cgPer100Bp[pa],y=rowMeans(betaAdj[pa,])*100,,col="firebrick1",pch=16,cex=.4)
lines(lowess((100*rowMeans(betaAdj[pa,]))~annoObj$cgPer100Bp[pa],f=2/5),col="firebrick1",lwd=3,lty=1)
lines(lowess((100*rowMeans(betaFinal[pa,]))~annoObj$cgPer100Bp[pa],f=2/5),col="firebrick2",lwd=3,lty=1)
lines(lowess((100*rowMeans(betaNorm[pa,]))~annoObj$cgPer100Bp[pa],f=2/5),col="firebrick",lwd=3,lty=1)
legend("topright",legend=c("Krebs et al. 2014 model","lowess adjusted beta","lowess original beta","lowess inferred normal"),
  lty=c(2,1,1,1),lwd=3,col=c("black","firebrick1","firebrick2","firebrick"),bty="n"
)
axis(1,lwd=2,font=2)
axis(2,lwd=2,font=2,las=2)

dev.off()
rm(mod,pa)

################################################################################
################################################################################
###Add CGI/Shore/Ocean

download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz",paste0(HOME,"/annotateFeatures/ucsc_hg38_cpgIslandExt.txt.gz"),mode="wb")

as.character( md5sum(paste0(HOME,"/annotateFeatures/ucsc_hg38_cpgIslandExt.txt.gz") ) )
#[1] "dac1e93e83aec3cc13c1b6ff5be5315e"

islands_df <- read.table(paste0(HOME,"/annotateFeatures/ucsc_hg38_cpgIslandExt.txt.gz"),sep = "\t",comment.char = "#",header = F,stringsAsFactors = FALSE, fill=TRUE)

str(islands_df)
#'data.frame':   31144 obs. of  11 variables:
# $ V1 : int  585 586 586 587 587 588 589 589 590 591 ...
# $ V2 : chr  "chr1" "chr1" "chr1" "chr1" ...
# $ V3 : int  28735 135124 199251 368792 381172 491107 597839 609358 778604 827036 ...
# $ V4 : int  29737 135563 200121 370063 382185 491546 598734 611269 779167 828065 ...
# $ V5 : chr  "CpG: 111" "CpG: 30" "CpG: 104" "CpG: 99" ...
# $ V6 : int  1002 439 870 1271 1013 439 895 1911 563 1029 ...
# $ V7 : int  111 30 104 99 84 29 94 171 60 111 ...
# $ V8 : int  731 295 643 777 734 295 570 1405 385 657 ...
# $ V9 : num  22.2 13.7 23.9 15.6 16.6 13.2 21 17.9 21.3 21.6 ...
# $ V10: num  73 67.2 73.9 61.1 72.5 67.2 63.7 73.5 68.4 63.8 ...
# $ V11: num  0.85 0.64 0.89 0.84 0.64 0.62 1.04 0.67 0.92 1.09 ...

colnames(islands_df)<-c("bin",
"chrom",
"chromStart",
"chromEnd",
"name",
"length",
"cpgNum",
"gcNum",
"perCpg",
"perGc",
"obsExp"
)

islands_df <- islands_df[islands_df$chrom %in% paste0("chr",c(1:22,"X","Y")),]

table(islands_df$chrom)
#  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
#  2535  1165  1371  1211   612   796   781  1501  1650   509  2492  1682   847 
# chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX  chrY 
#   446   737  1156  1025  1236  1255  1582  1040  1237   889   194 

##de-BED
islands_df$chromStart<-islands_df$chromStart+1

all(islands_df$length == (islands_df$chromEnd-islands_df$chromStart+1))
#[1] TRUE

write.table(islands_df,file=paste0(HOME,"/annotateFeatures/ucsc_hg38_cpgIslandExt_cannonicalChrFilter_deBed.txt"),sep = "\t",quote=F)

str(islands_df)
# 'data.frame':   27949 obs. of  11 variables:
#  $ bin       : int  585 586 586 587 587 588 589 589 590 591 ...
#  $ chrom     : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ chromStart: num  28736 135125 199252 368793 381173 ...
#  $ chromEnd  : int  29737 135563 200121 370063 382185 491546 598734 611269 779167 828065 ...
#  $ name      : chr  "CpG: 111" "CpG: 30" "CpG: 104" "CpG: 99" ...
#  $ length    : int  1002 439 870 1271 1013 439 895 1911 563 1029 ...
#  $ cpgNum    : int  111 30 104 99 84 29 94 171 60 111 ...
#  $ gcNum     : int  731 295 643 777 734 295 570 1405 385 657 ...
#  $ perCpg    : num  22.2 13.7 23.9 15.6 16.6 13.2 21 17.9 21.3 21.6 ...
#  $ perGc     : num  73 67.2 73.9 61.1 72.5 67.2 63.7 73.5 68.4 63.8 ...
#  $ obsExp    : num  0.85 0.64 0.89 0.84 0.64 0.62 1.04 0.67 0.92 1.09 ...

##create islands object
cgiObj<-GRanges(seqnames=islands_df$chrom,ranges=IRanges(start=islands_df$chromStart,end=islands_df$chromEnd))

##create shores object
shores<-GRanges(seqnames=islands_df$chrom,ranges=IRanges(start=islands_df$chromStart-2000,end=islands_df$chromEnd+2000))

shores<-reduce(shores)
shores<-setdiff(shores,cgiObj)

findOverlaps(shores,cgiObj)
# Hits object with 0 hits and 0 metadata columns:
#    queryHits subjectHits
#    <integer>   <integer>
#   -------
#   queryLength: 50647 / subjectLength: 27949

quantile(width(cgiObj))
  #  0%   25%   50%   75%  100% 
  # 201   328   575   969 45712 

quantile(width(shores))
  # 0%  25%  50%  75% 100% 
  # 64 2000 2000 2000 3998 

##get probe coordinates
# pAnno<-GRanges(seqnames=probesKeep$chr,
#   ranges=IRanges(start=probesKeep$start,
#     end=probesKeep$end)
# )
# names(pAnno)<-rownames(probesKeep)

# identical(names(pAnno),rownames(annoObj))
# #[1] TRUE

# ##give CpGs single coordinate - use 3' in CG - otherwise more complicated overlap calculations
# start(pAnno)<-end(pAnno)

length(pAnno)
#[1]  421368

##get overlaps
ol<-findOverlaps(pAnno,cgiObj)

ol2<-findOverlaps(pAnno,shores)

length(queryHits(ol))==length(unique(queryHits(ol)))
#[1] TRUE
length(queryHits(ol2))==length(unique(queryHits(ol2)))
#[1] TRUE

myCgi<-names(pAnno)[queryHits(ol)]
myShore<-names(pAnno)[queryHits(ol2)]

length(intersect(myCgi,myShore))
#[1] 0

##do CGIs
length(myCgi)
#[1] 135080

ol<-intersect(annoObj$illuminaID,myCgi)

length(ol)
#[1] 135080

annoObj$isCgi<-0
annoObj[ol,"isCgi"]<-1

##do Shores
length(myShore)
#[1] 100373

ol<-intersect(annoObj$illuminaID,myShore)

length(ol)
#[1] 100373

annoObj$isShore<-0
annoObj[ol,"isShore"]<-1

table(annoObj$isShore,annoObj$isCgi)
  #        0      1
  # 0 185915 135080
  # 1 100373      0

annoObj$isOcean<-1
annoObj$isOcean[annoObj$isShore==1 | annoObj$isCgi==1]<-0

table(annoObj$isOcean,annoObj$isShore)
  #        0      1
  # 0 135080 100373
  # 1 185915      0

table(annoObj$isOcean,annoObj$isCgi)
  #        0      1
  # 0 100373 135080
  # 1 185915      0

annoObj$cgiClass<-""
annoObj$cgiClass[annoObj$isOcean==1]<-"ocean"
annoObj$cgiClass[annoObj$isShore==1]<-"shore"
annoObj$cgiClass[annoObj$isCgi==1]<-"cgi"

chisq.test(table(annoObj$weberClass,annoObj$isShore))$res
  #             0         1
  # HCP -33.86372  60.55853
  # ICP -13.49352  24.13048
  # LCP  55.50627 -99.26193

##Plot CpG O/E distribution
pdf(file=paste0(HOME,"/annotateFeatures/dens_cpgOeSaxonovForCgiStatus.pdf"),width=12,height=12,useDingbats=FALSE)
par(fig=c(0,.5,.5,1),mar=c(5.1,4.1,4.1,2.1),font=2,font.sub=2,font.lab=2,font.axis=2,las=1)

plot( density(annoObj$saxonovOE[annoObj$isCgi==1]),xlim=c(0,1.2),ylim=c(0,8),
      col=brewer.pal(9,"Greens")[7],lwd=3,lty=2,axes=F,
      sub="",ylab="Density",xlab="CpG O/E",main="850K probe CpG density by feature class"
      )
lines( density(annoObj$saxonovOE[annoObj$weberClass=="HCP"]),col=brewer.pal(9,"Greens")[5],lwd=3)
lines( density(annoObj$saxonovOE[annoObj$saxonovClass=="HCG"]),col=brewer.pal(9,"Greens")[3],lwd=3)

lines( density(annoObj$saxonovOE[annoObj$isShore==1 ]),col=brewer.pal(9,"Oranges")[7],lwd=3,lty=2)
lines( density(annoObj$saxonovOE[annoObj$weberClass=="ICP"]),col=brewer.pal(9,"Oranges")[5],lwd=3)

lines( density(annoObj$saxonovOE[annoObj$isOcean==1 ]),col=brewer.pal(9,"Blues")[7],lwd=3,lty=2)
lines( density(annoObj$saxonovOE[annoObj$weberClass=="LCP"]),col=brewer.pal(9,"Blues")[5],lwd=3)
lines( density(annoObj$saxonovOE[annoObj$saxonovClass=="LCG"]),col=brewer.pal(9,"Blues")[3],lwd=3)

axis(1,lwd=2,las=1,font=2)
axis(2,lwd=2,las=2,font=2)
legend("topright",legend=c("CGI","Weber HCP","Saxonov HCG","Shore","Weber ICP","Ocean","Weber LCP","Saxonov LCG"),
      col=c(brewer.pal(9,"Greens")[c(7,5,3)],brewer.pal(9,"Oranges")[c(7,5)],brewer.pal(9,"Blues")[c(7,5,3)]),
      lty=c(2,1,1,2,1,2,1,1),bty="n",cex=1.25,lwd=3,text.font=2)

##add promoter data
par(fig=c(.5,1,.5,1),mar=c(5.1,4.1,4.1,2.1),font=2,font.sub=2,font.lab=2,font.axis=2,las=1,new=T)
plot( 1,xlim=c(0,1.2),ylim=c(0,6),type="n",
      col=brewer.pal(9,"Greens")[7],lwd=3,lty=2,axes=F,
      sub="",ylab="Density",xlab="CpG O/E",main="Promoter CpG density by feature class"
      )
##need to scale to height of peak
lines( density(objectsSaxonov$promOEsaxonov[objectsWeber$promClassWeber=="HCP"]),col=brewer.pal(9,"Greens")[5],lwd=3)
lines( density(objectsSaxonov$promOEsaxonov[objectsSaxonov$promClassSaxonov=="HCG"]),col=brewer.pal(9,"Greens")[3],lwd=3)

lines( density(objectsSaxonov$promOEsaxonov[objectsWeber$promClassWeber=="ICP"]),col=brewer.pal(9,"Oranges")[5],lwd=3)

lines( density(objectsSaxonov$promOEsaxonov[objectsWeber$promClassWeber=="LCP"]),col=brewer.pal(9,"Blues")[5],lwd=3)
lines( density(objectsSaxonov$promOEsaxonov[objectsSaxonov$promClassSaxonov=="LCG"]),col=brewer.pal(9,"Blues")[3],lwd=3)

axis(1,lwd=2,las=1,font=2)
axis(2,lwd=2,las=2,font=2)
legend("topright",legend=c("Weber HCP","Saxonov HCG","Weber ICP","Weber LCP","Saxonov LCG"),
      col=c(brewer.pal(9,"Greens")[c(5,3)],brewer.pal(9,"Oranges")[c(5)],brewer.pal(9,"Blues")[c(5,3)]),
      lty=c(2,1,1,2,1,2,1,1),bty="n",cex=1.25,lwd=3,text.font=2)

dev.off()

rm(ol,ol2,shores,cgiObj,islands_df,myCgi,myShore)

################################################################################
################################################################################
###Save new annotations

##add annotation for prom/body/non-genic
table(annoObj$isDistal,annoObj$hasPromOverlap)
  #        0      1
  # 0 128394 122709
  # 1 170265      0

table(annoObj$isGenic,annoObj$hasPromOverlap)
  #        0      1
  # 0  74819  43039
  # 1 223840  79670

table(annoObj$isGenic,annoObj$isDistalNonGenic)
  #        0      1
  # 0  74688  43170
  # 1 303510      0

table(annoObj$hasPromOverlap==1 , annoObj$hasUp5kbOverlap==1 | annoObj$hasDn5kbOverlap==1)
  #        FALSE   TRUE
  # FALSE 170265 128394
  # TRUE   82416  40293 ##some doubles.. Prom > 5kb

##build hierarchy of annotations
  ## body < 5kb up < 5kb dn < distal < prom
annoObj$featureClass<-""
annoObj$featureClass[annoObj$isGenic==1]<-"distal body"
annoObj$featureClass[annoObj$hasUp5kbOverlap==1 | annoObj$hasDn5kbOverlap==1]<-"proximal up"
annoObj$featureClass[annoObj$hasDn5kbOverlap==1]<-"proximal dn"
annoObj$featureClass[annoObj$isDistalNonGenic==1]<-"distal"
annoObj$featureClass[annoObj$hasPromOverlap==1]<-"promoter"

table(annoObj$featureClass)
     # distal distal body    promoter proximal dn proximal up 
     #  43170      127095      122709       86190       42204 

##also add for knownGene...
  ## body < 5kb up < 5kb dn < distal < prom
annoObj$featureClassUcsc<-""
annoObj$featureClassUcsc[annoObj$ucscKnownGeneIsGenic==1]<-"distal body"
annoObj$featureClassUcsc[annoObj$hasUCSCknownGeneUp5kbOverlap==1 | annoObj$hasUCSCknownGeneDn5kbOverlap==1]<-"proximal up"
annoObj$featureClassUcsc[annoObj$hasUCSCknownGeneDn5kbOverlap==1]<-"proximal dn"
annoObj$featureClassUcsc[annoObj$ucscKnownGeneIsDistalNonGenic==1]<-"distal"
annoObj$featureClassUcsc[annoObj$hasUCSCknownGenePromOverlap==1]<-"promoter"

table(annoObj$featureClassUcsc)
     # distal distal body    promoter proximal dn proximal up 
     #  77007      138020      105072       66025       35244 

table(annoObj$featureClass,annoObj$featureClassUcsc)
  #             distal distal body promoter proximal dn proximal up
  # distal       39931        2521      207         229         282
  # distal body  11546      114077      508         563         401
  # promoter      9137        5659   100972        4347        2594
  # proximal dn  10234        8342     2507       60286        4821
  # proximal up   6159        7421      878         600       27146

table(annoObj$featureClass,annoObj$featureClassUcsc,annoObj$hasAtacOverlap)
# , ,  = 0             
#               distal distal body promoter proximal dn proximal up
#   distal       32380        1969      111         164         200
#   distal body   9617       92262      258         411         310
#   promoter      4760        3391    37329        2134        1534
#   proximal dn   7396        6179     1120       37947        3518
#   proximal up   4750        5836      488         453       20737
# , ,  = 1
#               distal distal body promoter proximal dn proximal up
#   distal        7551         552       96          65          82
#   distal body   1929       21815      250         152          91
#   promoter      4377        2268    63643        2213        1060
#   proximal dn   2838        2163     1387       22339        1303
#   proximal up   1409        1585      390         147        6409

table(annoObj$featureClass,annoObj$cgiClass)
  #               cgi ocean shore
  # distal       6613 30095  6462
  # distal body 21428 85053 20614
  # promoter    71097 22403 29209
  # proximal dn 28036 28188 29966
  # proximal up  7906 20176 14122

##add CpG vs CpH identiifer
table( sub("(..).+","\\1",annoObj$illuminaID) )
#     cg     ch 
# 418858   2510 

annoObj$illuminaCpG_CpH_Probe<-sub("(..).+","\\1",annoObj$illuminaID)

annoObj$illuminaDesignType<-probesKeep$designType

save(annoObj,file=paste0(HOME,"/annotateFeatures/object_noChromFilter_",nrow(annoObj),"x",ncol(annoObj),"_expandedAnnotations.RData"))

################################################################################
################################################################################
##SCREEN regulatory genome annotations

##Search Candidate cis-Regulatory Elements
#https://screen.wenglab.org/

#This track displays the ENCODE Registry of candidate cis-Regulatory Elements (cCREs) in the human genome
# , a total of 926,535 elements identified and classified by the ENCODE Data Analysis Center according to biochemical signatures. 
# cCREs are the subset of representative DNase hypersensitive sites across ENCODE and Roadmap Epigenomics samples that are supported by either histone modifications (H3K4me3 and H3K27ac)
#  or CTCF-binding data. The Registry of cCREs is one of the core components of the integrative level of the ENCODE Encyclopedia of DNA Elements.


#ENCODE Project Consortium, Jill E. Moore, Michael J. Purcaro, Henry E. Pratt, Charles B. Epstein, Noam Shoresh, Jessika Adrian, et al. 2020. “Expanded Encyclopaedias of DNA Elements in the Human and Mouse Genomes.” Nature 583 (7818): 699–710.

#http://genome.ucsc.edu/cgi-bin/hgc?hgsid=1003912663_jB6Nd4sVqe6Xsd5kv0owpLXUdKoU&c=chr14&l=37588433&r=37590435&o=37588582&t=37588931&g=encodeCcreCombined&i=EH38E1709499
# red   prom  promoter-like signature       PLS
# orange  enhP  proximal enhancer-like signature  pELS
# yellow  enhD  distal enhancer-like signature    dELS
# pink    K4m3  DNase-H3K4me3           DNase-H3K4me3
# blue    CTCF  CTCF-only             CTCF-only

GET("https://api.wenglab.org/screen_v13/fdownloads/GRCh38-ccREs.bed",write_disk(paste0(HOME,"/annotateFeatures/GRCh38-ccREs.bed"), overwrite=TRUE))

as.character( md5sum( paste0(HOME,"/annotateFeatures/GRCh38-ccREs.bed") ) )
#[1] "bac32bbf3040f5dd233b3d14605d71b9"

cres<-readLines(paste0(HOME,"/annotateFeatures/GRCh38-ccREs.bed"))

head(cres)
# [1] "chr1\t181251\t181601\tEH38D2115304\tEH38E1310153\tpELS,CTCF-bound"
# [2] "chr1\t190865\t191071\tEH38D2115305\tEH38E1310154\tdELS,CTCF-bound"
# [3] "chr1\t778562\t778912\tEH38D2115324\tEH38E1310158\tPLS,CTCF-bound" 
# [4] "chr1\t779086\t779355\tEH38D2115325\tEH38E1310159\tPLS,CTCF-bound" 
# [5] "chr1\t779727\t780060\tEH38D2115326\tEH38E1310160\tpELS,CTCF-bound"
# [6] "chr1\t817080\t817403\tEH38D2115333\tEH38E1310166\tPLS,CTCF-bound" 

##coordinates + anno. some have comma-sep CTCF-annotation

cres<-lapply(cres,function(x) unlist(strsplit(x,"\t|,")))

table( unlist(lapply(cres,length)) )
#      6      7 
# 537721 388814 

cres<-lapply(cres,function(x) { if(length(x)!=7 ) { return(c(x,""))} else { return(x) } })

table( unlist(lapply(cres,length)) )
#      7 
# 926535 

cres<-do.call("rbind",cres)

colnames(cres)<-c("chr","start","end","id","id2","category","ctcf")
cres<-as.data.frame(cres,stringsAsFactors=FALSE)

cres$start<-as.integer(cres$start)
cres$end<-as.integer(cres$end)

table(cres$category,cres$ctcf)
  #                      CTCF-bound
  # CTCF-only          0      56766
  # dELS          448981     218618
  # DNase-H3K4me3  16737       8800
  # pELS           64421      77409
  # PLS             7582      27221

table(cres$chr)
#  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
# 83575 45532 45984 46113 25949 30498 29482 29435 37385 21934 26322 77347 24913 
# chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX  chrY 
# 11374 16784 62018 47225 53182 55781 48615 44459 38949 22671  1008 

##single coordinate position
cresObj<-GRanges(seqnames=cres$chr,ranges=IRanges(start=cres$start+1,end=cres$end)) ##+1 for de-BED

cresObj$category<-cres$category
cresObj$ctcf<-cres$ctcf

##any overlapping peaks?
sum( queryHits(findOverlaps(cresObj)) != subjectHits(findOverlaps(cresObj)) )
#[1] 0

##single coordinate position +/-500kb
cpgObj<-GRanges(seqnames=annoObj$chr,ranges=IRanges(start=annoObj$end,end=annoObj$end))

##annos
cpgObj$name<-annoObj$illuminaID
cpgObj$coord1<-annoObj$end

##overlap object defining comparisons
ol<-findOverlaps(cresObj,cpgObj)

ol
# Hits object with 207402 hits and 0 metadata columns:
#            queryHits subjectHits
#            <integer>   <integer>
#        [1]        13           3
#        [2]        15          12
#        [3]        16          13
#        [4]        18          25
#        [5]        18          26
#        ...       ...         ...
#   [207398]    926265      420902
#   [207399]    926265      420903
#   [207400]    926265      420904
#   [207401]    926270      421078
#   [207402]    926314      421181
#   -------
#   queryLength: 926535 / subjectLength: 421368

length(unique(queryHits(ol)))
#[1] 118240

length(unique(subjectHits(ol)))
#[1] 207402

ol2<-unlist(split(queryHits(ol),subjectHits(ol)))

head(ol2)
    #  3     12     13     15     22     23 
    # 13     15     16 534919 534922 534922 

table(cresObj$category[ol2])
    # CTCF-only          dELS DNase-H3K4me3          pELS           PLS 
    #      1283         60697          3550         65982         75890 

all(annoObj$illuminaID[as.integer(names(ol2))]==cpgObj$name[as.integer(names(ol2))])
#[1] TRUE

table(cresObj$ctcf[ol2])
     #       CTCF-bound 
     # 53499     153903 

annoObj$encodeCre<-0
annoObj[annoObj$illuminaID[as.integer(names(ol2))],"encodeCre"]<-1

annoObj$encodeCreCategory<-""
annoObj[annoObj$illuminaID[as.integer(names(ol2))],"encodeCreCategory"]<-cresObj$category[ol2]

annoObj$encodeCreCtcf<-0
annoObj[annoObj$illuminaID[as.integer(names(ol2))],"encodeCreCtcf"]<-as.integer(cresObj$ctcf[ol2]!="")

table(annoObj$encodeCreCategory,annoObj$hasAtacOverlap)
  #                    0      1
  #               179781  34185
  # CTCF-only       1083    200
  # dELS           39762  20935
  # DNase-H3K4me3   2940    610
  # pELS           34604  31378
  # PLS            17084  58806

table(annoObj$encodeCreCategory,annoObj$weberClass)
  #                 HCP   ICP   LCP
  #               51772 81854 80340
  # CTCF-only        58   593   632
  # dELS          12773 24109 23815
  # DNase-H3K4me3  1656  1568   326
  # pELS          46941 12206  6835
  # PLS           60701 10042  5147

table(annoObj$encodeCreCategory,annoObj$cgiClass)
  #                  cgi  ocean  shore
  #                39848 123441  50677
  # CTCF-only         63   1065    155
  # dELS           12198  39057   9442
  # DNase-H3K4me3   1455   1194    901
  # pELS           29495  10075  26412
  # PLS            52021  11083  12786

table(annoObj$encodeCreCategory,annoObj$hasPromOverlap)
  #                    0      1
  #               176912  37054
  # CTCF-only       1158    125
  # dELS           60121    576
  # DNase-H3K4me3   3291    259
  # pELS           42165  23817
  # PLS            15012  60878

table(annoObj$encodeCreCategory,annoObj$featureClass)
  #               distal distal body promoter proximal dn proximal up
  #                26257       82305    37054       42728       25622
  # CTCF-only        312         528      125         149         169
  # dELS           15411       31074      576        8492        5144
  # DNase-H3K4me3   1051         958      259         696         586
  # pELS              63        5904    23817       26283        9915
  # PLS               76        6326    60878        7842         768

chisq.test(table(annoObj$encodeCreCategory,annoObj$featureClass))$res               
  #                    distal distal body    promoter proximal dn proximal up
  #                 29.284098   69.939438 -101.178770   -4.963225   28.630432
  # CTCF-only       15.748301    7.168364  -12.862709   -7.002234    3.572286
  # dELS           116.570421   94.351209 -128.618399  -35.211711  -11.996607
  # DNase-H3K4me3   36.038689   -3.446174  -24.097774   -1.118694   12.220418
  # pELS           -81.452906  -99.223293   33.199125  110.063078   40.670591
  # PLS            -87.314524 -109.483162  260.844464  -61.650466  -78.375380


pdf(paste0(HOME,"/annotateFeatures/barplot_encodeCreOverlaps_noChrFilter.pdf"),width=12,height=12,useDingbats=FALSE)
par(fig=c(0,.5,.5,1),mar=c(6.1,3.1,3.1,1.1),font=2,font.main=2,font.sub=2,font.lab=2,font.axis=2,las=1,new=F)

ss<-table(factor(annoObj$encodeCreCategory,levels=c("PLS","pELS","dELS","DNase-H3K4me3","CTCF-only","")),paste0(sub(" body| up| dn","",annoObj$featureClass),c(""," ATAC")[1+annoObj$hasAtacOverlap]))
barplot(t( t(ss)/colSums(ss) ),col=c("firebrick2","orange","yellow","pink","dodgerblue","lightgrey"),las=2,
  ,axes=FALSE,sub="",xlab="",ylab="",main="CpG - ENCODE CRE overlaps",
  xlim=c(0.5,15),ylim=c(0,1),cex.main=1.5
  )
legend("topright",legend=c("PLS","pELS","dELS","DNase-H3K4me3","CTCF-only","none"),
  col=c("firebrick2","orange","yellow","pink","dodgerblue","lightgrey"),bty="n",
  cex=1.5,text.font=2,pch=15)
axis(2,cex.axis=1,font=2,las=1,at=seq(0,1,by=.25),lwd=2)

par(fig=c(.5,1,.5,1),mar=c(6.1,3.1,3.1,1.1),font=2,font.main=2,font.sub=2,font.lab=2,font.axis=2,las=1,new=T)

barplot(t( t(ss) ),col=c("firebrick2","orange","yellow","pink","dodgerblue","lightgrey"),las=2,
  ,axes=FALSE,sub="",xlab="",ylab="",main="CpG - ENCODE CRE overlaps",
  xlim=c(0.5,15),ylim=c(0,400e3),cex.main=1.5
  )
legend("topright",legend=c("PLS","pELS","dELS","DNase-H3K4me3","CTCF-only","none"),
  col=c("firebrick2","orange","yellow","pink","dodgerblue","lightgrey"),bty="n",
  cex=1.5,text.font=2,pch=15)
axis(2,cex.axis=1,font=2,las=1,at=seq(0,500e3,by=200e3),labels=c("0","200K","400K"),lwd=2)

par(fig=c(0,.5,0,.5),mar=c(6.1,3.1,3.1,1.1),font=2,font.main=2,font.sub=2,font.lab=2,font.axis=2,las=1,new=T)
ss<-table(factor(annoObj$encodeCreCategory,levels=c("PLS","pELS","dELS","DNase-H3K4me3","CTCF-only","")),paste0(annoObj$cgiClass,c(""," ATAC")[1+annoObj$hasAtacOverlap]))
barplot(t( t(ss)/colSums(ss) ),col=c("firebrick2","orange","yellow","pink","dodgerblue","lightgrey"),las=2,
  ,axes=FALSE,sub="",xlab="",ylab="",main="CpG - ENCODE CRE overlaps",
  xlim=c(0.5,15),ylim=c(0,1),cex.main=1.5
  )
legend("topright",legend=c("PLS","pELS","dELS","DNase-H3K4me3","CTCF-only","none"),
  col=c("firebrick2","orange","yellow","pink","dodgerblue","lightgrey"),bty="n",
  cex=1.5,text.font=2,pch=15)
axis(2,cex.axis=1,font=2,las=1,at=seq(0,1,by=.25),lwd=2)

par(fig=c(.5,1,0,.5),mar=c(6.1,3.1,3.1,1.1),font=2,font.main=2,font.sub=2,font.lab=2,font.axis=2,las=1,new=T)
plot(1,type="n",sub="",main="",xlab="",ylab="",axes=FALSE,xlim=c(0,1),ylim=c(0,1))
legend("center",legend=c("Promoter-like sequence","Proximal enhancer-like sequnece","Distal enhancer-like sequnece","DNase-H3K4me3","CTCF-only","none"),
  col=c("firebrick2","orange","yellow","pink","dodgerblue","lightgrey"),bty="n",
  cex=1.5,text.font=2,pch=15)

dev.off()

rm(ss)

save(annoObj,file=paste0(HOME,"/annotateFeatures/object_noChromFilter_",nrow(annoObj),"x",ncol(annoObj),"_expandedAnnotations.RData"))

################################################################################
################################################################################
##check ENCODE TFB across feature sets

##argument for specificity of signal

download.file("https://hgdownload.cse.ucsc.edu/goldenpath/hg38/encRegTfbsClustered/encRegTfbsClusteredWithCells.hg38.bed.gz",
  paste0(HOME,"/annotateFeatures/encRegTfbsClusteredWithCells.hg38.bed.gz"))

as.character(md5sum(paste0(HOME,"/annotateFeatures/encRegTfbsClusteredWithCells.hg38.bed.gz")))=="bc7220dd56257b107d7d56f22bc97cc0"
# [1]  TRUE 

# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/
#         * encRegTfbsClusteredInputs.txt.gz
#         * encRegTfbsClusteredSources.txt.gz

tfbs<-readLines(paste0(HOME,"/annotateFeatures/encRegTfbsClusteredWithCells.hg38.bed.gz"))

length(tfbs)
#[1] 10565630

head(tfbs)
# [1] "chr1\t921079\t921415\tTCF7\t429\tK562"    "chr1\t984163\t984499\tTCF7\t473\tK562"    "chr1\t995778\t996134\tTCF7\t239\tHepG2"  
# [4] "chr1\t1000694\t1001508\tTCF7\t354\tK562"  "chr1\t1013171\t1013708\tTCF7\t431\tK562"  "chr1\t1124686\t1125042\tTCF7\t284\tHepG2"

##coordinates + anno. some have comma-sep CTCF-annotation

tfbs<-lapply(tfbs,function(x) unlist(strsplit(x,"\t")))

table( unlist(lapply(tfbs,length)) )
#        6 
# 10565630 

tfbs.l<-unlist(lapply(tfbs,function(x) { length(unlist(strsplit(x[6],","))) }))

table( unlist(lapply(tfbs,length)) )
#      7 
# 926535 

table(tfbs.l)
# tfbs.l
#       1       2       3       4       5       6       7       8       9      10      11      12      13      14      15      16      17      18 
# 8623023 1141700  346297  162093   67846   44024   28644   32699   10120    9912    4403    5054    2865    2554    2272    2065    1962    1825 
#      19      20      21      22      23      24      25      26      27      28      29      30      31      32      33      34      35      36 
#    1659    1629    1533    1490    1369    1414    1301    1249    1278    1319    1284    1320    1322    1350    1443    1580    1845    2027 
#      37      38      39      40      41      42      43      44      45      46      47      48      49      50      51      52      53      54 
#    2425     600     579     599     539     560     556     494     554     524     493     467     424     462     466     427     415     433 
#      55      56      57      58      59      60      61      62      63      64      65      66      67      68      69      70      71      72 
#     416     420     434     411     390     396     389     380     390     341     380     367     371     373     374     364     319     401 
#      73      74      75      76      77      78      79      80      81      82      83      84      85      86      87      88      89      90 
#     352     368     384     406     370     384     355     368     397     388     364     408     397     382     383     390     410     404 
#      91      92      93      94      95      96      97      98      99     100     101     102     103     104     105     106     107     108 
#     430     442     387     472     477     471     477     533     566     592     597     621     660     692     824     934    1065    1274 
#     109     110     111     112     113     114 
#    1452    1807    2332    3560    4348       4 

tfbs<-do.call("rbind",tfbs)

head(tfbs)
#      [,1]   [,2]      [,3]      [,4]   [,5]  [,6]   
# [1,] "chr1" "921079"  "921415"  "TCF7" "429" "K562" 
# [2,] "chr1" "984163"  "984499"  "TCF7" "473" "K562" 
# [3,] "chr1" "995778"  "996134"  "TCF7" "239" "HepG2"
# [4,] "chr1" "1000694" "1001508" "TCF7" "354" "K562" 
# [5,] "chr1" "1013171" "1013708" "TCF7" "431" "K562" 
# [6,] "chr1" "1124686" "1125042" "TCF7" "284" "HepG2"

colnames(tfbs)<-c("chr","start","end","tf","score","cells")
tfbs<-as.data.frame(tfbs,stringsAsFactors=FALSE)

tfbs$start<-as.integer(tfbs$start)
tfbs$end<-as.integer(tfbs$end)
tfbs$score<-as.integer(tfbs$score)
tfbs$nCells<-tfbs.l

rm(tfbs.l)

str(tfbs)
# 'data.frame':   10565630 obs. of  7 variables:
#  $ chr   : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start : int  921079 984163 995778 1000694 1013171 1124686 1125284 1164818 1165845 1167766 ...
#  $ end   : int  921415 984499 996134 1001508 1013708 1125042 1125557 1165174 1166201 1168122 ...
#  $ tf    : chr  "TCF7" "TCF7" "TCF7" "TCF7" ...
#  $ score : int  429 473 239 354 431 284 1000 195 303 195 ...
#  $ cells : chr  "K562" "K562" "HepG2" "K562" ...
#  $ nCells: int  1 1 1 1 1 1 1 1 1 1 ...

tfbs<-tfbs[tfbs$chr %in% paste0("chr",c(1:22,"X")),]

table(tfbs$chr)
#    chr1   chr10   chr11   chr12   chr13   chr14   chr15   chr16   chr17   chr18   chr19    chr2   chr20   chr21   chr22    chr3    chr4    chr5 
# 1085445  472779  510155  512796  205293  317421  330182  415073  541692  188198  466013  826472  369444  128715  205542  640167  416762  533383 
#    chr6    chr7    chr8    chr9    chrX 
#  716419  580685  458675  403324  214822 

names(table(tfbs$tf))
#   [1] "AFF1"     "AGO1"     "AGO2"     "ARHGAP35" "ARID1B"   "ARID2"    "ARID3A"   "ARNT"     "ASH1L"    "ASH2L"    "ATF2"     "ATF3"     "ATF4"    
#  [14] "ATF7"     "ATM"      "BACH1"    "BATF"     "BCL11A"   "BCL3"     "BCOR"     "BHLHE40"  "BMI1"     "BRCA1"    "BRD4"     "BRD9"     "C11orf30"
#  [27] "CBFA2T2"  "CBFA2T3"  "CBFB"     "CBX1"     "CBX2"     "CBX3"     "CBX5"     "CBX8"     "CC2D1A"   "CCAR2"    "CDC5L"    "CEBPB"    "CHAMP1"  
#  [40] "CHD1"     "CHD4"     "CHD7"     "CLOCK"    "COPS2"    "CREB1"    "CREB3L1"  "CREBBP"   "CREM"     "CTBP1"    "CTCF"     "CUX1"     "DACH1"   
#  [53] "DEAF1"    "DNMT1"    "DPF2"     "E2F1"     "E2F4"     "E2F6"     "E2F7"     "E2F8"     "E4F1"     "EBF1"     "EED"      "EGR1"     "EHMT2"   
#  [66] "ELF1"     "ELF4"     "ELK1"     "EP300"    "EP400"    "ESR1"     "ESRRA"    "ETS1"     "ETV4"     "ETV6"     "EWSR1"    "EZH2"     "FIP1L1"  
#  [79] "FOS"      "FOSL1"    "FOSL2"    "FOXA1"    "FOXA2"    "FOXK2"    "FOXM1"    "FOXP1"    "FUS"      "GABPA"    "GABPB1"   "GATA1"    "GATA2"   
#  [92] "GATA3"    "GATA4"    "GATAD2A"  "GATAD2B"  "GMEB1"    "HCFC1"    "HDAC1"    "HDAC2"    "HDAC3"    "HDAC6"    "HES1"     "HMBOX1"   "HNF1A"   
# [105] "HNF4A"    "HNF4G"    "HNRNPH1"  "HNRNPK"   "HNRNPL"   "HNRNPLL"  "HNRNPUL1" "HSF1"     "IKZF1"    "IKZF2"    "IRF1"     "IRF2"     "IRF3"    
# [118] "IRF4"     "IRF5"     "JUN"      "JUNB"     "JUND"     "KAT2A"    "KAT2B"    "KAT8"     "KDM1A"    "KDM4A"    "KDM4B"    "KDM5A"    "KDM5B"   
# [131] "KLF16"    "KLF5"     "L3MBTL2"  "LCORL"    "LEF1"     "MAFF"     "MAFK"     "MAX"      "MBD2"     "MCM2"     "MCM3"     "MCM5"     "MCM7"    
# [144] "MEF2A"    "MEF2B"    "MEF2C"    "MEIS2"    "MGA"      "MIER1"    "MITF"     "MLLT1"    "MNT"      "MTA1"     "MTA2"     "MTA3"     "MXI1"    
# [157] "MYB"      "MYBL2"    "MYC"      "MYNN"     "NANOG"    "NBN"      "NCOA1"    "NCOA2"    "NCOA3"    "NCOA4"    "NCOA6"    "NCOR1"    "NEUROD1" 
# [170] "NFATC1"   "NFATC3"   "NFE2"     "NFE2L2"   "NFIB"     "NFIC"     "NFRKB"    "NFXL1"    "NFYA"     "NFYB"     "NR0B1"    "NR2C1"    "NR2C2"   
# [183] "NR2F1"    "NR2F2"    "NR2F6"    "NR3C1"    "NRF1"     "NUFIP1"   "PAX5"     "PAX8"     "PBX3"     "PCBP1"    "PCBP2"    "PHB2"     "PHF20"   
# [196] "PHF21A"   "PHF8"     "PKNOX1"   "PLRG1"    "PML"      "POLR2A"   "POLR2G"   "POU2F2"   "PRDM10"   "PRPF4"    "PTBP1"    "PYGO2"    "RAD21"   
# [209] "RAD51"    "RB1"      "RBBP5"    "RBFOX2"   "RBM14"    "RBM15"    "RBM17"    "RBM22"    "RBM25"    "RBM34"    "RBM39"    "RCOR1"    "RELB"    
# [222] "REST"     "RFX1"     "RFX5"     "RLF"      "RNF2"     "RUNX1"    "RUNX3"    "RXRA"     "SAFB"     "SAFB2"    "SAP30"    "SETDB1"   "SIN3A"   
# [235] "SIN3B"    "SIRT6"    "SIX4"     "SIX5"     "SKI"      "SKIL"     "SMAD1"    "SMAD2"    "SMAD5"    "SMARCA4"  "SMARCA5"  "SMARCB1"  "SMARCC2" 
# [248] "SMARCE1"  "SMC3"     "SNRNP70"  "SOX13"    "SOX6"     "SP1"      "SPI1"     "SREBF1"   "SREBF2"   "SRF"      "SRSF4"    "SRSF7"    "SRSF9"   
# [261] "STAT1"    "STAT2"    "STAT3"    "STAT5A"   "SUPT20H"  "SUZ12"    "TAF1"     "TAF15"    "TAF7"     "TAF9B"    "TAL1"     "TBL1XR1"  "TBP"     
# [274] "TBX21"    "TBX3"     "TCF12"    "TCF7"     "TCF7L2"   "TEAD4"    "TFAP4"    "THAP1"    "THRA"     "TRIM22"   "TRIM24"   "TRIM28"   "TRIP13"  
# [287] "U2AF1"    "U2AF2"    "UBTF"     "USF1"     "USF2"     "WHSC1"    "WRNIP1"   "XRCC3"    "XRCC5"    "YY1"      "ZBED1"    "ZBTB1"    "ZBTB11"  
# [300] "ZBTB2"    "ZBTB33"   "ZBTB40"   "ZBTB5"    "ZBTB7A"   "ZBTB7B"   "ZBTB8A"   "ZEB1"     "ZEB2"     "ZFP91"    "ZFX"      "ZHX1"     "ZHX2"    
# [313] "ZKSCAN1"  "ZMIZ1"    "ZMYM3"    "ZNF143"   "ZNF184"   "ZNF207"   "ZNF217"   "ZNF24"    "ZNF263"   "ZNF274"   "ZNF280A"  "ZNF282"   "ZNF316"  
# [326] "ZNF318"   "ZNF384"   "ZNF407"   "ZNF444"   "ZNF507"   "ZNF512B"  "ZNF574"   "ZNF579"   "ZNF592"   "ZNF639"   "ZNF687"   "ZNF8"     "ZNF830"  
# [339] "ZSCAN29"  "ZZZ3"    

##make object
tfbsObj<-GRanges(seqnames=tfbs$chr,ranges=IRanges(start=tfbs$start+1,end=tfbs$end)) #+1 for de-BED

tfbsObj$nCells<-tfbs$nCells
tfbsObj$tf<-tfbs$tf

##single coordinate position - have used 3' before
cpgObj<-GRanges(seqnames=annoObj$chr,ranges=IRanges(start=annoObj$end,end=annoObj$end))

##some end coordinates also but do not need to correct as these cannot find match..
cpgObj$name<-annoObj$illuminaID
cpgObj$coord1<-annoObj$end

##overlap object defining comparisons
ol<-findOverlaps(tfbsObj,cpgObj)

ol
# Hits object with 10787155 hits and 0 metadata columns:
#              queryHits subjectHits
#              <integer>   <integer>
#          [1]         1          43
#          [2]         1          44
#          [3]         2         169
#          [4]         3         173
#          [5]         4         187
#          ...       ...         ...
#   [10787151]  10539229      397695
#   [10787152]  10539229      397696
#   [10787153]  10539229      397697
#   [10787154]  10539248      398082
#   [10787155]  10539314      399205
#   -------
#   queryLength: 10539457 / subjectLength: 421368

length(unique(subjectHits(ol)))
#[1] 330644

ol2<-unlist(lapply(split(tfbsObj$nCells[queryHits(ol)],cpgObj$name[subjectHits(ol)]),mean))

head(ol2)
# cg00000029 cg00000108 cg00000165 cg00000289 cg00000292 cg00000321 
#   2.068966   1.125000   4.000000   3.000000   1.250000   3.500000 

annoObj$encodeChipMeanNcellPeaks<-0
annoObj[names(ol2),"encodeChipMeanNcellPeaks"]<-round(ol2,3)

quantile(annoObj$encodeChipMeanNcellPeaks)
  #    0%     25%     50%     75%    100% 
  # 0.000   1.000   1.556   2.774 113.000 

##number of unique TFs
ol2<-unlist(lapply(split(tfbsObj$tf[queryHits(ol)],cpgObj$name[subjectHits(ol)]),function(x) length(unique(x)) ))

head(ol2)
# cg00000029 cg00000108 cg00000165 cg00000289 cg00000292 cg00000321 
#         29          8          2          1          4          2 

quantile(ol2)
  # 0%  25%  50%  75% 100% 
  #  1    3   12   47  258 

annoObj$encodeChipUniqueTfPeaks<-0
annoObj[names(ol2),"encodeChipUniqueTfPeaks"]<-ol2

quantile(annoObj$encodeChipUniqueTfPeaks)
  # 0%  25%  50%  75% 100% 
  #  0    1    6   33  258 

pdf(paste0(HOME,"/annotateFeatures/barplot_encodeTfbsOverlaps_noChrFilter.pdf"),width=12,height=12,useDingbats=FALSE)

##nCells - context
a<-unlist(lapply(split(annoObj$encodeChipMeanNcellPeaks,paste0(sub(" body| up| dn","",annoObj$featureClass),c(""," ATAC")[1+annoObj$hasAtacOverlap])),median))
b<-unlist(lapply(split(annoObj$encodeChipMeanNcellPeaks,paste0(sub(" body| up| dn","",annoObj$featureClass),c(""," ATAC")[1+annoObj$hasAtacOverlap])),quantile,.25))
c<-unlist(lapply(split(annoObj$encodeChipMeanNcellPeaks,paste0(sub(" body| up| dn","",annoObj$featureClass),c(""," ATAC")[1+annoObj$hasAtacOverlap])),quantile,.75))

par(fig=c(0,.5,.5,1),mar=c(8.5,4.1,2.1,0.1),cex.main=1.5,font=2,font.main=2,font.sub=2,font.lab=2,font.axis=2,las=1,new=F)
plot(1,type="n",main = "CpG - ChIP-seq peak overlaps",sub="" , xlab = "", ylab = "",
  xlim=c(0.5,7),ylim=c(0,4),axes=F
  )   
cols<-c(brewer.pal(9,"Oranges")[5],brewer.pal(9,"Oranges")[5],brewer.pal(9,"Blues")[7],brewer.pal(9,"Blues")[7],brewer.pal(9,"Blues")[5],brewer.pal(9,"Blues")[5])
for(i in 1:6) {
  rect(xleft=i-.4,xright=i+.4,ytop=c[i],ybottom=b[i],lwd=3,col=cols[i]) 
  lines(x=rep(i,2)-c(.4,-.4),,y=rep(a[i],2),lwd=3,col=1)
}
mtext(text="Number of cell lines (median)",side=2,las=3,font=2,line=2.5,cex=1.5)
axis(2,cex.axis=1.5,font=2,las=1,at=seq(0,4,by=1),lwd=3)
axis(1,cex.axis=1.5,font=2,las=2,at=seq(1,6,by=1),lwd=0,labels=c("Distal","Distal ATAC","Promoter","Prom ATAC","Proximal","Prox ATAC"))

##nTFs - context
a<-unlist(lapply(split(annoObj$encodeChipUniqueTfPeaks,paste0(sub(" body| up| dn","",annoObj$featureClass),c(""," ATAC")[1+annoObj$hasAtacOverlap])),median))
b<-unlist(lapply(split(annoObj$encodeChipUniqueTfPeaks,paste0(sub(" body| up| dn","",annoObj$featureClass),c(""," ATAC")[1+annoObj$hasAtacOverlap])),quantile,.25))
c<-unlist(lapply(split(annoObj$encodeChipUniqueTfPeaks,paste0(sub(" body| up| dn","",annoObj$featureClass),c(""," ATAC")[1+annoObj$hasAtacOverlap])),quantile,.75))

par(fig=c(.5,1,.5,1),mar=c(8.5,4.1,2.1,0.1),cex.main=1.5,font=2,font.main=2,font.sub=2,font.lab=2,font.axis=2,las=1,new=T)
plot(1,type="n",main = "CpG - ChIP-seq peak overlaps",sub="" , xlab = "", ylab = "",
  xlim=c(0.5,7),ylim=c(0,120),axes=F
  )   
cols<-c(brewer.pal(9,"Oranges")[5],brewer.pal(9,"Oranges")[5],brewer.pal(9,"Blues")[7],brewer.pal(9,"Blues")[7],brewer.pal(9,"Blues")[5],brewer.pal(9,"Blues")[5])
for(i in 1:6) {
  rect(xleft=i-.4,xright=i+.4,ytop=c[i],ybottom=b[i],lwd=3,col=cols[i]) 
  lines(x=rep(i,2)-c(.4,-.4),,y=rep(a[i],2),lwd=3,col=1)
}
mtext(text="Number of TFs (median)",side=2,las=3,font=2,line=2.5,cex=1.5)
axis(2,cex.axis=1.5,font=2,las=1,at=seq(0,120,by=30),lwd=3)
axis(1,cex.axis=1.5,font=2,las=2,at=seq(1,6,by=1),lwd=0,labels=c("Distal","Distal ATAC","Promoter","Prom ATAC","Proximal","Prox ATAC"))

##nCells - cgi
a<-unlist(lapply(split(annoObj$encodeChipMeanNcellPeaks,paste0(annoObj$cgiClass,c(""," ATAC")[1+annoObj$hasAtacOverlap])),median))
b<-unlist(lapply(split(annoObj$encodeChipMeanNcellPeaks,paste0(annoObj$cgiClass,c(""," ATAC")[1+annoObj$hasAtacOverlap])),quantile,.25))
c<-unlist(lapply(split(annoObj$encodeChipMeanNcellPeaks,paste0(annoObj$cgiClass,c(""," ATAC")[1+annoObj$hasAtacOverlap])),quantile,.75))

par(fig=c(0,.5,0,.5),mar=c(8.5,4.1,2.1,0.1),cex.main=1.5,font=2,font.main=2,font.sub=2,font.lab=2,font.axis=2,las=1,new=T)
plot(1,type="n",main = "CpG - ChIP-seq peak overlaps",sub="" , xlab = "", ylab = "",
  xlim=c(0.5,7),ylim=c(0,4),axes=F
  )   
cols<-c(brewer.pal(9,"Greens")[7],brewer.pal(9,"Greens")[7],brewer.pal(9,"Oranges")[7],brewer.pal(9,"Oranges")[7],brewer.pal(9,"Blues")[7],brewer.pal(9,"Blues")[7])
for(i in 1:6) {
  rect(xleft=i-.4,xright=i+.4,ytop=c[i],ybottom=b[i],lwd=3,col=cols[i]) 
  lines(x=rep(i,2)-c(.4,-.4),,y=rep(a[i],2),lwd=3,col=1)
}
mtext(text="Number of cell lines (median)",side=2,las=3,font=2,line=2.5,cex=1.5)
axis(2,cex.axis=1.5,font=2,las=1,at=seq(0,4,by=1),lwd=3)
axis(1,cex.axis=1.5,font=2,las=2,at=seq(1,6,by=1),lwd=0,labels=c("CGI","CGI ATAC","Ocean","Ocean ATAC","Shore","Shore ATAC"))

##nTFs - cgi
a<-unlist(lapply(split(annoObj$encodeChipUniqueTfPeaks,paste0(annoObj$cgiClass,c(""," ATAC")[1+annoObj$hasAtacOverlap])),median))
b<-unlist(lapply(split(annoObj$encodeChipUniqueTfPeaks,paste0(annoObj$cgiClass,c(""," ATAC")[1+annoObj$hasAtacOverlap])),quantile,.25))
c<-unlist(lapply(split(annoObj$encodeChipUniqueTfPeaks,paste0(annoObj$cgiClass,c(""," ATAC")[1+annoObj$hasAtacOverlap])),quantile,.75))

par(fig=c(.5,1,0,.5),mar=c(8.5,4.1,2.1,0.1),cex.main=1.5,font=2,font.main=2,font.sub=2,font.lab=2,font.axis=2,las=1,new=T)
plot(1,type="n",main = "CpG - ChIP-seq peak overlaps",sub="" , xlab = "", ylab = "",
  xlim=c(0.5,7),ylim=c(0,120),axes=F
  )   
cols<-c(brewer.pal(9,"Greens")[7],brewer.pal(9,"Greens")[7],brewer.pal(9,"Oranges")[7],brewer.pal(9,"Oranges")[7],brewer.pal(9,"Blues")[7],brewer.pal(9,"Blues")[7])
for(i in 1:6) {
  rect(xleft=i-.4,xright=i+.4,ytop=c[i],ybottom=b[i],lwd=3,col=cols[i]) 
  lines(x=rep(i,2)-c(.4,-.4),,y=rep(a[i],2),lwd=3,col=1)
}
mtext(text="Number of TFs (median)",side=2,las=3,font=2,line=2.5,cex=1.5)
axis(2,cex.axis=1.5,font=2,las=1,at=seq(0,120,by=30),lwd=3)
axis(1,cex.axis=1.5,font=2,las=2,at=seq(1,6,by=1),lwd=0,labels=c("CGI","CGI ATAC","Ocean","Ocean ATAC","Shore","Shore ATAC"))

dev.off()

rm(i,cols,a,b,c)

##save results
save(annoObj,file=paste0(HOME,"/annotateFeatures/object_noChromFilter_",nrow(annoObj),"x",ncol(annoObj),"_expandedAnnotations.RData"))

################################################################################
################################################################################
##Get names of Overlapping TFs for each CpG

##make object
tfbsObj<-GRanges(seqnames=tfbs$chr,ranges=IRanges(start=tfbs$start+1,end=tfbs$end)) #+1 for de-BED

tfbsObj$nCells<-tfbs$nCells
tfbsObj$tf<-tfbs$tf

##single coordinate position
cpgObj<-GRanges(seqnames=annoObj$chr,ranges=IRanges(start=annoObj$end,end=annoObj$end))

##some end coordinates also but do not need to correct as these cannot find match..
cpgObj$name<-annoObj$illuminaID
cpgObj$coord1<-annoObj$end

##overlap object defining comparisons
ol<-findOverlaps(tfbsObj,cpgObj)

ol
# Hits object with 10787155 hits and 0 metadata columns:
#              queryHits subjectHits
#              <integer>   <integer>
#          [1]         1          43
#          [2]         1          44
#          [3]         2         169
#          [4]         3         173
#          [5]         4         187
#          ...       ...         ...
#   [10787151]  10539229      397695
#   [10787152]  10539229      397696
#   [10787153]  10539229      397697
#   [10787154]  10539248      398082
#   [10787155]  10539314      399205
#   -------
#   queryLength: 10539457 / subjectLength: 421368

length(unique(subjectHits(ol)))
#[1] 330644

ol2<-unlist(lapply(split(tfbsObj$tf[queryHits(ol)],cpgObj$name[subjectHits(ol)]),function(x) paste(unique(x),collapse="|") ))

head(ol2)
#                                                                                                                                                                       cg00000029 
# "EP300|CTCF|NRF1|CLOCK|KDM5A|GABPA|ZNF143|ZBTB33|ATF2|AGO1|CREB3L1|RELB|HCFC1|ZNF217|ZFP91|NFIC|IKZF1|ZNF592|ZBTB8A|SMAD1|RBFOX2|RXRA|ETV4|ZNF407|TBX21|KLF5|NR2C1|HNRNPLL|CBFB" 
#                                                                                                                                                                       cg00000108 
#                                                                                                                                     "FOXA2|ZBTB33|HNF4A|YY1|SP1|NR2F2|ATF3|RXRA" 
#                                                                                                                                                                       cg00000165 
#                                                                                                                                                                      "EZH2|RNF2" 
#                                                                                                                                                                       cg00000289 
#                                                                                                                                                                          "CEBPB" 
#                                                                                                                                                                       cg00000292 
#                                                                                                                                                      "POLR2A|RBM39|ZBTB7A|KDM5B" 
#                                                                                                                                                                       cg00000321 
#                                                                                                                                                                     "EZH2|HNF4A" 

tfMat<-matrix(0,nrow=nrow(annoObj),ncol=length(names(table(tfbs$tf))))
colnames(tfMat)<-names(table(tfbs$tf))
rownames(tfMat)<-rownames(annoObj)

for( i in unique(tfbsObj$tf)) {
  cat(i,"\n")
  sel<-match(names(ol2[grep(paste0(paste0("^",i,"$"),"|",
  paste0("^",i,"\\|"),"|",
  paste0("\\|",i,"\\|"),"|",
  paste0("\\|",i,"$")),ol2)]),rownames(tfMat))
  if( sum(tfbsObj$tf[queryHits(ol)]==i) != length(sel) ) { stop("something wrong!") }
  tfMat[sel,i]<-1
}

quantile(rowSums(tfMat))
  # 0%  25%  50%  75% 100% 
  #  0    1    6   33  258 

quantile(colSums(tfMat))
    #    0%       25%       50%       75%      100% 
    # 12.00   9491.50  24419.50  46051.25 192718.00 

tfMat<-as.data.frame(tfMat,stringsAsFactors=FALSE)

##save results
save(tfMat,file=paste0(HOME,"/annotateFeatures/object_noChromFilter_",nrow(tfMat),"x",ncol(tfMat),"_encodeTfbsAnnotations.RData"))

################################################################################
################################################################################
###Add repmasker elements

##hg38 - should be 5,607,738 rows
tmp<-tempfile()
download.file("https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz",
  paste0(HOME,"/annotateFeatures/rmsk.hg38.bed.gz")
  )

as.character( md5sum( paste0(HOME,"/annotateFeatures/rmsk.hg38.bed.gz") ) )
#[1] "40874e5f37923660429bfa3fe599c5e0"

rmsk<-readLines( paste0(HOME,"/annotateFeatures/rmsk.hg38.bed.gz") )

length(rmsk)
#[1] 5607738

head(rmsk)

tmp<-c("bin","swScore","milliDiv","milliDel","milliIns","genoName","genoStart","genoEnd","genoLeft","strand","repName","repClass","repFamily","repStart","repEnd","repLeft","id")
 
rmsk<-lapply(rmsk,function(x) { unlist(strsplit(x,"\t")) }) 

unique(unlist(lapply(rmsk,length)))
#[1] 17

rmsk<-do.call("rbind",rmsk)
colnames(rmsk)<-tmp
rm(tmp)

str(rmsk)
 
rmsk38<-rmsk
rm(rmsk)

##Use single-coordinate probe locations - have used 3'
grProbe<-GRanges(seqnames=annoObj$chr,ranges=IRanges(start=annoObj$end,end=annoObj$end)) ##bed
grProbe$illuminaID<-annoObj$illuminaID

grRep<-GRanges(seqnames=rmsk38[,"genoName"],ranges=IRanges(start=as.integer(rmsk38[,"genoStart"])+1,end=as.integer(rmsk38[,"genoEnd"]) ) ) ##is bed

sum(width(grRep)) - sum( as.integer(rmsk38[,"genoEnd"]) - as.integer(rmsk38[,"genoStart"]) )
#[1] 0

ol<-findOverlaps(grProbe,grRep)

ol
# Hits object with 57216 hits and 0 metadata columns:
#           queryHits subjectHits
#           <integer>   <integer>
#       [1]         9        1464
#       [2]        10        1464
#       [3]       108        1508
#       [4]       109        1511
#       [5]       110        1511
#       ...       ...         ...
#   [57212]    421349     5510633
#   [57213]    421350     5510633
#   [57214]    421351     5510634
#   [57215]    421358     5511543
#   [57216]    421362     5511552
#   -------
#   queryLength: 421368 / subjectLength: 5607738

tmp<-split(subjectHits(ol),grProbe$illuminaID[queryHits(ol)])
table(unlist(lapply(tmp,length)))
#     1     2 
# 57028    94 

##some multihits..
is.multi<-names(unlist(lapply(tmp,length)))[unlist(lapply(tmp,length))>1]
is.uni<-names(unlist(lapply(tmp,length)))[unlist(lapply(tmp,length))==1]

length(is.uni)
#[1] 57028

length(is.multi)
#[1] 94

##probes with any repeat overlap
unames<-unique(names(tmp))

tabsav<-table(rmsk38[unlist(tmp),"repFamily"])

tabsav
     # 5S-Deu-L2            Alu          centr            CR1            DNA           DNA?        Dong-R4 
     #        53           7715            206            776             43             30              4 
     #      ERV1          ERV1?           ERVK           ERVL      ERVL-MaLR          ERVL?          Gypsy 
     #      3719             15            150           2490           3489             38            206 
     #    Gypsy?            hAT         hAT-Ac  hAT-Blackjack    hAT-Charlie     hAT-Tip100    hAT-Tip100? 
     #        76             79             66            205           3266            659             43 
     #      hAT?       Helitron      Helitron?             L1             L2 Low_complexity            LTR 
     #        20             26              8           7402           8095           1294             14 
     #      LTR?         Merlin            MIR      MULE-MuDR       Penelope       PiggyBac      PiggyBac? 
     #        61              2           9376             19              9             17             14 
     #       RNA           rRNA       RTE-BovB          RTE-X      Satellite          scRNA  Simple_repeat 
     #        10             64            149            163            187             10           4866 
     #     SINE?          snRNA         srpRNA            SVA  TcMar-Mariner      TcMar-Tc2   TcMar-Tigger 
     #         1             41             23             83             74             76           1420 
     #    TcMar?           telo           tRNA       tRNA-Deu       tRNA-RTE        Unknown 
     #         7             52            157              9             85             54 

##remove those with multioverlap completely..
tmp<-tmp[unlist(lapply(tmp,length))==1]

table(unlist(lapply(tmp,length)))
#     1 
# 57028 

table(rmsk38[unlist(tmp),"repFamily"])
     # 5S-Deu-L2            Alu          centr            CR1            DNA           DNA?        Dong-R4 
     #        53           7712            204            773             39             30              4 
     #      ERV1          ERV1?           ERVK           ERVL      ERVL-MaLR          ERVL?          Gypsy 
     #      3717             15            150           2488           3480             38            206 
     #    Gypsy?            hAT         hAT-Ac  hAT-Blackjack    hAT-Charlie     hAT-Tip100    hAT-Tip100? 
     #        76             79             65            205           3260            659             43 
     #      hAT?       Helitron      Helitron?             L1             L2 Low_complexity            LTR 
     #        20             26              8           7395           8070           1283             14 
     #      LTR?         Merlin            MIR      MULE-MuDR       Penelope       PiggyBac      PiggyBac? 
     #        61              2           9330             19              9             17             14 
     #       RNA           rRNA       RTE-BovB          RTE-X      Satellite          scRNA  Simple_repeat 
     #        10             64            143            163            179             10           4823 
     #     SINE?          snRNA         srpRNA            SVA  TcMar-Mariner      TcMar-Tc2   TcMar-Tigger 
     #         1             41             23             80             74             76           1417 
     #    TcMar?           telo           tRNA       tRNA-Deu       tRNA-RTE        Unknown 
     #         7             52            155              9             83             54 

##lost entries?
all(names(tabsav)==names(table(rmsk38[unlist(tmp),"repFamily"])))
#[1] TRUE

tabsav-table(rmsk38[unlist(tmp),"repFamily"])
     # 5S-Deu-L2            Alu          centr            CR1            DNA           DNA?        Dong-R4 
     #         0              3              2              3              4              0              0 
     #      ERV1          ERV1?           ERVK           ERVL      ERVL-MaLR          ERVL?          Gypsy 
     #         2              0              0              2              9              0              0 
     #    Gypsy?            hAT         hAT-Ac  hAT-Blackjack    hAT-Charlie     hAT-Tip100    hAT-Tip100? 
     #         0              0              1              0              6              0              0 
     #      hAT?       Helitron      Helitron?             L1             L2 Low_complexity            LTR 
     #         0              0              0              7             25             11              0 
     #      LTR?         Merlin            MIR      MULE-MuDR       Penelope       PiggyBac      PiggyBac? 
     #         0              0             46              0              0              0              0 
     #       RNA           rRNA       RTE-BovB          RTE-X      Satellite          scRNA  Simple_repeat 
     #         0              0              6              0              8              0             43 
     #     SINE?          snRNA         srpRNA            SVA  TcMar-Mariner      TcMar-Tc2   TcMar-Tigger 
     #         0              0              0              3              0              0              3 
     #    TcMar?           telo           tRNA       tRNA-Deu       tRNA-RTE        Unknown 
     #         0              0              2              0              2              0 

rm(tabsav)

##get rows in Rep for each Cpg
tmp2<-rmsk38[unlist(tmp),"repClass"]
names(tmp2)<-names(tmp)

tmp3<-rmsk38[unlist(tmp),"repFamily"]
names(tmp3)<-names(tmp)

tmp4<-paste(tmp2,tmp3,sep="|")
names(tmp4)<-names(tmp)

# ##remove those with too few entries
# kill<-table( tmp4 )
# kill<-names(kill)[kill>=1000]

tmp5<-rmsk38[unlist(tmp),"repName"]
names(tmp5)<-names(tmp)

tmp4<-paste(tmp4,tmp5,sep="|")
names(tmp4)<-names(tmp)

length(tmp4)
#[1] 57028

head(tmp4)
#              cg00000108              cg00000109              cg00000884              cg00001099 
#         "SINE|MIR|MIRb"         "LINE|L1|L1MC4"       "LTR|ERVK|MER9a3"       "LTR|ERV1|MER67D" 
#              cg00001269              cg00001534 
# "DNA|hAT-Charlie|MER20"           "LINE|L2|L2b" 

head(sub("(.+\\|.+)\\|.+$","\\1",tmp4) )
       # cg00000108        cg00000109        cg00000884        cg00001099        cg00001269        cg00001534 
       # "SINE|MIR"         "LINE|L1"        "LTR|ERVK"        "LTR|ERV1" "DNA|hAT-Charlie"         "LINE|L2" 

#tmp4<-tmp4[ sub("(.+\\|.+)\\|.+$","\\1",tmp4) %in% kill]

length(tmp4)
#[1] 57028

##add annotations
annoObj$hasAnyRepeatOverlap<-0
annoObj[unique(c(is.multi,is.uni)),"hasAnyRepeatOverlap"]<-1

annoObj$hasMultiRepeatOverlap<-0
annoObj[is.multi,"hasMultiRepeatOverlap"]<-1

annoObj$hasOneRepeatOverlap<-0
annoObj[is.uni,"hasOneRepeatOverlap"]<-1

annoObj$hasRepeatOverlap<-0
annoObj[names(tmp4),"hasRepeatOverlap"]<-1

annoObj$repeatClass<-""
annoObj[names(tmp4),"repeatClass"]<-sub("(.+)\\|.+\\|.+$","\\1",tmp4)

annoObj$repeatFamily<-""
annoObj[names(tmp4),"repeatFamily"]<-sub(".+\\|(.+)\\|.+$","\\1",tmp4)

annoObj$repeatName<-""
annoObj[names(tmp4),"repeatName"]<-sub(".+\\|.+\\|(.+)$","\\1",tmp4)

str(annoObj)
# 'data.frame':   421368 obs. of  62 variables:
#  $ illuminaID                   : chr  "cg21870274" "cg08258224" "cg16619049" "cg18147296" ...
#  $ chr                          : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start                        : num  69591 864703 870161 877159 898803 ...
#  $ end                          : num  69592 864704 870162 877160 898804 ...
#  $ hasUCSCknownGeneOverlap      : num  1 0 1 1 0 0 0 0 0 0 ...
#  $ nameUCSCknownGeneOverlap     : chr  "OR4F5" "" "FAM41C" "FAM41C" ...
#  $ numberUCSCknownGeneOverlap   : num  1 0 1 1 0 0 0 0 0 0 ...
#  $ hasUCSCknownGenePromOverlap  : num  0 0 0 1 0 0 0 0 0 0 ...
#  $ hasUCSCknownGeneUp5kbOverlap : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hasUCSCknownGeneDn5kbOverlap : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ hasUCSCknownTxPromOverlap    : num  0 0 1 1 0 0 0 0 0 0 ...
#  $ hasUCSCknownTxUp5kbOverlap   : num  0 0 1 0 0 0 0 1 1 1 ...
#  $ hasUCSCknownTxDn5kbOverlap   : num  1 1 0 1 0 0 0 0 0 0 ...
#  $ ucscKnownGeneIsDistal        : num  0 1 1 0 1 1 1 1 1 1 ...
#  $ ucscKnownGeneIsGenic         : num  1 0 1 1 0 0 0 0 0 0 ...
#  $ ucscKnownGeneIsDistalNonGenic: num  0 1 0 0 1 1 1 1 1 1 ...
#  $ ucscKnownGeneIsGeneBody      : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ hasGeneOverlap               : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ nameGeneOverlap              : chr  "ENSG00000186092|OR4F5" "" "ENSG00000230368|FAM41C" "" ...
#  $ numberGeneOverlap            : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ hasUp5kbOverlap              : num  0 0 1 0 0 0 0 1 1 1 ...
#  $ nameUp5kbOverlap             : chr  "" "" "ENSG00000234711|TUBB8P11" "" ...
#  $ numberUp5kbOverlap           : num  0 0 1 0 0 0 0 1 1 1 ...
#  $ hasDn5kbOverlap              : num  1 0 0 1 0 0 0 0 0 0 ...
#  $ nameDn5kbOverlap             : chr  "ENSG00000186092|OR4F5" "" "" "ENSG00000234711|TUBB8P11" ...
#  $ numberDn5kbOverlap           : num  1 0 0 1 0 0 0 0 0 0 ...
#  $ hasPromOverlap               : num  0 0 0 1 0 0 0 0 0 0 ...
#  $ namePromOverlap              : chr  "" "" "" "ENSG00000230368|FAM41C" ...
#  $ numberPromOverlap            : num  0 0 0 1 0 0 0 0 0 0 ...
#  $ hasAtacOverlap               : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ nameAtacOverlap              : chr  "" "" "chr1:869670-870169" "" ...
#  $ numberAtacOverlap            : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ isDistal                     : num  0 1 0 0 1 1 1 0 0 0 ...
#  $ isGenic                      : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ isDistalNonGenic             : num  0 1 0 0 1 1 1 0 0 0 ...
#  $ isGeneBody                   : num  1 0 1 0 0 0 0 0 0 0 ...
#  $ weberOE                      : num  0.248 0.249 0.637 0.163 0.565 ...
#  $ weberClass                   : chr  "LCP" "LCP" "HCP" "LCP" ...
#  $ saxonovOE                    : num  0.238 0.268 0.473 0.202 0.441 ...
#  $ saxonovClass                 : chr  "LCG" "LCG" "HCG" "LCG" ...
#  $ cgCount300Bp                 : int  6 2 27 4 9 13 16 8 18 16 ...
#  $ cgPer100Bp                   : num  2 0.667 9 1.333 3 ...
#  $ isCgi                        : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ isShore                      : num  0 1 0 0 0 0 0 0 1 1 ...
#  $ isOcean                      : num  1 0 0 1 1 1 1 1 0 0 ...
#  $ cgiClass                     : chr  "ocean" "shore" "cgi" "ocean" ...
#  $ featureClass                 : chr  "proximal dn" "distal" "proximal up" "promoter" ...
#  $ featureClassUcsc             : chr  "proximal dn" "distal" "distal body" "promoter" ...
#  $ illuminaCpG_CpH_Probe        : chr  "cg" "cg" "cg" "cg" ...
#  $ illuminaDesignType           : chr  "II" "II" "I" "II" ...
#  $ encodeCre                    : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ encodeCreCategory            : chr  "" "" "PLS" "" ...
#  $ encodeCreCtcf                : num  0 0 1 0 0 0 0 0 0 0 ...
#  $ encodeChipMeanNcellPeaks     : num  0 0 4.73 0 1 ...
#  $ encodeChipUniqueTfPeaks      : num  0 0 30 0 9 7 7 1 1 0 ...
#  $ hasAnyRepeatOverlap          : num  0 0 0 0 0 0 0 0 1 1 ...
#  $ hasMultiRepeatOverlap        : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hasOneRepeatOverlap          : num  0 0 0 0 0 0 0 0 1 1 ...
#  $ hasRepeatOverlap             : num  0 0 0 0 0 0 0 0 1 1 ...
#  $ repeatClass                  : chr  "" "" "" "" ...
#  $ repeatFamily                 : chr  "" "" "" "" ...
#  $ repeatName                   : chr  "" "" "" "" ...

rm(tmp,tmp2,tmp3,tmp4,tmp5,is.uni,is.multi)

save(annoObj,file=paste0(HOME,"/annotateFeatures/object_noChromFilter_",nrow(annoObj),"x",ncol(annoObj),"_expandedAnnotations.RData"))

q("no")
###END