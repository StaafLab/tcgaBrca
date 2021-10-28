#####======================================================================#####
### Correct TCGA BRCA beta values for infiltration
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

if(!requireNamespace("flexmix", quietly = TRUE)) {
  install.packages("flexmix") }

library("flexmix")

if(!requireNamespace("doParallel", quietly = TRUE)) {
	install.packages("doParallel") }

library(doParallel)

if(!requireNamespace("parallel", quietly = TRUE)) {
	install.packages("parallel") }

library(parallel)

if(!requireNamespace("tools", quietly = TRUE)) {
	install.packages("tools") }

library(tools)

################################################################################
##Get core tumor set from methylation data

load(file=paste0(HOME,"/me/object_flexmixAdjustmentObject.RData"))

ls()
# [1] "adjustBeta" "betaNew"    "fracTum"    "HOME"       "MANIFEST"  
# [6] "TUMOR_TYPE" "uids"    

################################################################################
###match beta data to uids

str(betaNew)
 # num [1:421368, 1:630] 0.718 0.766 0.189 0.925 0.868 0.642 0.95 1 0.986 1 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:421368] "cg21870274" "cg08258224" "cg16619049" "cg18147296" ...
 #  ..$ : chr [1:630] "TCGA-3C-AAAU-01A" "TCGA-3C-AALI-01A" "TCGA-3C-AALJ-01A" "TCGA-3C-AALK-01A" ...

################################################################################
###using multicore

no_cores <- detectCores(logical = TRUE)

cat("using", no_cores-3,"cores","\n")

cl <- makeCluster(no_cores-3)  
registerDoParallel(cl)  

##estimated runtime depends on number of cores mainly..

##home i7-6700K 4core/8thread

##work threadripper 16core/32thread

clusterEvalQ(cl, {
  library("flexmix")
})
# [[1]]
# [1] "flexmix"   "lattice"   "stats"     "graphics"  "grDevices" "utils"     "datasets" 
# [8] "methods"   "base"     

# [[2]]
# [1] "flexmix"   "lattice"   "stats"     "graphics"  "grDevices" "utils"     "datasets" 
# [8] "methods"   "base"     

# [[3]]
# [1] "flexmix"   "lattice"   "stats"     "graphics"  "grDevices" "utils"     "datasets" 
# [8] "methods"   "base"     

# [[4]]
# [1] "flexmix"   "lattice"   "stats"     "graphics"  "grDevices" "utils"     "datasets" 
# [8] "methods"   "base"     

# [[5]]
# [1] "flexmix"   "lattice"   "stats"     "graphics"  "grDevices" "utils"     "datasets" 
# [8] "methods"   "base"     

##add a unique seed to each row to make run deterministic
	##ensures all runs get unique seed
betaRun<-cbind(seed=1:nrow(betaNew),betaNew)
betaNames<-colnames(betaNew)

##sync all parallel sessions
	##will use same start seed for all sessions -> could be problematic as small beta noise factor will repeat across sessions
#clusterSetRNGStream(cl, 20200918) 
system.time(
	results<-parRapply(cl = cl, betaRun, adjustBeta,purity=fracTum,snames=betaNames,seed=TRUE)
)
  #   user   system  elapsed 
  # 145.89    30.34 46493.75 

table(unlist(lapply(results,function(x) x$n.groups)))
 #     1      2      3 
 # 69291 100990 251087 

rm(betaRun,betaNames)

save(results,file=paste0(HOME,"/me/object_flexmixAdjustmentOutput.RData"))

################################################################################

str(results[[1]])
# List of 10
#  $ y.norm          : Named num [1:649] 0.585 0.491 0.644 0.812 0.605 0.695 0.687 0.613 0.632 0.485 ...
#   ..- attr(*, "names")= chr [1:649] "9993943013_R04C01" "9993943013_R01C02" "9993943005_R02C02" "9993943017_R01C01" ...
#  $ y.tum           : Named num [1:649] 0.749 0.416 0.808 0.976 0.769 0.62 0.852 0.777 0.557 0.41 ...
#   ..- attr(*, "names")= chr [1:649] "9993943013_R04C01" "9993943013_R01C02" "9993943005_R02C02" "9993943017_R01C01" ...
#  $ y.orig          : Named num [1:649] 0.718 0.443 0.731 0.912 0.682 0.655 0.824 0.69 0.599 0.437 ...
#   ..- attr(*, "names")= chr [1:649] "9993943013_R04C01" "9993943013_R01C02" "9993943005_R02C02" "9993943017_R01C01" ...
#  $ groups          : int [1:649] 1 2 1 1 1 2 1 1 2 2 ...
#  $ n.groups        : int 2
#  $ med.norm        : num 0.651
#  $ glob.cor        : num 0.037
#  $ avg.betaDiff    : num -0.03
#  $ model.intercepts: num [1:2] 0.867 0.472
#  $ model.slopes    : num [1:2] -0.165 0.075

x<-unlist(lapply(results,function(x) x$model.intercepts))

y<-unlist(lapply(results,function(x) x$model.slopes))

z<-unlist(lapply(results,function(x) x$glob.cor))

k<-unlist(lapply(results,function(x) x$avg.betaDiff))

sum(is.na(x))
#[1] 0
sum(is.na(y))
#[1] 0
sum(is.na(z))
#[1] 0
sum(is.na(k))
#[1] 0

pdf(file=paste0(HOME,"/me/plot_flexmixAdjustment_coefficientsAndAdjustmentStats.pdf"),width=8,height=8,useDingbats=FALSE)

set.seed(12345)
r<-sample(1:length(x),15000)
plot(x[r],y[r],pch=16,cex=.3,main="intercept vs slope",xlab="intercept",ylab="slope",sub="15000 random intercept")

hist(y[r],breaks=51,main="adjustment slopes",sub="15000 random slopes",xlab="slope")

hist(x[r],breaks=51,main="adjustment intercepts",sub="15000 random intercept",xlab="intercept")

plot(abs(z[r]),abs(k[r]),pch=16,cex=.3,
		main="",
		xlab="abs correlation Meth-Purity",ylab="abs delta Meth post-pre correction",
		sub="15000 random CpGs")

dev.off()

rm(x,y,z,k,r)

################################################################################
###Extract corrected matrix

##adjusted
betaAdj<-do.call("rbind",lapply(results,function(x) x$y.tum))

##unadjusted
betaNorm<-do.call("rbind",lapply(results,function(x) x$y.norm))

##unadjusted
betaOrig<-do.call("rbind",lapply(results,function(x) x$y.orig))

##save three main objects

##corrected betas - tumor
con<-gzfile(paste0(HOME,"/me/data450k_",nrow(betaAdj),"x",ncol(betaAdj),"_minfiNormalized_ringnerAdjusted_purityAdjusted_tumorBetaValues.txt.gz"),"w")
write.table(betaAdj,file=con,sep="\t",quote=F)
close(con)
rm(con)

md5sum(paste0(HOME,"/me/data450k_",nrow(betaAdj),"x",ncol(betaAdj),"_minfiNormalized_ringnerAdjusted_purityAdjusted_tumorBetaValues.txt.gz"))

##corrected betas - normal
con<-gzfile(paste0(HOME,"/me/data450k_",nrow(betaAdj),"x",ncol(betaNorm),"_minfiNormalized_ringnerAdjusted_purityAdjusted_normalBetaValues.txt.gz"),"w")
write.table(betaNorm,file=con,sep="\t",quote=F)
close(con)
rm(con)

md5sum(paste0(HOME,"/me/data450k_",nrow(betaAdj),"x",ncol(betaNorm),"_minfiNormalized_ringnerAdjusted_purityAdjusted_normalBetaValues.txt.gz"))

##corrected betas - normal
con<-gzfile(paste0(HOME,"/me/data450k_",nrow(betaAdj),"x",ncol(betaOrig),"_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.txt.gz"),"w")
write.table(betaOrig,file=con,sep="\t",quote=F)
close(con)
rm(con)

md5sum(paste0(HOME,"/me/data450k_",nrow(betaAdj),"x",ncol(betaOrig),"_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.txt.gz"))

##as R-obj
save(betaAdj,file=paste0(HOME,"/me/data450k_",nrow(betaAdj),"x",ncol(betaAdj),"_minfiNormalized_ringnerAdjusted_purityAdjusted_tumorBetaValues.RData"))
save(betaNorm,file=paste0(HOME,"/me/data450k_",nrow(betaAdj),"x",ncol(betaNorm),"_minfiNormalized_ringnerAdjusted_purityAdjusted_normalBetaValues.RData"))
save(betaOrig,file=paste0(HOME,"/me/data450k_",nrow(betaAdj),"x",ncol(betaOrig),"_minfiNormalized_ringnerAdjusted_purityAdjusted_originalBetaValues.RData"))

write.table(cbind(sample=colnames(betaNew),tumor_fraction=fracTum),
	file=paste0(HOME,"/me/samplePurityVectorForFlexmixCorrection.txt"),
	sep="\t",quote=F,row.names=FALSE)

##get groupsings per CpG
betaAdjGroups<-do.call("rbind",lapply(results,function(x) x$groups))
colnames(betaAdjGroups)<-colnames(betaAdj)
rownames(betaAdjGroups)<-rownames(betaAdj)

save(betaAdjGroups,file=paste0(HOME,"/me/data450k_",nrow(betaAdj),"x",ncol(betaAdjGroups),"_minfiNormalized_ringnerAdjusted_purityAdjusted_tumorBetaGroups.RData"))

q("no")
###END