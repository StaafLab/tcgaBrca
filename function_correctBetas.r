#####======================================================================#####
### Function for correcting Illumina betas for Lund TNBC data set
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Affiliation: Lund University / Oncology & Pathology

################################################################################
##source from Staaf-lab git-clone

#work
# dir.create("~/hdd1/luTnbc")
# HOME<-"~/hdd1/luTnbc"
# GIT<-"~/Documents/luTnbc"
#home
# dir.create("I:/data/luTnbc")
# HOME<-"I:/data/luTnbc"
# GIT<-"F:/gitProjects/luTnbc"

##source
#source(paste0(GIT,"/function_correctBetas.r"))

################################################################################
##load required packages

if(!requireNamespace("flexmix", quietly = TRUE)) {
  install.packages("flexmix") }

if(! "flexmix" %in% names(sessionInfo()$otherPkgs) ) {
  library("flexmix") } 

################################################################################
##load required packages

##define function with input = betas and purity estimate
  ##output = line parameters for L1/L2/L3
adjustBeta<-function(methylation=NULL,purity=NULL,snames=NULL,nmax=3,nrep=3,seed=TRUE) {
  #grab unique row-seed from vector slot 1
    #delete seed afterwards
  if(seed) {
  set.seed( as.integer(methylation[1]) )
  methylation<-methylation[-1]
  }
  #define variables
  x<-as.numeric(purity)
  x2<-1-as.numeric(purity)
  y<-as.numeric(methylation)
  #calculate global corr
  gl.corr<-suppressWarnings(cor(x,y))
  gl.corr[is.na(gl.corr)]<-0
  gl.corr<-round(gl.corr,3)
  #add small gaussian noise to x - problem with large number of zero samples
  y2<-y+rnorm(length(y),mean=0,sd=.005)  
  #do modeling 
  model <- stepFlexmix(y2 ~ x,k = 1:nmax, nrep = nrep,verbose = FALSE)
  model <- getModel(model, "BIC")
  #get clusters
  cl<-clusters(model)

  #make sure clusters are numbered 1 to 3, odd cases exist where one pop has zero members from flexMix
    #can rename because flexmix object not used after this
  cl<-as.integer(factor(cl))

  ##get line parameters for each pop - calculate from data - use original y
  res.norm<-unlist(lapply(1:nmax,function(z) { 
    if(z %in% cl) {
      m<-lm(y[cl==z]~x[cl==z])
      r<-coefficients(m)[1]+residuals(m)
      names(r)<-snames[cl==z]
      r
    } else { NULL }
  }))
  res.norm<-res.norm[snames]
  ##get line parameters for each pop - calculate from data - use original y
  res.tum<-unlist(lapply(1:nmax,function(z) { 
    if(z %in% cl) {
      m<-lm(y[cl==z]~x2[cl==z])
      r<-coefficients(m)[1]+residuals(m)
      names(r)<-snames[cl==z]
      r
    } else { NULL }
  }))
  res.tum<-res.tum[snames]

  ##in very rare instances flexmix calls 3 populations but one groups has zero members. -> has parameters for non-existant pop in output?!
  #res.int<-round(as.numeric(unlist(lapply(slot(model,"components"),function(z) slot(z[[1]],"parameters")$coef[1]))),3) 
  #res.slope<-round(as.numeric(unlist(lapply(slot(model,"components"),function(z) slot(z[[1]],"parameters")$coef[2]))),3)

  ##get line parameters for each pop - calculate from data - use original y and 1-purity
  res.int<-unlist(lapply(1:nmax,function(z) { 
    if(z %in% cl) {
      m<-lm(y[cl==z]~x2[cl==z])
      r<-coefficients(m)[1]
      round(as.numeric(r),3)
    } else { NULL }
  }))
  ##intercept is never zero but if only one member in flexMix population slope will be NA - fix
    ##if pop only has one member then both tumor and normal estimate will be original beta value
      ##could fix by assigning to nearest pop with >1 member but not certain this any better..
  res.slope<-unlist(lapply(1:nmax,function(z) { 
    if(z %in% cl) {
      m<-lm(y[cl==z]~x2[cl==z])
      r<-coefficients(m)[2]
      if(is.na(r)) { r<-0 } ##fix for small number of NAs - set slope to zero
      round(as.numeric(r),3)
    } else { NULL }
  }))

  ##cap at 0 and 1
  res.tum[res.tum > 1] <- 1
  res.tum[res.tum < 0] <- 0
  res.norm[res.norm > 1] <- 1
  res.norm[res.norm < 0] <- 0

  #round
  res.tum<-round(res.tum,3)
  res.norm<-round(res.norm,3)
  res.orig<-round(methylation,3)

  ##return some stats
  return( list(y.norm=res.norm,
    y.tum=res.tum,
    y.orig=res.orig,
    groups=cl,
    n.groups=length(levels(factor(cl))),
    med.norm=round(median(res.norm),3),
    glob.cor=gl.corr,
    avg.betaDiff=round(mean(res.orig-res.tum),3),
    model.intercepts=res.int,
    model.slopes=res.slope   
    )
  )
}

###END