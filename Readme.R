
## A working example for runing the KL-meta and YETI methods

mypath ="~/Downloads/yetiexample/“
setwd(mypath)

library(Matrix) ## require package when install "metafor" package
#install.packages(“metafor”)
library(metafor) ## use a function "rm" when run YETI method


##==================================================================================================================================
## Load a simualted SNP dataset
## 
## The info about "SNPdata": 
##  (i) L=200, K=2, # of subjects in cases=5000, # of subjects in controls=5000
## (ii) the dominant model is considered, ie, SNPs may take values 0 for two copies of a putative non-risk allele, or 1 otherwise.
##==================================================================================================================================

SNPdata = read.table("SNPdata.txt",header=TRUE)
head(SNPdata)

##====================================================================================================================================
## Run the KL-meta method (only apply for dominant model):
##
##  Input: data (e.g., ID, Y, Z, SNP1,...,SNPL), L (# of markers), K (# of populations/studies), pre.alpha (tuning parameter), alpha
## Output: SNP1 (id for the 1st interacting SNP), SNP2 (id for the 2nd interacting SNP), bet3.overall.est, beta3.overall.se, pvalue   
##====================================================================================================================================

## source the KLmeta funciton
source("KLmeta_func.R")
kl = KLmeta.func(mydata=SNPdata, L=200, K=2, pre.alpha=0.01, alpha=0.05)
kl
#rm(list=ls())


##====================================================================================================================================
## Run the YETI method (only apply for the dominant model):
##
##  Input: data (e.g., ID, Y, Z, SNP1,...,SNPL), L (# of markers), K (# of populations/studies), pre.alpha (tuning parameter), alpha
## Output: SNP1 (id for the 1st interacting SNP), SNP2 (id for the 2nd interacting SNP), bet3.overall.est, beta3.overall.se, pvalue   
##====================================================================================================================================

## source the KLmeta funciton
source("YETI_func.R")
yeti = YETI.func(mydata=SNPdata, L=200, K=2, pre.alpha=0.01, alpha=0.05)
yeti
#rm(list=ls())


##====================================================================================================================================
## Run the YETI2 method (only apply for the dominant model):
##
##  Input: data (e.g., ID, Y, Z, SNP1,...,SNPL), L (# of markers), K (# of populations/studies), pre.alpha (tuning parameter), alpha
## Output: SNP1 (id for the 1st interacting SNP), SNP2 (id for the 2nd interacting SNP), bet3.overall.est, beta3.overall.se, pvalue   
##====================================================================================================================================

## source the YETI2 funciton
source("YETI2_func.R")
yeti2 = YETI2.func(mydata=SNPdata, L=200, K=2, pre.alpha=0.01, alpha=0.05)
yeti2
#rm(list=ls())
