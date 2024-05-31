#' Run this script via Rscript with two arguments
#' This script runs in 2 parts, chosen by the first argument (integer 1,2)
#' Part 1: Simulate SNPs under Kingman coalescent w.
#' exp. growth by sampling from expected SFS
#' Part 2: Run ABC parameter estimation for growth rate
#' Both parts use parallelisation
#' Second argument (integer): Sets a seed for replication of random parts                        

args1 <- commandArgs(TRUE) 
part <- as.integer(args1[1])
coh <- as.integer(args1[2]) #used for seed setting


#' R libraries
if (part==1){
#parallelisation 
library(future)
library(future.apply)
  }
if (part==2){
library(abcrf)
  }


ncores1 <- 3
#' number of cores used
mc1 <- switch(part,1,ncores1,10)


#' window size for measuring and simulating genetic diversity
winsize <- 250000
stepsize <- 50000



if (part==1){
#' Simulate via sampling from expected SFS
#' set simulation parameters

nsims <- 30000 #number of sims
#' get parameters as observed from the data
load(paste0("data/mitis_syn_obsdiv",winsize/1000,"K.RData"))


svec <- obs_stats$S[-1] #first entry is genomewide

n <- obs_stats$n[1]
set.seed(coh+14)

#' Sample typical number of mutations: uniform within the counts of seg sites of 
#' two randomly drawn windows
#' currently not used, use controlled by boolean twowinsamp
#' currently used: sample a window at random w. replacement, use its seg sites
twowinsamp <- FALSE
if (twowinsamp){
sampf <- function(){w2 <- sample(length(svec),2,replace = FALSE)
return(sample(svec[w2[1]]:svec[w2[2]],1))}} else {
sampf <- function(){sample(svec,1)}  
}

drawS <- replicate(nsims,sampf())

#' get expected genome-wide SFS to draw SNPs from
#' pregenerated files (the Phi tables from 
#' https://doi.org/10.1371/journal.pgen.1010677,
#' corresponding code repository 
#' https://github.com/fabfreund/usfs_mmc

nsfs <- read.delim("data/MinAlpha-2-MaxAlpha-2-NoStepAlpha-0-MinRho-0-MaxRho-25-NoStepRho-250-noSamples-55_expectedPhi.txt", header=FALSE)
#' remove the first column (total length of coalescent)
nsfs <- nsfs[,-1]
#' growth rates are as in the pregenerated Phi table - discrete!
growth_rates <- seq(0,25,0.1)
#' We just draw uniformly
g_params <- sample(growth_rates,size = nsims,replace = TRUE)
#' we name the rows of the Phi table by the growth rates they are 
#' describing
rownames(nsfs) <- growth_rates

#' function to simulate SNP sets
seq_sim <-  function(g,s1){
  snp_counts <- sample(x = 1:ncol(nsfs),
                       size = s1,
                       replace = TRUE,
                       prob = nsfs[as.character(g),])
  snp_counts <- table(snp_counts)
  seq1 <- NULL
  for (char_class in names(snp_counts)){
    temp <- replicate(snp_counts[char_class],
                      sample(c(rep(1,as.integer(char_class)),
                               rep(0,n-as.integer(char_class)))))   
    seq1 <- cbind(seq1,temp)
    seq1 <- seq1[,sample(ncol(seq1))]
  }
  return(seq1)
}
#' Genetic diversity statistics
#' Get 19 minor allele frequency quantiles
allele_freqs <- function(seq1,
                         quant_v = seq(.05,.95,.05)){
  if (!any(is.na(seq1))){
  if (!(is.matrix(seq1))){as.matrix(seq1)}  
  temp1 <- colSums(seq1)/nrow(seq1)
  temp1 <- sapply(temp1,function(x){min(x,1-x)})
  temp2 <- quantile(temp1,quant_v)
  names(temp2) <- paste0("AF_qu",quant_v)
  return(temp2)
}
}

#' Get pairwise Hamming distances, but normalized to max Hamming observed
hammfun <- function(seq1){
  if (is.matrix(seq1)){
    dist1 <- dist(seq1,method="manhattan")
    dist1 <- dist1/max(dist1) #to counter-act issues w. 
  } else {dist1 <- 0}
  quant_v = c(.1,.3,.5,.7,.9)
  gendist <- quantile(dist1,quant_v)
  names(gendist) <- paste0("hamm_q",seq(1,9,2))
  return(gendist)
}


#' Get Tajima's D, from https://raw.githubusercontent.com/sdwfrost/popseq/master/R/sfsR.R
sfsR=function(hapmatrix){
  n=dim(hapmatrix)[1]
  S=dim(hapmatrix)[2]
  I=1:(n-1)
  a1=sum(1/I)
  a2=sum(1/I^2)
  b1=(n+1)/3/(n-1)
  b2=2*(n^2+n+3)/9/n/(n-1)
  c1=b1-1/a1
  c2=b2-(n+2)/a1/n+a2/a1^2
  e1=c1/a1
  e2=c2/(a1^2+a2)
  Dnum=sum(hamming.distance(hapmatrix))/n/(n-1)-S/a1
  Dden=sqrt(e1*S+e2*S*(S-1))
  return(c("TajD"=Dnum/Dden))
}

require(e1071)
tajd <- function(seq1){sfsR(seq1+1)}



stats1 <- function(seq1){c("n"=nrow(seq1),"S"=ncol(seq1),allele_freqs(seq1),hammfun(seq1),
                           tajd(seq1))}

#' Plan for parallelisation
plan(strategy = multisession, workers=mc1)

#' Simulating in parallel
sims1 <- future_sapply((1:nsims),
                function(x){cat("iteration",x,"\n")
                            stats1(seq_sim(g = g_params[x],
                                           s1=drawS[x]))},
                future.scheduling=TRUE,future.seed=TRUE)

save(sims1,g_params,drawS,
     file = paste0("sims_mitis_PRF",coh,"_",winsize/1000,"K.RData"),
     compress = "xz")
}

if (part == 2){
#' Load observed and simulated diversity 
load(paste0("data/mitis_syn_obsdiv",winsize/1000,"K.RData"))
obs_stats <- obs_stats[-1,] #kick out genomewide, no clear drop-off in windows (there are some with clearly smaller S, though)
obs_stats <- obs_stats[,-1] #kick out sample size


load(paste0("sims_mitis_PRF",coh,"_",winsize/1000,"K.RData"))
sims2 <- t(sims1)
g_params2 <- g_params

sims2 <- sims2[,-1]

#' Parameter estimation Kingman+exp via ABCrf
sumstat1 <- data.frame(g_params2,sims2)
growth_rf <- regAbcrf(g_params2~.,
                      data=sumstat1,
                      paral = TRUE,ncores=mc1)
growthfit_rf <- predict(growth_rf,as.data.frame(obs_stats),
                        sumstat1,paral=TRUE,ncores=mc1,
                        quantiles = seq(0,1,0.025)) 

save(growthfit_rf,growth_rf,file=paste0("res_100winsmitisPRF_",winsize/1000,"K.RData"),compress = "xz")
}
