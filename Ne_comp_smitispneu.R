#' Parameters for both species 
load("data/snps_in_neut.RData") #from alignment
mode1 <- c("neutral","all")[2] #Base calculations on all SNPs or only non-synonymous
 
n <- list("spneu"=315,"smitis"=55) #sample sizes
coregenomesize <- list("spneu"=snps_lengths$spneu$corel[[mode1]],
                       "smitis"=snps_lengths$smitis$corel[[mode1]]) #only ne
mu_year <- list("spneu"=1.57*10^(-6),"smitis"=1.57*10^(-6)) # Per-site, per-year mut rate 
gen_py <- list("spneu"=14,"smitis"=14) #Number of generations per year

#Number of SNPs 
S_allcore <- list("spneu"=snps_lengths$spneu$noSNPs[[mode1]],
                  "smitis"=snps_lengths$smitis$noSNPs[[mode1]]) 

#Function to estimate Ne
compute_Ne <- function(mu_pg,EL,S1){S1/(mu_pg*EL)}

#' Compute for Spneu
mu_gen <- mu_year$spneu*coregenomesize$spneu/gen_py$spneu

#Read in ABC results
spneu_res_table <- read.csv("data/abc_results_spneu.csv")
g_estims <- spneu_res_table$median_PRF

g <- c("median"=median(g_estims),"lowmed"=min(g_estims),
       "highmed"=max(g_estims),"q25"=quantile(g_estims,.25),
        "q75"=quantile(g_estims,.75))

#Use recursive approach to compute expected branch lengths
#See in-script documentation in ext_fun_CREB.R
source("ext_fun_CREB.R")
ELn <- sapply(g,function(x){sum(REBiepg(n$spneu,x))})

Ne_spneu <- sapply(ELn,function(l1){compute_Ne(mu_pg=mu_gen,l1,
                                               S1=S_allcore$spneu)})
#' Compute for Smitis
mu_gen <- mu_year$smitis*coregenomesize$smitis/gen_py$smitis

#' Read in APC results
smitis_res_table <- read.csv("data/abc_results_smitis.csv")
g_estims <- smitis_res_table$median_PRF
g <- c("median"=median(g_estims),"lowmed"=min(g_estims),
       "highmed"=max(g_estims),"q25"=quantile(g_estims,.25),
       "q75"=quantile(g_estims,.75))
#source("ext_fun_CREB.R")
ELn <- sapply(g,function(x){sum(REBiepg(n$smitis,x))})

Ne_smitis <- compute_Ne(mu_pg=mu_gen,EL = ELn,S1=S_allcore$smitis)

Nes_both <- list(Ne_smitis=Ne_smitis,Ne_spneu=Ne_spneu)
save(Nes_both,file = paste0("Ne_estimates_",mode1,".RData"))
