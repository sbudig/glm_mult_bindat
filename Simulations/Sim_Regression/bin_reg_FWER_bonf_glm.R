# Simulation of the FWER in a regression model with the GLM and the Bonferroni method

# Setup -------------------------------------------------------------------

# Functions
source("bin_reg_source.R")

# required packages -------------------------------------------------------

library("tidyverse")

# Simulation Function -----------------------------------------------------
# ntrt: sample size per dose/treatment group
# pis: probabilities of occurrences fo each endpoint
# nresp: number of endpoints
# size: cluster size (=1 for binary data)
# corr: correlation between endpoints
# di: dose groups
# alpha: type 1 error
# This function simulates correlated binary data with the corresponding parameters. 
# A GLM is fitted to the simulated data for each endpoint and 
# then a correction for the p-values of the parameters of interest is performed 
# using bonferroni. The parameters used and the corresponding results to obtain 
# the FWER are returned in a data frame. 
sim_mmm_bin_fwer_bonf_glm <- function(nsim, 
                        ntrt, 
                        pis, 
                        nresp, 
                        size, 
                        corr, 
                        di = c(0, 2.5, 5, 10), 
                        alpha = 0.05){
  
  NTRT <- rep(ntrt, times = length(di))
  names(NTRT) <- di
  
  mutrt <- mi0ed(mi=logit(pis), ndos=length(di), nresp=nresp)
  MUTRT <- t(mutrt)
  
  NRESP <- nrow(mutrt)
  NAMS <- paste("s",1:NRESP, sep="")
  NAMF <- paste("f",1:NRESP, sep="")
  FORMVEC <- paste("cbind(", NAMS, ",",  NAMF, ") ~ dose", sep="")
  
  FORMLIST <- as.list(FORMVEC)
  names(FORMLIST) <- NAMS
  
  SIGMA <- rep(1, NRESP)
  
  VCOV <- equicorrvc(sigma=SIGMA, corr=corr)
  
  print(paste0("pis ", pis, " nresp ", nresp, " ntrt ", ntrt, " corr ", corr))
  
  nrej <- 0
  nfail <- 0
  nreji <- 0
  nerror <- 0
  
  for(i in 1:nsim)
  {
    dati <- reg_simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT, covmat=VCOV, size=size)
    
    
    FITLIST <- FORMLIST
    for(k in 1:NRESP){FITLIST[[k]] <- glm(FORMLIST[[k]], data=dati, family=binomial)}
    
    tryCatch({
      
      SUMLIST <- lapply(FITLIST, summary)
      COEFLIST <- lapply(SUMLIST, coef)
      COEFLIST <- do.call(rbind.data.frame,(COEFLIST))[,4]
      
      SEQ <- seq(from = 2, to = length(COEFLIST), by = 2 )
      pwd <- p.adjust(COEFLIST[SEQ], method = "bonferroni" )
      
      nrej <- nrej + sum(na.omit(pwd<alpha))
      nreji <- nreji + as.integer(any(na.omit(pwd<alpha)))
      nfail <- nfail + sum(is.na(pwd<alpha))
    }, 
    error = function(e){
      nerror <<- nerror + 1
    })
  }
  
  return(data.frame(
    method = "bonf",
    model = "glm",
    nsim = nsim, pis = pis, nresp = nresp, 
    ntrt = ntrt, corr = corr, 
    nfail = nfail, nerror = nerror, 
    nreji = nreji, nrej = nrej))
}

system.time(sim_mmm_bin_fwer_bonf_glm_res <- do.call(rbind,apply(as.matrix(simdat_reg_FWER), 1, 
                                                   function(x){sim_mmm_bin_fwer_bonf_glm(nsim=unname(x[1]),pis=unname(x[2]),nresp=unname(x[3]), 
                                                                           ntrt=unname(x[4]),size=unname(x[5]),corr=unname(x[6]))})))

write.csv(sim_mmm_bin_fwer_bonf_glm_res, ".\\intermediate_results\\bin_reg_FWER_bonf_glm_1.csv", row.names = FALSE)



