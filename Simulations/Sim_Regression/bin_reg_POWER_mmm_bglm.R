# Simulation of the Power in a regression model with the BayesGLM and the MMM method

# Setup -------------------------------------------------------------------

# Functions
source("bin_reg_source.R")

# required packages -------------------------------------------------------

library("multcomp")
library("tidyverse")
library("arm") # bayesglm

# Simulation Function -----------------------------------------------------
# ntrt: sample size per dose/treatment group
# pi0: intercept
# nep_h0: number of endpoints under H0
# nep_ha: number of endpoints under HA
# slope: slope parameter of regression
# size: cluster size (=1 for binary data)
# corr: correlation between endpoints
# di: dose groups
# alpha: type 1 error
# This function simulates correlated binary data with the corresponding parameters. 
# A BayesGLM is fitted to the simulated data for each endpoint and 
# then a correction for the p-values of the parameters of interest is performed 
# using the MMM approach. The parameters used and the corresponding 
# results to obtain the Power are returned in a data frame. 
sim_mmm_bin_power_mmm_bglm <- function(nsim, 
                        ntrt, 
                        pi0, 
                        nep_h0, 
                        nep_ha, 
                        slope, 
                        size, 
                        corr, 
                        di = c(0, 2.5, 5, 10), 
                        alpha = 0.05){
  
  NTRT <- rep(ntrt, times = length(di))
  names(NTRT) <- di
  
  eta <- logit(pi0)+slope*di
  
  eta_h0 <- matrix(rep(logit(pi0), each = length(di), times = nep_h0), ncol = length(di))
  eta_ha <- matrix(rep(eta, each= nep_ha ), ncol = length(di))
  
  MUTRT_h0 <- t(eta_h0)
  MUTRT_ha <- t(eta_ha)
  
  NRESP <- nep_h0+nep_ha
  NAMS <- paste("s",1:NRESP, sep="")
  NAMF <- paste("f",1:NRESP, sep="")
  FORMVEC <- paste("cbind(", NAMS, ",",  NAMF, ") ~ dose", sep="")
  
  FORMLIST <- as.list(FORMVEC)
  names(FORMLIST) <- NAMS
  
  SIGMA_h0 <- rep(1, nep_h0)
  SIGMA_ha <- rep(1, nep_ha)
  
  VCOV_h0 <- equicorrvc(sigma=SIGMA_h0, corr=corr)
  VCOV_ha <- equicorrvc(sigma=SIGMA_ha, corr=corr)
  
  print(paste0("slope ", slope, " pi0 ", pi0, " nep_h0 ", nep_h0, " nep_ha ", nep_ha, " ntrt ", ntrt, " corr ", corr))
  
  seq_h0 <- seq(from = 1, to = nep_h0, by = 1)
  seq_ha <- seq(from = nep_h0+1, to = nep_h0+nep_ha, by = 1)
  
  nrej <- 0   # How many rejected per Dataset
  nreji <- 0  # At least one rejected per Dataset
  nrejcor <- 0 # how often all h0 are true
  
  npow <- 0   # How many differences found per Dataset
  npowi <- 0  # At least one difference found per Dataset
  npowcor <- 0 # How often find both differences
  
  nfail <- 0  
  nerror <- 0
  
  for(i in 1:nsim)
  {
    
    dat_h0 <- reg_pow_simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_h0, covmat=VCOV_h0, size=size, hyp = "h0", neph0 = nep_h0)
    dat_ha <- reg_pow_simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_ha, covmat=VCOV_ha, size=size, hyp = "ha", neph0 = nep_h0)
    
    dati <- cbind(dat_h0,dat_ha)
    #dplyr::select(dati, starts_with("s"))
    
    
    FITLIST <- FORMLIST
    for(k in 1:NRESP){FITLIST[[k]] <- bayesglm(FORMLIST[[k]], data=dati, family=binomial)}
    
    class(FITLIST) <- "mmm"
    
    tryCatch({
      XT <- seq(from=2, to= 2*NRESP, by = 2)
      K <- diag(length(coef(FITLIST)))[XT,]
      rownames(K) <- names(coef(FITLIST))[XT]
      
      GF <- glht(FITLIST, linfct = K)
      SGF <- summary(GF)
      
      nrej <- nrej + sum(na.omit(SGF$test$pvalues[seq_h0]<alpha))
      nreji <- nreji + as.integer(any(na.omit(SGF$test$pvalues[seq_h0]<alpha)))
      nrejcor <- nrejcor + as.integer((sum(na.omit(SGF$test$pvalues[seq_h0]<alpha)) == 0 ))
      
      npow <- npow + sum(na.omit(SGF$test$pvalues[seq_ha]<alpha))
      npowi <- npowi + as.integer(any(na.omit(SGF$test$pvalues[seq_ha]<alpha)))
      npowcor <- npowcor + as.integer((sum(na.omit(SGF$test$pvalues[seq_ha]<alpha)) == nep_ha))
      
      nfail <- nfail + sum(is.na(SGF$test$pvalues))
    }, 
    error = function(e){
      nerror <<- nerror + 1
    })
  }
  
  return(data.frame(
    method = "mmm",
    model = "bglm",
    nsim = nsim, ntrt = ntrt, pi0 = pi0,
    nep_h0 = nep_h0, nep_ha = nep_ha, 
    slope = slope, corr = corr, 
    nfail = nfail, nerror = nerror, 
    nrej = nrej, nreji = nreji, nrejcor = nrejcor, 
    npow = npow, npowi = npowi, npowcor = npowcor))
}

simdat <- filter(simdat_reg_POWER, (nep_h0+nep_ha)==10)


system.time(sim_mmm_bin_power_mmm_bglm_res <- do.call(rbind,apply(as.matrix(simdat), 1, 
                                                   function(x){sim_mmm_bin_power_mmm_bglm(nsim=unname(x[1]),ntrt=unname(x[2]),pi0=unname(x[3]), 
                                                                           nep_h0=unname(x[4]),nep_ha=unname(x[5]),slope=unname(x[6]),
                                                                           size=unname(x[7]),corr=unname(x[8]))})))

write.csv(sim_mmm_bin_power_mmm_bglm_res, ".\\intermediate_results\\bin_reg_POWER_mmm_bglm_1.csv", row.names = FALSE)


