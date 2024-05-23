# Simulation of the FWER in a regression model with the BayesGLM and the MMM method

# Setup -------------------------------------------------------------------

# Functions
source("bin_reg_source.R")

# required packages -------------------------------------------------------

library("multcomp") # MMM method
library("tidyverse")
library("arm") # for bayesglm

# Simulation Function -----------------------------------------------------
# ntrt: sample size per dose/treatment group
# pis: probabilities of occurrences fo each endpoint
# nresp: number of endpoints
# size: cluster size (=1 for binary data)
# corr: correlation between endpoints
# di: dose groups
# alpha: type 1 error

sim_mmm_bin <- function(nsim, 
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
    for(k in 1:NRESP){FITLIST[[k]] <- bayesglm(FORMLIST[[k]], data=dati, family=binomial)}
    
    class(FITLIST) <- "mmm"
    
    tryCatch({
      XT <- seq(from=2, to= 2*NRESP, by = 2)
      K <- diag(length(coef(FITLIST)))[XT,]
      rownames(K) <- names(coef(FITLIST))[XT]
      
      GF <- glht(FITLIST, linfct = K)
      SGF <- summary(GF)
      
      
      nrej <- nrej + sum(na.omit(SGF$test$pvalues<alpha))
      nreji <- nreji + as.integer(any(na.omit(SGF$test$pvalues<alpha)))
      nfail <- nfail + sum(is.na(SGF$test$pvalues))
    }, 
    error = function(e){
      nerror <<- nerror + 1
    })
  }
  
  return(data.frame(
    method = "mmm",
    model = "bglm",
    nsim = nsim, pis = pis, nresp = nresp, 
    ntrt = ntrt, corr = corr, 
    nfail = nfail, nerror = nerror, 
    nreji = nreji, nrej = nrej))
}

system.time(sim_mmm_bin_res <- do.call(rbind,apply(as.matrix(simdat_reg_FWER), 1, 
                                                   function(x){sim_mmm_bin(nsim=unname(x[1]),pis=unname(x[2]),nresp=unname(x[3]), 
                                                                           ntrt=unname(x[4]),size=unname(x[5]),corr=unname(x[6]))})))

write.csv(sim_mmm_bin_res, ".\\intermediate_results\\bin_reg_FWER_mmm_bglm_1.csv", row.names = FALSE)



