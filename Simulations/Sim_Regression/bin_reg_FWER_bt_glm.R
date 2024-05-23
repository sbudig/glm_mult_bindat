# Simulation of the FWER in a regression model with the GLM and the bootstrap method

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

sim_mmm_bin <- function(nsim, 
                        ntrt, 
                        pis, 
                        nresp, 
                        size, 
                        corr, 
                        nboot,
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
  
  # for bootstrap DF
  ltrt <- names(NTRT)
  nvectrt <- as.vector(NTRT)
  ftrt <- factor(rep(ltrt, times=nvectrt), levels=ltrt)

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
    
    SUMLIST <- lapply(FITLIST, function(x) summary(x)$coefficients[8])
    PVEC <- as.vector(unlist(SUMLIST))
    
    BPM <- matrix(nrow = 0, ncol = nresp)
    
   tryCatch({
      
      for (b in 1:nboot) {
        
        RSDAT <- apply(dati[seq(from = 1, to = nresp, by = 1)], 2, sample, replace = TRUE)
        RFDAT <- size-RSDAT
        colnames(RSDAT) <- paste("s",1:nresp, sep="")
        colnames(RFDAT) <- paste("f",1:nresp, sep="")
        
        RSDF <- cbind(as.data.frame(RSDAT), as.data.frame(RFDAT), "dose"=as.numeric(as.character(ftrt)))
        
        RSFITLIST <- FORMLIST
        for(k in 1:NRESP){RSFITLIST[[k]] <- glm(FORMLIST[[k]], data=RSDF, family=binomial)}
        
        RSSUMLIST <- lapply(RSFITLIST, function(x) summary(x)$coefficients[8])
        
        RSPVEC <- as.vector(unlist(RSSUMLIST))
        
        BPM <- rbind(BPM,as.integer(min(RSPVEC) <= PVEC))
        
      }
      
      nrej <- nrej + sum(na.omit((colSums(BPM)/nboot)<alpha))
      nreji <- nreji + as.integer((any(na.omit((colSums(BPM)/nboot)<alpha))))
      nfail <- nfail + sum(is.na(colSums(BPM)/nboot))
      
      
      
   },
  error = function(e){
   nerror <<- nerror + 1
  })
  }
  
  return(data.frame(
    method = "bt",
    model = "glm",
    nsim = nsim, pis = pis, nresp = nresp, 
    ntrt = ntrt, corr = corr, 
    nfail = nfail, nerror = nerror, 
    nreji = nreji, nrej = nrej))
}

system.time(sim_mmm_bin_res <- do.call(rbind,apply(as.matrix(simdat_reg_FWER), 1, 
                                                   function(x){sim_mmm_bin(nsim=unname(x[1]),pis=unname(x[2]),nresp=unname(x[3]), 
                                                                           ntrt=unname(x[4]),size=unname(x[5]),corr=unname(x[6]), 
                                                                           nboot=unname(x[7]))})))

write.csv(sim_mmm_bin_res, ".\\intermediate_results\\bin_reg_FWER_bt_glm_1.csv", row.names = FALSE)

