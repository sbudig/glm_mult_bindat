# Setup -------------------------------------------------------------------

# takes the number from the batch file as an input to save the csv file 
# with the respective number extension
args <- commandArgs(trailingOnly = TRUE)

# Functions
source("bin_gc_source.R")

# required packages -------------------------------------------------------

library("multcomp")
library("tidyverse")
library("arm")

# Simulation Function -----------------------------------------------------
# method: 1=Bonferrroni, 2=MMM, 3=Bootstrap
# model: 1=GLM, 2=BayesGLM
# nsim: number of simulations
# nboot: number of bootstrap samples
# ntrt: sample size per dose/treatment group
# corr: correlation between endpoints
# alpha: type 1 error
# size: cluster size (=1 for binary data)
# di: dose groups
sim_mmm_bin <- function(method, 
                        model, 
                        nsim, 
                        nboot,
                        ntrt, 
                        corr,
                        alpha=0.05, 
                        size = 1, 
                        di = c(0,2,5,10)){

  
  # m2 <- matrix(c(0.075, 0.075, 0.075, 0.075,
  #                0.1,   0.2,   0.3,   0.4,
  #                0.1,   0.1,   0.1,   0.1,
  #                0.1,   0.1,   0.1,   0.1,
  #                0.15,  0.15,  0.15,  0.15,
  #                0.125, 0.125, 0.125, 0.125,
  #                0.2,   0.3,   0.4,   0.5,
  #                0.3,   0.3,   0.3,   0.3,
  #                0.6,   0.6,   0.6,   0.6,
  #                0.8,   0.8,   0.8,   0.8),
  #              nrow = 4, ncol = 10)  
  
  
  method <- case_when(method == 1 ~ "bonf",
                      method == 2 ~ "mmm",
                      method == 3 ~ "bt")
  model <- case_when(model == 1 ~ "glm",
                     model == 2 ~ "bglm")  
  
  
  NTRT <- rep(ntrt, times = length(di))
  names(NTRT) <- di
  
  m1 <- matrix(c( 0.01,   0.01,   0.01,   0.01,
                 0.01,   0.01,   0.01,   0.01),
               nrow = 4, ncol = 2)
  
  m2 <- matrix(c( 0.25,   0.25,   0.25,   0.25,
                  0.25,   0.25,   0.25,   0.25),
               nrow = 4, ncol = 2)
  
  m3 <- matrix(c( 0.5,   0.5,   0.5,   0.5,
                   0.5,   0.5,   0.5,   0.5),
                nrow = 4, ncol = 2)
  
  m4 <- matrix(c( 0.125,   0.125,   0.125,   0.125,
                  0.125,   0.125,   0.125,   0.125),
               nrow = 4, ncol = 2)
  
  m5 <- matrix(c( 0.375,   0.375,   0.375,   0.375,
                  0.375,   0.375,   0.375,   0.375),
               nrow = 4, ncol = 2)
  



  MUTRT_m1 <- logit(m1)
  MUTRT_m2 <- logit(m2)
  MUTRT_m3 <- logit(m3)
  MUTRT_m4 <- logit(m4)
  MUTRT_m5 <- logit(m5)

  NRESP <- ncol(m1)+ncol(m2)+ncol(m3)+ncol(m4)+ncol(m5)
  NAMS <- paste("s",1:NRESP, sep="")
  NAMF <- paste("f",1:NRESP, sep="")
  FORMVEC <- paste("cbind(", NAMS, ",",  NAMF, ") ~ dose", sep="")
  
  FORMLIST <- as.list(FORMVEC)
  names(FORMLIST) <- NAMS
  
  SIGMA_m1 <- rep(1, times = ncol(m1))
  SIGMA_m2 <- rep(1, times = ncol(m2))
  SIGMA_m3 <- rep(1, times = ncol(m3))
  SIGMA_m4 <- rep(1, times = ncol(m4))
  SIGMA_m5 <- rep(1, times = ncol(m5))
  
  VCOV_m1 <- equicorrvc(sigma=SIGMA_m1, corr=corr)
  VCOV_m2 <- equicorrvc(sigma=SIGMA_m2, corr=corr)
  VCOV_m3 <- equicorrvc(sigma=SIGMA_m3, corr=corr)
  VCOV_m4 <- equicorrvc(sigma=SIGMA_m4, corr=corr)
  VCOV_m5 <- equicorrvc(sigma=SIGMA_m5, corr=corr)
  
  print(paste0(method, model, " ntrt ", ntrt, " corr ", corr))
  
  h0havec <- as.numeric(apply(cbind(m1,m2,m3,m4,m5), 2, function(x) length(unique(x)) == 1))
  
  seq_h0 <- seq(from = 1, to = sum(h0havec)*3, by = 1)
  nep_ha <- sum(1-h0havec)
  if (nep_ha>0){
  seq_ha <- seq(from = max(seq_h0)+1, to = max(seq_h0)+sum(1-h0havec)*3, by = 1)
  }
  
  ltrt <- names(NTRT)
  nvectrt <- as.vector(NTRT)
  ftrt <- factor(rep(ltrt, times=nvectrt), levels=ltrt)
  
  nrej <- 0   # How many rejected per Dataset
  nreji <- 0  # At least one rejected per Dataset
  nrejcor <- 0 # how often all h0 are true
  
  npow <- 0   # How many differences found per Dataset
  npowi <- 0  # At least one difference found per Dataset
  npowcor <- 0 # How often find both differences
  
  nfail <- 0  
  nerror <- 0  
  
  nh <- seq(from = 1, to = (NRESP)*4)
  nh <- nh[-seq(1, length(nh),4)]
  
  # Bonferroni Simulation
  if (method == "bonf"){
    
    for(i in 1:nsim)
    {
      
      dat_m1 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m1, covmat=VCOV_m1, size=size, stage = 1, NRESP = NRESP)
      dat_m2 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m2, covmat=VCOV_m2, size=size, stage = ncol(m1) + 1, NRESP = NRESP)
      dat_m3 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m3, covmat=VCOV_m3, size=size, stage = ncol(m1) + ncol(m2) + 1, NRESP = NRESP)
      dat_m4 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m4, covmat=VCOV_m4, size=size, stage = ncol(m1) + ncol(m2) + ncol(m3) + 1, NRESP = NRESP)
      dat_m5 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m5, covmat=VCOV_m5, size=size, stage = ncol(m1) + ncol(m2) + ncol(m3) + ncol(m4) + 1, NRESP = NRESP)

      dati <- cbind(dat_m1,dat_m2,dat_m3,dat_m4,dat_m5)
      dati$dose <- as.factor(ftrt)
      
      FITLIST <- FORMLIST
      if (model == "glm"){
        for(k in 1:NRESP){FITLIST[[k]] <- glm(FORMLIST[[k]], data=dati, family=binomial)}
      }
      if (model == "bglm"){
        for(k in 1:NRESP){FITLIST[[k]] <- bayesglm(FORMLIST[[k]], data=dati, family=binomial)}
      }
      
      tryCatch({
        
        s_coef <- lapply(FITLIST, function(x) summary(x)$coefficients[,4])
        p_nh <- as.vector(unlist(s_coef))[nh]
        pwd <- p.adjust(p_nh, method = "bonferroni" )
        
        nrej <- nrej + sum(na.omit(pwd[seq_h0]<alpha))
        nreji <- nreji + as.integer(any(na.omit(pwd[seq_h0]<alpha)))
        nrejcor <- nrejcor + as.integer((sum(na.omit(pwd[seq_h0]<alpha)) == 0 ))
        
        if (sum(1-h0havec)>0){
        npow <- npow + sum(na.omit(pwd[seq_ha]<alpha))
        npowi <- npowi + as.integer(any(na.omit(pwd[seq_ha]<alpha)))
        npowcor <- npowcor + as.integer((sum(na.omit(pwd[seq_ha]<alpha)) == nep_ha))
        }
        nfail <- nfail + sum(is.na(pwd))
        
      }, 
      error = function(e){
        nerror <<- nerror + 1
      })
    }
  }
  
  # MMM Simulation
  if (method == "mmm"){
    
    for(i in 1:nsim)
    {
      
      dat_m1 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m1, covmat=VCOV_m1, size=size, stage = 1, NRESP = NRESP)
      dat_m2 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m2, covmat=VCOV_m2, size=size, stage = ncol(m1) + 1, NRESP = NRESP)
      dat_m3 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m3, covmat=VCOV_m3, size=size, stage = ncol(m1) + ncol(m2) + 1, NRESP = NRESP)
      dat_m4 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m4, covmat=VCOV_m4, size=size, stage = ncol(m1) + ncol(m2) + ncol(m3) + 1, NRESP = NRESP)
      dat_m5 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m5, covmat=VCOV_m5, size=size, stage = ncol(m1) + ncol(m2) + ncol(m3) + ncol(m4) + 1, NRESP = NRESP)
      
      dati <- cbind(dat_m1,dat_m2,dat_m3,dat_m4,dat_m5)
      dati$dose <- as.factor(ftrt)
      
      FITLIST <- FORMLIST
      if (model == "glm"){
        for(k in 1:NRESP){FITLIST[[k]] <- glm(FORMLIST[[k]], data=dati, family=binomial)}
      }
      if (model == "bglm"){
        for(k in 1:NRESP){FITLIST[[k]] <- bayesglm(FORMLIST[[k]], data=dati, family=binomial)}
      }
      
      tryCatch({
        
        class(FITLIST) <- "mmm"
        
        K <- diag(length(coef(FITLIST)))[nh,]
        rownames(K) <- names(coef(FITLIST))[nh]
        
        GF <- glht(FITLIST, linfct = K)
        SGF <- summary(GF)
        
        nrej <- nrej + sum(na.omit(SGF$test$pvalues[seq_h0]<alpha))
        nreji <- nreji + as.integer(any(na.omit(SGF$test$pvalues[seq_h0]<alpha)))
        nrejcor <- nrejcor + as.integer((sum(na.omit(SGF$test$pvalues[seq_h0]<alpha)) == 0 ))
        
        if (sum(1-h0havec)>0){
        npow <- npow + sum(na.omit(SGF$test$pvalues[seq_ha]<alpha))
        npowi <- npowi + as.integer(any(na.omit(SGF$test$pvalues[seq_ha]<alpha)))
        npowcor <- npowcor + as.integer((sum(na.omit(SGF$test$pvalues[seq_ha]<alpha)) == nep_ha))
        }
        
        nfail <- nfail + sum(is.na(SGF$test$pvalues))
        
      }, 
      error = function(e){
        nerror <<- nerror + 1
      })
    }
  }
  
  # Bootstrap Simulation
  if (method == "bt"){
    for(i in 1:nsim)
    {
      
      dat_m1 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m1, covmat=VCOV_m1, size=size, stage = 1, NRESP = NRESP)
      dat_m2 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m2, covmat=VCOV_m2, size=size, stage = ncol(m1) + 1, NRESP = NRESP)
      dat_m3 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m3, covmat=VCOV_m3, size=size, stage = ncol(m1) + ncol(m2) + 1, NRESP = NRESP)
      dat_m4 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m4, covmat=VCOV_m4, size=size, stage = ncol(m1) + ncol(m2) + ncol(m3) + 1, NRESP = NRESP)
      dat_m5 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m5, covmat=VCOV_m5, size=size, stage = ncol(m1) + ncol(m2) + ncol(m3) + ncol(m4) + 1, NRESP = NRESP)
      
      dati <- cbind(dat_m1,dat_m2,dat_m3,dat_m4,dat_m5)
      dati$dose <- as.factor(ftrt)
      
      FITLIST <- FORMLIST
      if (model == "glm"){
        for(k in 1:NRESP){FITLIST[[k]] <- glm(FORMLIST[[k]], data=dati, family=binomial)}
      }
      if (model == "bglm"){
        for(k in 1:NRESP){FITLIST[[k]] <- bayesglm(FORMLIST[[k]], data=dati, family=binomial)}
      }
      
      SUMLIST <- lapply(FITLIST, function(x) summary(x)$coefficients[,4])
      PVEC <- as.vector(unlist(SUMLIST))[nh]
      
      BPM <- matrix(nrow = 0, ncol = length(nh))
      
      tryCatch({
        
        for (b in 1:nboot) {
          
          RSDF <- dati[sample(1:nrow(dati),nrow(dati), replace = TRUE),]
          RSDF$dose <- as.factor(as.character(ftrt))
          
          RSFITLIST <- FORMLIST
          for(k in 1:NRESP){RSFITLIST[[k]] <- glm(FORMLIST[[k]], data=RSDF, family=binomial)}
          
          RSSUMLIST <- lapply(RSFITLIST, function(x) summary(x)$coefficients[,4])
          
          RSPVEC <- as.vector(unlist(RSSUMLIST))[nh]
          
          BPM <- rbind(BPM,as.integer(min(RSPVEC) <= PVEC))
          
        }
        
        nrej <- nrej + sum(na.omit((colSums(BPM)/nboot)[seq_h0]<alpha))
        nreji <- nreji + as.integer(any(na.omit((colSums(BPM)/nboot)[seq_h0]<alpha)))
        nrejcor <- nrejcor + as.integer((sum(na.omit((colSums(BPM)/nboot)[seq_h0]<alpha)) == 0 ))
        
        if (sum(1-h0havec)>0){
        npow <- npow + sum(na.omit((colSums(BPM)/nboot)[seq_ha]<alpha))
        npowi <- npowi + as.integer(any(na.omit((colSums(BPM)/nboot)[seq_ha]<alpha)))
        npowcor <- npowcor + as.integer((sum(na.omit((colSums(BPM)/nboot)[seq_ha]<alpha)) == nep_ha))
        }
        
        nfail <- nfail + sum(is.na((colSums(BPM)/nboot)))
        
      }, 
      error = function(e){
        nerror <<- nerror + 1
      })
    }
  }
  
  return(data.frame(
    method = method,
    model = model,
    nsim = nsim, 
    ntrt = ntrt, 
    corr = corr, 
    nfail = nfail, nerror = nerror, 
    nrej = nrej, nreji = nreji, nrejcor = nrejcor,
    npow = npow, npowi = npowi, npowcor = npowcor)
  )
}

sim_mmm_bin_res <- do.call(rbind,apply(as.matrix(simdat), 1, 
                                                   function(x){sim_mmm_bin(method=unname(x[1]),model=unname(x[2]),
                                                                           nsim=unname(x[3]),nboot=unname(x[4]),
                                                                           ntrt=unname(x[5]),
                                                                           corr=unname(x[6]))}))
  
write.csv(sim_mmm_bin_res, paste0(".\\intermediate_results\\bin_gc_ss_2x001_2x0125_2x025_2x0375_2x05_",args[1],".csv"), row.names = FALSE)
