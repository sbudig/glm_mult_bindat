
unlink(".RData")
unlink(".Rhistory")

args <- commandArgs(trailingOnly = TRUE)
# required packages -------------------------------------------------------

library("multcomp")
library("tidyverse")
library("arm")

# required helper functions -----------------------------------------------

## logit function
logit <- function(p){
  log(p/(1-p))
}

## backtransform logit in p(event)
expit <- function(x){
  exp(x)/(1 + exp(x))
}

## fill matrix with specific dimension and mi
mi0ed <- function(mi=10, ndos=4, nresp=3){
  mimat <- matrix(rep(mi, times=ndos*nresp), ncol=ndos)
  return(mimat)
}

## helper function for parameter input; 
## build a list of binomial distribution parameters
paralistbin <- function(size=1, prob){
  sv <- rep(size, length.out=length(prob))
  outlist <- list()
  for(i in 1:length(prob)){outlist[[i]] <- list(size=sv[i], prob=prob[i])}
  return(outlist)
} 
# paralistbin(size=10, prob=c(0.6,0.2,0.2))

## define an correlated matrix
equicorrvc <- function(sigma, corr){
  nresp <- length(sigma)
  ECM <- matrix(rep(corr, times=nresp^2),ncol=nresp)
  diag(ECM) <- rep(1,nresp)
  VCM <- diag(sigma)%*%ECM%*%diag(sigma) 
  return(VCM)
}

# helper function for correlated binary data
loc.min <- function(my.mat,d){
  w=is.matrix(my.mat)
  if (w==F){
    stop("This is not a matrix!\n")
  }
  if (nrow(my.mat)!=ncol(my.mat)){
    stop("This is not a square matrix!\n")
  }
  n=nrow(my.mat)
  my.vec=as.vector(t(my.mat))
  my.vec[my.vec==0]=999
  my.index=min((1:length(my.vec))[my.vec==min(my.vec)])
  row.index=floor((my.index-1)/n)+1
  col.index=my.index-d*floor((my.index-1)/n)
  c(row.index,col.index)
}

# Generate correlated binary data function --------------------------------
simind_corbin_prv <- function(ntrt, mutrt, covmat, size, stage, NRESP){
  
  noend <- ncol(covmat) 
  
  probtrt <- expit(mutrt)
  
  DAT <- NULL
  
  for(i in 1:length(ntrt)){
    
    DATi <- draw_correlated_binary(no.row=ntrt[i], d=noend, prop.vec=probtrt[i,], corr.mat=covmat)
    #DATi <- draw.correlated.binary(no.row=ntrt[i], d=noend, prop.vec=probtrt[i,], corr.mat=covmat)
    DAT <- rbind(DAT, DATi)
    
  }
  
  SDAT <- DAT
  FDAT <- size-DAT
  colnames(SDAT) <- paste("s",stage:(stage+noend-1), sep="")
  colnames(FDAT) <- paste("f",stage:(stage+noend-1), sep="")

  
  OUTDAT <- cbind(as.data.frame(SDAT), as.data.frame(FDAT))
  return(OUTDAT)
  
}

# obtain correlated binary data

draw_correlated_binary<- function (no.row, d, prop.vec, corr.mat) 
{
  alpha = matrix(0, d, d)
  
  for (i in 1:d) {
    for (j in 1:d) {
      alpha[i, j] = log(1 + corr.mat[i, j] * sqrt((1 -
                                                     prop.vec[i]) * (1 - prop.vec[j])/(prop.vec[i] * 
                                                                                         prop.vec[j])))
    }
  }
  
  if (d==2) {
    
    x = matrix(0, no.row, d)
    y = matrix(0, no.row, d)
    pois = numeric(d)
    sump = numeric(d)
    for (k in 1:no.row) {
      
      pois[1] = rpois(1, alpha[1,1] - alpha[1,2])  
      pois[2] = rpois(1, alpha[2,2] - alpha[1,2]) 
      pois[3] = rpois(1, alpha[1,2])
      
      sump[1] = pois[1] + pois[3]
      sump[2] = pois[2] + pois[3]
      
      
      x[k, ] = sump
    }
    y[x == 0] = 1
    y[x != 0] = 0
    y
    
  } else {
    
    beta = matrix(0, d, d * d)
    summ = 1
    counter = 0
    while ( summ > 0 ) {
      counter = counter + 1
      minloc = loc.min(alpha, d)
      w = matrix(1, d, d)
      mat.min = apply(alpha, 2, min)
      pos = c(1:d)
      zero.pos = which(mat.min==0)
      if ( length(zero.pos) == 0 ) {
        nonzero.pos = pos
      } else {
        nonzero.pos = pos[-zero.pos]
      }
      
      my.min = apply(matrix(alpha[, -minloc], d, d - length(unique(minloc))), 
                     2, min)
      if (length(my.min) == 1) {
        w[, -minloc][my.min == 0] = 0
        w[-minloc, ][my.min == 0] = 0
      }
      if (length(my.min) > 1) {
        w[, -minloc][, my.min == 0] = 0
        w[-minloc, ][my.min == 0, ] = 0
        w[alpha == 0] = 0
      }
      for (i in 1:d) {
        # if ( sum( mat.min != 0) =  d ) {
        #   T = rep(1,d)
        # }
        beta[i, counter] = alpha[minloc[1], minloc[2]] * 
          1 * ((minloc[1] == i) | (minloc[2] == i) | (i %in% nonzero.pos))
        
      }
      
      if ( 0 %in% diag(alpha)[minloc] ) {
        stop("The method won't work in this parameter setting\n")
      }
      
      alpha = alpha - alpha[minloc[1], minloc[2]] * w
      summ = sum(alpha)
      
    }
    
    tbeta = t(beta)
    w = (tbeta != 0)
    x = matrix(0, no.row, d)
    y = matrix(0, no.row, d)
    pois = numeric(nrow(tbeta))
    sump = numeric(d)
    for (k in 1:no.row) {
      for (j in 1:nrow(tbeta)) {
        pois[j] = rpois(1, max(tbeta[j, ]))
      }
      for (i in 1:d) {
        sump[i] = sum(pois * w[, i])
      }
      x[k, ] = sump
    }
    y[x == 0] = 1
    y[x != 0] = 0
    y
  }
}


# Simulation Function -----------------------------------------------------

sim_mmm_bin <- function(method, model, nsim, nboot,
                          ntrt, corr, 
                          alpha=0.05, size = 1, di = c(0,2,5,10)){

  
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
  
  m1 <- matrix(rep(logit(0.4), each = length(di), times = 2), ncol = length(di))
  m2 <-  matrix(rep(c(logit(0.1),logit(0.1),logit(0.1),logit(0.2)), each = 8), ncol = length(di))
  
  MUTRT_m1 <- t(m1)
  MUTRT_m2 <- t(m2)

  NRESP <- ncol(MUTRT_m1)+ncol(MUTRT_m2)
  NAMS <- paste("s",1:NRESP, sep="")
  NAMF <- paste("f",1:NRESP, sep="")
  FORMVEC <- paste("cbind(", NAMS, ",",  NAMF, ") ~ dose", sep="")
  
  FORMLIST <- as.list(FORMVEC)
  names(FORMLIST) <- NAMS
  
  SIGMA_m1 <- rep(1, times = ncol(MUTRT_m1))
  SIGMA_m2 <- rep(1, times = ncol(MUTRT_m2))
  
  VCOV_m1 <- equicorrvc(sigma=SIGMA_m1, corr=corr)
  VCOV_m2 <- equicorrvc(sigma=SIGMA_m2, corr=corr)
  
  print(paste0(method, model, " ntrt ", ntrt, " corr ", corr))
  
  h0havec <- as.numeric(apply(cbind(MUTRT_m1,MUTRT_m2), 2, function(x) length(unique(x)) == 1))
  
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
      dat_m2 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m2, covmat=VCOV_m2, size=size, stage = ncol(MUTRT_m1) + 1, NRESP = NRESP)

      dati <- cbind(dat_m1,dat_m2)
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
      dat_m2 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m2, covmat=VCOV_m2, size=size, stage = ncol(MUTRT_m1) + 1, NRESP = NRESP)
      
      dati <- cbind(dat_m1,dat_m2)
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
      dat_m2 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_m2, covmat=VCOV_m2, size=size, stage = ncol(MUTRT_m1) + 1, NRESP = NRESP)
      
      dati <- cbind(dat_m1,dat_m2)
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

# Parameter Settings Grid

simdat <- expand.grid(
  method = c(1,2,3),
  model = c(1,2),
  nsim = 10000,
  nboot = 1000,
  ntrt = c(25, 50),
  corr = c(0, 0.5, 0.9)
  )

system.time(sim_mmm_bin_res <- do.call(rbind,apply(as.matrix(simdat), 1, 
                                                   function(x){sim_mmm_bin(method=unname(x[1]),model=unname(x[2]),
                                                                           nsim=unname(x[3]),nboot=unname(x[4]),
                                                                           ntrt=unname(x[5]),
                                                                           corr=unname(x[6]))})))
  
write.csv(sim_mmm_bin_res, paste0("bin_gc_ss_h02_04_ha8_01_d01_",args[1],".csv"), row.names = FALSE)
