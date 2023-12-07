# Code for simulation of FWER and Power of correlated binary data

# Setup -------------------------------------------------------------------

unlink(".RData")
unlink(".Rhistory")

# required packages -------------------------------------------------------

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

## define an correlated matrix
equicorrvc <- function(sigma, corr){
  nresp <- length(sigma)
  ECM <- matrix(rep(corr, times=nresp^2),ncol=nresp)
  diag(ECM) <- rep(1,nresp)
  VCM <- diag(sigma)%*%ECM%*%diag(sigma) 
  return(VCM)
}

# Generate correlated binary data function --------------------------------

simind_corbin_prv <- function(ntrt, mutrt, covmat, size){
  
  noend <- ncol(covmat) 
  ltrt <- names(ntrt)
  nvectrt <- as.vector(ntrt)
  ftrt <- factor(rep(ltrt, times=nvectrt), levels=ltrt)
  
  probtrt <- expit(mutrt)
  
  DAT <- NULL
  
  for(i in 1:length(ntrt)){
    
    DATi <- draw_correlated_binary(no.row=ntrt[i], d=noend, prop.vec=probtrt[i,], corr.mat=covmat)
    DAT <- rbind(DAT, DATi)
    
  }
  
  SDAT <- DAT
  FDAT <- size-DAT
  colnames(SDAT) <- paste("s",1:noend, sep="")
  colnames(FDAT) <- paste("f",1:noend, sep="")
  
  OUTDAT <- cbind(as.data.frame(SDAT), as.data.frame(FDAT), n=size, "dose"=as.numeric(as.character(ftrt)))
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

# helper function for correlated bin
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



# Simulation Function -----------------------------------------------------

sim_mmm_bin <- function(nsim, ntrt, pis, nresp, size, corr, nboot, di = c(0, 2.5, 5, 10), alpha = 0.05){
  
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
    
    dati <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT, covmat=VCOV, size=size)
    
    
    FITLIST <- FORMLIST
    for(k in 1:NRESP){FITLIST[[k]] <- bayesglm(FORMLIST[[k]], data=dati, family=binomial)}
    
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
    model = "bglm",
    nsim = nsim, pis = pis, nresp = nresp, 
    ntrt = ntrt, corr = corr, 
    nfail = nfail, nerror = nerror, 
    nreji = nreji, nrej = nrej))
}


simdat <- expand.grid(
  nsim = 10000,
  pis = c(0.001, 0.0625, 0.125, 0.25, 0.375, 0.5),
  nresp = c(3, 5, 10),
  ntrt = c(25, 50, 100),
  size = 1,
  corr = c(0, 0.5, 0.9),
  nboot = 10000
)

system.time(sim_mmm_bin_res <- do.call(rbind,apply(as.matrix(simdat), 1, 
                                                   function(x){sim_mmm_bin(nsim=unname(x[1]),pis=unname(x[2]),nresp=unname(x[3]), 
                                                                           ntrt=unname(x[4]),size=unname(x[5]),corr=unname(x[6]), 
                                                                           nboot=unname(x[7]))})))

write.csv(sim_mmm_bin_res, "bin_reg_FWER_bt_bglm_1.csv", row.names = FALSE)



