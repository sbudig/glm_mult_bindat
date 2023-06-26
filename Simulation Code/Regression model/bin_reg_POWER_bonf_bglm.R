
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
# paralistbin(size=10, prob=c(0.6,0.2,0.2))

## define an correlated matrix
equicorrvc <- function(sigma, corr){
  nresp <- length(sigma)
  ECM <- matrix(rep(corr, times=nresp^2),ncol=nresp)
  diag(ECM) <- rep(1,nresp)
  VCM <- diag(sigma)%*%ECM%*%diag(sigma) 
  return(VCM)
}


# Generate correlated binary data function --------------------------------
simind_corbin_prv <- function(ntrt, mutrt, covmat, size, hyp, neph0){
  
  noend <- ncol(covmat) 
  ltrt <- names(ntrt)
  nvectrt <- as.vector(ntrt)
  ftrt <- factor(rep(ltrt, times=nvectrt), levels=ltrt)
  
  probtrt <- expit(mutrt)
  
  DAT <- NULL
  
  for(i in 1:length(ntrt)){
    
    DATi <- draw_correlated_binary(no.row=ntrt[i], d=noend, prop.vec=probtrt[i,], corr.mat=covmat)
    #DATi <- draw.correlated.binary(no.row=ntrt[i], d=noend, prop.vec=probtrt[i,], corr.mat=covmat)
    DAT <- rbind(DAT, DATi)
    
  }
  
  SDAT <- DAT
  FDAT <- size-DAT
  if(hyp == "h0"){
    colnames(SDAT) <- paste("s",1:noend, sep="")
    colnames(FDAT) <- paste("f",1:noend, sep="")
  } else if(hyp == "ha"){
    colnames(SDAT) <- paste("s",(neph0+1):(neph0+noend), sep="")
    colnames(FDAT) <- paste("f",(neph0+1):(neph0+noend), sep="")  
  }
  
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

sim_mmm_bin <- function(nsim, ntrt, pi0, nep_h0, nep_ha, slope, size, corr, 
                        di = c(0, 2.5, 5, 10), alpha = 0.05){
  
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
    
    dat_h0 <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_h0, covmat=VCOV_h0, size=size, hyp = "h0", neph0 = nep_h0)
    dat_ha <- simind_corbin_prv(ntrt=NTRT, mutrt=MUTRT_ha, covmat=VCOV_ha, size=size, hyp = "ha", neph0 = nep_h0)
    
    dati <- cbind(dat_h0,dat_ha)
    #dplyr::select(dati, starts_with("s"))
    
    
    FITLIST <- FORMLIST
    for(k in 1:NRESP){FITLIST[[k]] <- bayesglm(FORMLIST[[k]], data=dati, family=binomial)}
    
    tryCatch({
    
      SUMLIST <- lapply(FITLIST, function(x) summary(x)$coefficients[8])
      PVEC <- as.vector(unlist(SUMLIST))
      
      pwd <- p.adjust(PVEC, method = "bonferroni" )
      
      nrej <- nrej + sum(na.omit(pwd[seq_h0]<alpha))
      nreji <- nreji + as.integer(any(na.omit(pwd[seq_h0]<alpha)))
      nrejcor <- nrejcor + as.integer((sum(na.omit(pwd[seq_h0]<alpha)) == 0 ))
      
      npow <- npow + sum(na.omit(pwd[seq_ha]<alpha))
      npowi <- npowi + as.integer(any(na.omit(pwd[seq_ha]<alpha)))
      npowcor <- npowcor + as.integer((sum(na.omit(pwd[seq_ha]<alpha)) == nep_ha))
      
      nfail <- nfail + sum(is.na(pwd))
    }, 
    error = function(e){
      nerror <<- nerror + 1
    })
  }
  
  return(data.frame(
    method = "bonf",
    model = "glm",
    nsim = nsim, ntrt = ntrt, pi0 = pi0,
    nep_h0 = nep_h0, nep_ha = nep_ha, 
    slope = slope, corr = corr, 
    nfail = nfail, nerror = nerror, 
    nrej = nrej, nreji = nreji, nrejcor = nrejcor, 
    npow = npow, npowi = npowi, npowcor = npowcor))
}

simdat <- expand.grid(
  nsim = 10000,
  ntrt = c(25, 50, 100),
  pi0 = 0.1,
  nep_h0 = c(8,2),
  nep_ha=c(2,8),
  slope = c(0.001, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2),
  size = 1,
  corr = c(0, 0.5, 0.9)
)

simdat <- filter(simdat, (nep_h0+nep_ha)==10)

system.time(sim_mmm_bin_res <- do.call(rbind,apply(as.matrix(simdat), 1, 
                                                   function(x){sim_mmm_bin(nsim=unname(x[1]),ntrt=unname(x[2]),pi0=unname(x[3]), 
                                                                           nep_h0=unname(x[4]),nep_ha=unname(x[5]),slope=unname(x[6]),
                                                                           size=unname(x[7]),corr=unname(x[8]))})))

write.csv(sim_mmm_bin_res, "bin_reg_POWER_bonf_bglm_1.csv", row.names = FALSE)


