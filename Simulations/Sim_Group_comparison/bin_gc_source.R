
# Parameter Settings ------------------------------------------------------

# nsim: number of simulation runs (Manuskript nsim = 10000)
nsim <- 10000
# nboot: number of bootstrap runs (Manuskript nboot = 1000)
nboot <- 1000

# obtain grid for simulation parameters
# method: 1=Bonferrroni, 2=MMM, 3=Bootstrap
# model: 1=GLM, 2=BayesGLM
# ntrt: sample size per dose/treatment group
# corr: correlation between endpoints

simdat <- expand.grid(
  method = c(1,2,3),
  model = c(1,2),
  nsim = nsim,
  nboot = nboot,
  ntrt = c(25, 50),
  corr = c(0, 0.5, 0.9)
)


# Functions ---------------------------------------------------------------

# Helper functions -----------------------------------------------

## logit function
# p is probability(event)
logit <- function(p){
  log(p/(1-p))
}

## backtransform logit in p(event)
expit <- function(x){
  exp(x)/(1 + exp(x))
}

## fill matrix with specific dimension and mi
# ndos: ncolumns
# nresp: nrows
# value to fill each cell
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
# 1 on the diagonal, correlation on offdiagonal
equicorrvc <- function(sigma, corr){
  nresp <- length(sigma)
  ECM <- matrix(rep(corr, times=nresp^2),ncol=nresp)
  diag(ECM) <- rep(1,nresp)
  VCM <- diag(sigma)%*%ECM%*%diag(sigma) 
  return(VCM)
}

# helper function for correlated bin
# Converts the matrix my.mat to a vector, replacing any zeros with 999.
# Finds the indices of the minimum value in this transformed vector.
# Converts these vector indices back to row and column indices of the original matrix.
# Returns the row and column indices as a vector.
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

# generate correlated binary data with respective parameters for group comparison models
# ntrt: sample size of group
# mutrt: probability of group on logitscale
# covmat: correlation matrix of endpoints
# size: cluster size (binary data = 1)
# stage: defines for which Endpoint
# NRESP: number of endpoints
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

# internal function of simind_corbin_prv to obtain respective binary data for each group
# no.row: for which group
# d: number of endpoints
# prop.vec: probability vector of respective group
# corr.mat: correlation matrix of endpoints

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

