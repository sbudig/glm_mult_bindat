unlink(".RData")
unlink(".Rhistory")

#set appropriate wokring directory
#setwd(".\Code_and_Data\Simulation_results\Park_Generation_of_correlated_binary_data")

# required packages -------------------------------------------------------

#install.packages("multcomp")

library("multcomp")
library("tidyverse")

# required helper functions -----------------------------------------------

## logit function
logit <- function(p) {
  log(p / (1 - p))
}

## backtransform logit in p(event)
expit <- function(x) {
  exp(x) / (1 + exp(x))
}

## fill matrix with specific dimension and mi
mi0ed <- function(mi = 10,
                  ndos = 4,
                  nresp = 3) {
  mimat <- matrix(rep(mi, times = ndos * nresp), ncol = ndos)
  return(mimat)
}

## helper function for parameter input;
## build a list of binomial distribution parameters
paralistbin <- function(size = 1, prob) {
  sv <- rep(size, length.out = length(prob))
  outlist <- list()
  for (i in 1:length(prob)) {
    outlist[[i]] <- list(size = sv[i], prob = prob[i])
  }
  return(outlist)
}

# paralistbin(size=10, prob=c(0.6,0.2,0.2))

# Generate correlated binary data function --------------------------------

simind_corbin_prv <- function(ntrt, mutrt, covmat, size) {
  noend <- ncol(covmat)
  ltrt <- names(ntrt)
  nvectrt <- as.vector(ntrt)
  ftrt <- factor(rep(ltrt, times = nvectrt), levels = ltrt)
  
  probtrt <- expit(mutrt)
  
  DAT <- NULL
  
  for (i in 1:length(ntrt)) {
    DATi <-
      draw_correlated_binary(
        no.row = ntrt[i],
        d = noend,
        prop.vec = probtrt[i, ],
        corr.mat = covmat
      )
    DAT <- rbind(DAT, DATi)
    
  }
  
  SDAT <- DAT
  FDAT <- size - DAT
  colnames(SDAT) <- paste("s", 1:noend, sep = "")
  colnames(FDAT) <- paste("f", 1:noend, sep = "")
  
  OUTDAT <-
    cbind(
      as.data.frame(SDAT),
      as.data.frame(FDAT),
      n = size,
      "dose" = as.numeric(as.character(ftrt))
    )
  return(OUTDAT)
}

# obtain correlated binary data

draw_correlated_binary <- function (no.row, d, prop.vec, corr.mat)
{
  alpha = matrix(0, d, d)
  
  for (i in 1:d) {
    for (j in 1:d) {
      alpha[i, j] = log(1 + corr.mat[i, j] * sqrt((1 -
                                                     prop.vec[i]) * (1 - prop.vec[j]) /
                                                    (prop.vec[i] *
                                                       prop.vec[j])))
    }
  }
  
  if (d == 2) {
    x = matrix(0, no.row, d)
    y = matrix(0, no.row, d)
    pois = numeric(d)
    sump = numeric(d)
    for (k in 1:no.row) {
      pois[1] = rpois(1, alpha[1, 1] - alpha[1, 2])
      pois[2] = rpois(1, alpha[2, 2] - alpha[1, 2])
      pois[3] = rpois(1, alpha[1, 2])
      
      sump[1] = pois[1] + pois[3]
      sump[2] = pois[2] + pois[3]
      
      
      x[k,] = sump
    }
    y[x == 0] = 1
    y[x != 0] = 0
    y
    
  } else {
    beta = matrix(0, d, d * d)
    summ = 1
    counter = 0
    while (summ > 0) {
      counter = counter + 1
      minloc = loc.min(alpha, d)
      w = matrix(1, d, d)
      mat.min = apply(alpha, 2, min)
      pos = c(1:d)
      zero.pos = which(mat.min == 0)
      if (length(zero.pos) == 0) {
        nonzero.pos = pos
      } else {
        nonzero.pos = pos[-zero.pos]
      }
      
      my.min = apply(matrix(alpha[,-minloc], d, d - length(unique(minloc))),
                     2, min)
      if (length(my.min) == 1) {
        w[,-minloc][my.min == 0] = 0
        w[-minloc,][my.min == 0] = 0
      }
      if (length(my.min) > 1) {
        w[,-minloc][, my.min == 0] = 0
        w[-minloc,][my.min == 0,] = 0
        w[alpha == 0] = 0
      }
      for (i in 1:d) {
        # if ( sum( mat.min != 0) =  d ) {
        #   T = rep(1,d)
        # }
        beta[i, counter] = alpha[minloc[1], minloc[2]] *
          1 * ((minloc[1] == i) |
                 (minloc[2] == i) | (i %in% nonzero.pos))
        
      }
      
      if (0 %in% diag(alpha)[minloc]) {
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
        pois[j] = rpois(1, max(tbeta[j,]))
      }
      for (i in 1:d) {
        sump[i] = sum(pois * w[, i])
      }
      x[k,] = sump
    }
    y[x == 0] = 1
    y[x != 0] = 0
    y
  }
}

loc.min <- function(my.mat, d) {
  w = is.matrix(my.mat)
  if (w == F) {
    stop("This is not a matrix!\n")
  }
  if (nrow(my.mat) != ncol(my.mat)) {
    stop("This is not a square matrix!\n")
  }
  n = nrow(my.mat)
  my.vec = as.vector(t(my.mat))
  my.vec[my.vec == 0] = 999
  my.index = min((1:length(my.vec))[my.vec == min(my.vec)])
  row.index = floor((my.index - 1) / n) + 1
  col.index = my.index - d * floor((my.index - 1) / n)
  c(row.index, col.index)
}

# Simulation Function -----------------------------------------------------

sim_mmm_bin <-
  function(nsim,
           di = c(0, 2, 5, 10),
           ntrt,
           pi0,
           size,
           alpha = 0.05) {
    NTRT <- rep(ntrt, times = length(di))
    names(NTRT) <- di
    
    mutrt <- mi0ed(mi = logit(pi0),
                   ndos = length(di),
                   nresp = 3)
    MUTRT <- t(mutrt)
    
    NRESP <- nrow(mutrt)
    NAMS <- paste("s", 1:NRESP, sep = "")
    NAMF <- paste("f", 1:NRESP, sep = "")
    FORMVEC <-
      paste("cbind(", NAMS, ",",  NAMF, ") ~ dose", sep = "")
    
    FORMLIST <- as.list(FORMVEC)
    names(FORMLIST) <- NAMS
    
    VCOV <- matrix(c(1, 0.9, 0.5,
                     0.9,   1, 0.5,
                     0.5, 0.5,   1),
                   nrow = 3,
                   ncol = 3)
    
    print("----------------------------")
    print(paste0("ntrt", ntrt))
    print(paste0("pi0", pi0))
  
    nerror <- 0
    
    results <-
      data.frame(
        cor_s_s1s2 = numeric(),
        cor_s_s1s3 = numeric(),
        cor_s_s2s3 = numeric(),
        cor_m_s1s2 = numeric(),
        cor_m_s1s3 = numeric(),
        cor_m_s2s3 = numeric(),
        cor_p_s1s2 = numeric(),
        cor_p_s1s3 = numeric(),
        cor_p_s2s3 = numeric(),
        t1 = numeric(),
        t2 = numeric(),
        t3 = numeric(),
        ntrt = numeric(),
        pi0 = numeric(),
        nsim = numeric()
      )
    
    for (i in 1:nsim)
    {
      dati <-
        simind_corbin_prv(
          ntrt = NTRT,
          mutrt = MUTRT,
          covmat = VCOV,
          size = size
        )
      
      tryCatch({
        FITLIST <- FORMLIST
        for (k in 1:NRESP) {
          FITLIST[[k]] <- glm(FORMLIST[[k]], data = dati, family = binomial)
        }
        
        class(FITLIST) <- "mmm"
        
        
        XT <- seq(from = 2,
                  to = 2 * NRESP,
                  by = 2)
        K <- diag(length(coef(FITLIST)))[XT,]
        rownames(K) <- names(coef(FITLIST))[XT]
        
        GF <- glht(FITLIST, linfct = K)
        
        covm <- vcov(GF)
        cor_mmm <- cov2cor(covm)
        
        cor_pear <- cor(dati[, c(1, 2, 3)])
        
        betahat <- coef(GF)
        m <- coef(GF, rhs = TRUE)
        
        ses <- sqrt(diag(covm))
        tstat <- (betahat - m) / ses
        
        
        temp <-
          data.frame(
            cor_s_s1s2 =  VCOV[1, 2],
            cor_s_s1s3 =  VCOV[1, 3],
            cor_s_s2s3 =  VCOV[2, 3],
            cor_m_s1s2 = cor_mmm[1, 2],
            cor_m_s1s3 = cor_mmm[1, 3],
            cor_m_s2s3 = cor_mmm[2, 3],
            cor_p_s1s2 = cor_pear[1, 2],
            cor_p_s1s3 = cor_pear[1, 3],
            cor_p_s2s3 = cor_pear[2, 3],
            t1 =  as.numeric(tstat[1]),
            t2 = as.numeric(tstat[2]),
            t3 = as.numeric(tstat[3]),
            ntrt = ntrt,
            pi0 = pi0,
            nsim = nsim
          )
        
        results <- rbind(results, temp)
        
        
      },
      error = function(e) {
        nerror <<- nerror + 1
      })
    }
    
    return(results)
  }


simdat <- expand.grid(
  nsim = 5000,
  ntrt = c(25, 100),
  pi0 = c(0.01),
  size = 1
)

system.time(sim_mmm_bin_res <-
              do.call(rbind, apply(as.matrix(simdat), 1,
                                   function(x) {
                                     sim_mmm_bin(
                                       nsim = unname(x[1]),
                                       ntrt = unname(x[2]),
                                       pi0 = unname(x[3]),
                                       size =
                                         unname(x[4])
                                     )
                                   })))

write.csv(sim_mmm_bin_res, "bin_cor_endpoints.csv", row.names = FALSE)



