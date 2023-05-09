
# In a bioassay by the national toxicology program (NTP, 2000), the toxicity of methyleugenol
# was studied using female mice. Methyleugonol is used as an additive in many food and
# cosmetic products. Groups of 50 mice were given methyleugenol in 0.5% methylcellulose
# at doses of 0, 37, 75 or 150 mg/kg for 105 weeks and tumor incidence was assessed at a
# total of 89 tumor sites. We restricted the analysis to the 10 tumor classifications where a
# tumor occurred more than 5 times in total across all mice. The presence or absence for the
# 10 selected tumor classifications at the 4 dose levels is shown graphically in black/white
# in figure 0.1. It is apparent that the tumor classifications t26, t27 and t28, t29 negatively
# correlate with each other, which is due to the fact that they are both subclasses of the
# same turmor.


# Toxicology and carcinogenesis studies of methyleugenol (cas no. 93-15-2) in
# f344/n rats and b6c3f1 mice (gavage studies). National Toxicology Program technical
# report series, 491:1-412

# Setup -------------------------------------------------------------------

#setwd("")

## read NTP bioassay on methyleugenol data
load(file = "miceF.rda")

## prepare data for analysis
## group variable changed into dose variable with respective dose levels
miceF$dose <- miceF$group
miceF$dose[miceF$dose==0] <-0
miceF$dose[miceF$dose==1] <-37
miceF$dose[miceF$dose==2] <-75
miceF$dose[miceF$dose==3] <-150

## tumor classifications with more than 5 occurences are selected
cst <- colSums(miceF[, 4:92])
tt5 <-as.data.frame(miceF[, names(cst[cst>5])])
tt5$id <- seq(1, 200, 1)
tt5$dose <- miceF$dose

# Generalized linear models -----------------------------------------------

## generalized linear models are fitted to each of the tumor class 
t24_glm <- glm(cbind(t24,1-t24)~dose, data = tt5, family = binomial())
t26_glm <- glm(cbind(t26,1-t26)~dose, data = tt5, family = binomial())
t27_glm <- glm(cbind(t27,1-t27)~dose, data = tt5, family = binomial())
t28_glm <- glm(cbind(t28,1-t28)~dose, data = tt5, family = binomial())
t29_glm <- glm(cbind(t29,1-t29)~dose, data = tt5, family = binomial())
t36_glm <- glm(cbind(t36,1-t36)~dose, data = tt5, family = binomial())
t41_glm <- glm(cbind(t41,1-t41)~dose, data = tt5, family = binomial())
t42_glm <- glm(cbind(t42,1-t42)~dose, data = tt5, family = binomial())
t59_glm <- glm(cbind(t59,1-t59)~dose, data = tt5, family = binomial())
t71_glm <- glm(cbind(t71,1-t71)~dose, data = tt5, family = binomial())

## glms are combined into list for further analysis
fitlist_glm <- list(t24_glm, t26_glm, t27_glm, t28_glm, t29_glm, 
                      t36_glm, t41_glm, t42_glm, t59_glm, t71_glm)


# Bayesian generalized linear models --------------------------------------

## read "arm" library for bayesian glm model fitting
#install.packages("arm")
library(arm)

## Bayesian generalized linear models are fitted to each of the tumor class 
t24_bglm <- glm(cbind(t24,1-t24)~dose, data = tt5, family = binomial())
t26_bglm <- glm(cbind(t26,1-t26)~dose, data = tt5, family = binomial())
t27_bglm <- glm(cbind(t27,1-t27)~dose, data = tt5, family = binomial())
t28_bglm <- glm(cbind(t28,1-t28)~dose, data = tt5, family = binomial())
t29_bglm <- glm(cbind(t29,1-t29)~dose, data = tt5, family = binomial())
t36_bglm <- glm(cbind(t36,1-t36)~dose, data = tt5, family = binomial())
t41_bglm <- glm(cbind(t41,1-t41)~dose, data = tt5, family = binomial())
t42_bglm <- glm(cbind(t42,1-t42)~dose, data = tt5, family = binomial())
t59_bglm <- glm(cbind(t59,1-t59)~dose, data = tt5, family = binomial())
t71_bglm <- glm(cbind(t71,1-t71)~dose, data = tt5, family = binomial())

## Bglms are combined into list for further analysis
fitlist_bglm <- list(t24_bglm, t26_bglm, t27_bglm, t28_bglm, t29_bglm, 
                      t36_bglm, t41_bglm, t42_bglm, t59_bglm, t71_bglm)


# Bonferroni adjustment (GLM) ---------------------------------------------

## p-values are extracted from summary 
coeflist_glm <- lapply(lapply(fitlist_glm, summary), coef)
coeflist_glm <- do.call(rbind.data.frame,(coeflist_glm))[,4]

## seq is used to just extract slope p-values
seq <- seq(from = 2, to = length(coeflist_glm), by = 2)

## adjusted p-values are obtained with p.adjust function 'method = "bonferroni"'
p_adj_glm_bonf <- p.adjust(coeflist_glm[seq], method = "bonferroni")
p_adj_glm_bonf


# Bonferroni adjustment (BGLM) --------------------------------------------

## p-values are extracted from summary
coeflist_bglm <- lapply(lapply(fitlist_bglm, summary), coef)
coeflist_bglm <- do.call(rbind.data.frame,(coeflist_bglm))[,4]

## seq is used to just extract slope p-values
seq <- seq(from = 2, to = length(coeflist_bglm), by = 2)

## adjusted p-values are obtained with p.adjust function 'method = "bonferroni"'
p_adj_bglm_bonf <- p.adjust(coeflist_bglm[seq], method = "bonferroni")
p_adj_bglm_bonf


# MMM adjustment (GLM) ----------------------------------------------------

## read "multcomp" library for application of 'Multiple marginal models'
#install.packages("multcomp")
library(multcomp)

## combine respective models into list of class 'mmm'
mmm_glm <- mmm(t24_glm, t26_glm, t27_glm, t28_glm, t29_glm, 
                t36_glm, t41_glm, t42_glm, t59_glm, t71_glm)

## specify the matrix with the hypotheses to be tested
seq_mmm <- seq(from=2, to= 2*length(mmm_glm), by = 2)
K <- diag(length(coef(mmm_glm)))[seq_mmm,]
rownames(K) <- names(coef(mmm_glm))[seq_mmm]

## use function 'glht' with the matrix and 'summary' to obtain the adjusted p-values
mmm_glht_glm <- glht(mmm_glm, linfct = K)
p_adj_glm_mmm <- summary(mmm_glht_glm)
p_adj_glm_mmm


# MMM adjustment (BGLM) ---------------------------------------------------

## read "multcomp" library for application of 'Multiple marginal models'
#install.packages("multcomp")
library(multcomp)

## combine respective models into list of class 'mmm'
mmm_bglm <- mmm(t24_bglm, t26_bglm, t27_bglm, t28_bglm, t29_bglm, 
               t36_bglm, t41_bglm, t42_bglm, t59_bglm, t71_bglm)

## specify the matrix with the hypotheses to be tested
seq_mmm <- seq(from=2, to= 2*length(mmm_bglm), by = 2)
K <- diag(length(coef(mmm_bglm)))[seq_mmm,]
rownames(K) <- names(coef(mmm_bglm))[seq_mmm]

## use function 'glht' with the matrix and 'summary' to obtain the adjusted p-values
mmm_glht_bglm <- glht(mmm_bglm, linfct = K)
p_adj_bglm_mmm <- summary(mmm_glht_bglm)
p_adj_bglm_mmm


# Bootstrap adjustment (GLM) ----------------------------------------------

## extract the p-values of interest with the summary function and store into vector
sumlist_bt <- lapply(fitlist_glm, function(x) summary(x)$coefficients[8])
p_vec <- as.vector(unlist(sumlist_bt))

## J is the number of endpoints
J <- length(p_vec)

## specify the matrix where bootstrap results are stored and set number of bootstraps
BPM <- matrix(nrow = 0, ncol = length(p_vec))
nboot <- 10000

for (b in 1:nboot) {
  
  ## resample data for each of the tumor vectors
  RSDAT <- apply(tt5[seq(from = 1, to = J, by = 1)], 2, sample)
  RFDAT <- size-RSDAT
  
  ## create new data.frame with the resampled data and success/occurence and failure/non-occurence
  colnames(RSDAT) <- paste("s",1:J, sep="")
  colnames(RFDAT) <- paste("f",1:J, sep="")
  RSDF <- cbind(as.data.frame(RSDAT), as.data.frame(RFDAT), "dose"=tt5$dose)
  
  ## create appropriate lists, so glms can be fitted automatically to all endpoints
  NAMS <- paste("s",1:J, sep="")
  NAMF <- paste("f",1:J, sep="")
  FORMVEC <- paste("cbind(", NAMS, ",",  NAMF, ") ~ dose", sep="")
  
  FORMLIST <- as.list(FORMVEC)
  names(FORMLIST) <- NAMS
  RSFITLIST <- FORMLIST
  
  ## fit glms for the respective endpoints
  for(k in 1:J){RSFITLIST[[k]] <- glm(FORMLIST[[k]], data=RSDF, family=binomial)}
  
  ## extract p-values of interest and store into vector
  RSSUMLIST <- lapply(RSFITLIST, function(x) summary(x)$coefficients[8])
  RSPVEC <- as.vector(unlist(RSSUMLIST))
  
  ## extract the minimum p-value and check whether it is smaller 
  ## than the respective original unadjusted p-values.  
  ## Note the result for the respective endpoint
  BPM <- rbind(BPM,as.integer(min(RSPVEC) <= p_vec))
  
}

## calculate how many times the minimum p value from the respective bootstrap run was smaller 
## than the respective unadjusted p-value

p_adj_glm_bt <- colSums(BPM)/10000
p_adj_glm_bt


# Bootstrap adjustment (BGLM) ---------------------------------------------

## extract the p-values of interest with the summary function and store into vector
sumlist_bt <- lapply(fitlist_bglm, function(x) summary(x)$coefficients[8])
p_vec <- as.vector(unlist(sumlist_bt))

## J is the number of endpoints
J <- length(p_vec)

## specify the matrix where bootstrap results are stored and set number of bootstraps
BPM <- matrix(nrow = 0, ncol = length(p_vec))
nboot <- 10000

for (b in 1:nboot) {
  
  ## resample data for each of the tumor vectors
  RSDAT <- apply(tt5[seq(from = 1, to = J, by = 1)], 2, sample)
  RFDAT <- size-RSDAT
  
  ## create new data.frame with the resampled data and success/occurence and failure/non-occurence
  colnames(RSDAT) <- paste("s",1:J, sep="")
  colnames(RFDAT) <- paste("f",1:J, sep="")
  RSDF <- cbind(as.data.frame(RSDAT), as.data.frame(RFDAT), "dose"=tt5$dose)
  
  ## create appropriate lists, so glms can be fitted automatically to all endpoints
  NAMS <- paste("s",1:J, sep="")
  NAMF <- paste("f",1:J, sep="")
  FORMVEC <- paste("cbind(", NAMS, ",",  NAMF, ") ~ dose", sep="")
  
  FORMLIST <- as.list(FORMVEC)
  names(FORMLIST) <- NAMS
  RSFITLIST <- FORMLIST
  
  ## fit glms for the respective endpoints
  for(k in 1:J){RSFITLIST[[k]] <- glm(FORMLIST[[k]], data=RSDF, family=binomial)}
  
  ## extract p-values of interest and store into vector
  RSSUMLIST <- lapply(RSFITLIST, function(x) summary(x)$coefficients[8])
  RSPVEC <- as.vector(unlist(RSSUMLIST))
  
  ## extract the minimum p-value and check whether it is smaller 
  ## than the respective original unadjusted p-values.  
  ## Note the result for the respective endpoint
  BPM <- rbind(BPM,as.integer(min(RSPVEC) <= p_vec))
  
}

## calculate how many times the minimum p value from the respective bootstrap run was smaller 
## than the respective unadjusted p-value

p_adj_bglm_bt <- colSums(BPM)/10000
p_adj_bglm_bt





