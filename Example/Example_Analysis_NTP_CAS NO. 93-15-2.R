
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

# Short Overview:
# First the data is imported and prepared for further analysis.
# Then the data set is first analysed using regression models.
# Appropariate GLMs or BayesGLMs are fitted for each tumour and combined into a list.
# Then the p-value adjustments are performed with the respective method.
# The data set was also analysed using group comparison models.

# Setup -------------------------------------------------------------------

# Set seed for reproducibility
set.seed(12345)

## read NTP bioassay on methyleugenol data
load(file = ".\\data\\miceF.rda")

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

# Regression models -------------------------------------------------------

## Generalized linear models ----------------------------------------------

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


## Bayesian generalized linear models -------------------------------------

## read "arm" library for bayesian glm model fitting
#install.packages("arm")
library(arm)

## Bayesian generalized linear models are fitted to each of the tumor class 
t24_bglm <- bayesglm(cbind(t24,1-t24)~dose, data = tt5, family = binomial())
t26_bglm <- bayesglm(cbind(t26,1-t26)~dose, data = tt5, family = binomial())
t27_bglm <- bayesglm(cbind(t27,1-t27)~dose, data = tt5, family = binomial())
t28_bglm <- bayesglm(cbind(t28,1-t28)~dose, data = tt5, family = binomial())
t29_bglm <- bayesglm(cbind(t29,1-t29)~dose, data = tt5, family = binomial())
t36_bglm <- bayesglm(cbind(t36,1-t36)~dose, data = tt5, family = binomial())
t41_bglm <- bayesglm(cbind(t41,1-t41)~dose, data = tt5, family = binomial())
t42_bglm <- bayesglm(cbind(t42,1-t42)~dose, data = tt5, family = binomial())
t59_bglm <- bayesglm(cbind(t59,1-t59)~dose, data = tt5, family = binomial())
t71_bglm <- bayesglm(cbind(t71,1-t71)~dose, data = tt5, family = binomial())

## Bglms are combined into list for further analysis
fitlist_bglm <- list(t24_bglm, t26_bglm, t27_bglm, t28_bglm, t29_bglm, 
                      t36_bglm, t41_bglm, t42_bglm, t59_bglm, t71_bglm)


## Bonferroni adjustment (GLM) --------------------------------------------

## p-values for each model are extracted from summary 
p_vec_glm <- lapply(fitlist_glm, function(x) summary(x)$coefficients[8])

## adjusted p-values are obtained with p.adjust function 'method = "bonferroni"'
p_adj_glm_bonf <- p.adjust(p_vec_glm, method = "bonferroni")

## The sorting of the p-values corresponds to the order of the models in the list.
p_adj_glm_bonf


## Bonferroni adjustment (BGLM) -------------------------------------------

## p-values for each model are extracted from summary 
p_vec_bglm <- lapply(fitlist_bglm, function(x) summary(x)$coefficients[8])

## adjusted p-values are obtained with p.adjust function 'method = "bonferroni"'
p_adj_bglm_bonf <- p.adjust(p_vec_bglm, method = "bonferroni")

## The sorting of the p-values corresponds to the order of the models in the list.
p_adj_bglm_bonf


## MMM adjustment (GLM) ---------------------------------------------------

## read "multcomp" library for application of 'Multiple marginal models'
#install.packages("multcomp")
library(multcomp)

## combine respective models into list of class 'mmm'
mmm_glm <- mmm(t24_glm, t26_glm, t27_glm, t28_glm, t29_glm, 
                t36_glm, t41_glm, t42_glm, t59_glm, t71_glm)

## specify the matrix with the hypotheses to be tested
## we just want to test every second coefficient (only the slope, not the intercept)
seq_mmm <- seq(from=2, to= 2*length(mmm_glm), by = 2)
K <- diag(length(coef(mmm_glm)))[seq_mmm,]
rownames(K) <- names(coef(mmm_glm))[seq_mmm]

## use function 'glht' with the matrix and 'summary' to obtain the adjusted p-values
mmm_glht_glm <- glht(mmm_glm, linfct = K)
p_adj_glm_mmm <- summary(mmm_glht_glm)

p_adj_glm_mmm


## MMM adjustment (BGLM) --------------------------------------------------

## read "multcomp" library for application of 'Multiple marginal models'
#install.packages("multcomp")
library(multcomp)

## combine respective models into list of class 'mmm'
mmm_bglm <- mmm(t24_bglm, t26_bglm, t27_bglm, t28_bglm, t29_bglm, 
               t36_bglm, t41_bglm, t42_bglm, t59_bglm, t71_bglm)

## specify the matrix with the hypotheses to be tested
## we just want to test every second coefficient (only the slope, not the intercept)
seq_mmm <- seq(from=2, to= 2*length(mmm_bglm), by = 2)
K <- diag(length(coef(mmm_bglm)))[seq_mmm,]
rownames(K) <- names(coef(mmm_bglm))[seq_mmm]

## use function 'glht' with the matrix and 'summary' to obtain the adjusted p-values
mmm_glht_bglm <- glht(mmm_bglm, linfct = K)
p_adj_bglm_mmm <- summary(mmm_glht_bglm)

p_adj_bglm_mmm


## Bootstrap adjustment (GLM) ---------------------------------------------

## extract the p-values of interest with the summary function and store into vector
l_pvalues_glm <- lapply(fitlist_glm, function(x) summary(x)$coefficients[8])
v_pvalues_glm <- as.vector(unlist(l_pvalues_glm))

## J is the number of endpoints
J <- length(v_pvalues_glm)

## specify the matrix where bootstrap results are stored
BPM <- matrix(nrow = 0, ncol = length(v_pvalues_glm))

# set number of bootstrap iterations
nboot <- 10000

for (b in 1:nboot) {
  
  ## resample data for each of the tumor vectors
  RSDAT <- apply(tt5[seq(from = 1, to = J, by = 1)], 2, sample, replace = TRUE)
  RFDAT <- 1-RSDAT
  
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
  BPM <- rbind(BPM,as.integer(min(RSPVEC) <= v_pvalues_glm))
  
}

## calculate how many times the minimum p value from the respective bootstrap run was smaller 
## than the respective unadjusted p-value

p_adj_glm_bt <- colSums(BPM)/nboot

## The sorting of the p-values corresponds to the order of the models in the list.
p_adj_glm_bt


## Bootstrap adjustment (BGLM) --------------------------------------------

## extract the p-values of interest with the summary function and store into vector
l_pvalues_bglm <- lapply(fitlist_bglm, function(x) summary(x)$coefficients[8])
v_pvalues_bglm <- as.vector(unlist(l_pvalues_bglm))

## J is the number of endpoints
J <- length(v_pvalues_bglm)

## specify the matrix where bootstrap results are stored
BPM <- matrix(nrow = 0, ncol = length(v_pvalues_bglm))

# set number of bootstrap iterations
nboot <- 10000

for (b in 1:nboot) {
  
  ## resample data for each of the tumor vectors
  RSDAT <- apply(tt5[seq(from = 1, to = J, by = 1)], 2, sample, replace = TRUE)
  RFDAT <- 1-RSDAT
  
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
  for(k in 1:J){RSFITLIST[[k]] <- bayesglm(FORMLIST[[k]], data=RSDF, family=binomial)}
  
  ## extract p-values of interest and store into vector
  RSSUMLIST <- lapply(RSFITLIST, function(x) summary(x)$coefficients[8])
  RSPVEC <- as.vector(unlist(RSSUMLIST))
  
  ## extract the minimum p-value and check whether it is smaller 
  ## than the respective original unadjusted p-values.  
  ## Note the result for the respective endpoint
  BPM <- rbind(BPM,as.integer(min(RSPVEC) <= v_pvalues_bglm))
  
}

## calculate how many times the minimum p value from the respective bootstrap run was smaller 
## than the respective unadjusted p-value

p_adj_bglm_bt <- colSums(BPM)/nboot

## The sorting of the p-values corresponds to the order of the models in the list.
p_adj_bglm_bt


## Dataframe - Adjusted p-values (Table 1) --------------------------------

# Generation of the data frame containing the adjusted p values of all methods 
# for the regression models. This data frame is Table 1 in the manuscript

df_regression_pvalues <- 
  data.frame(Endpoint = colnames(tt5)[1:10],
           Bonferroni = round(p_adj_glm_bonf, 4),
           B_Bonferroni = round(p_adj_bglm_bonf, 4),
           MMM = round(p_adj_glm_mmm$test$pvalues, 4),
           B_MMM = round(p_adj_bglm_mmm$test$pvalues, 4),
           Bootstrap = p_adj_glm_bt,
           B_Bootstrap = p_adj_bglm_bt)

write.csv(df_regression_pvalues, 
          ".\\results\\Table_1_Example_Analysis_NTP_Regression_pvalues.csv", 
          row.names = FALSE)


## Dataframe - Effect Sizes and SEs --------------------------------------

# Generate data frame containing the effect sizes and standard errors 
# for all regression models

reg_glm_list <- lapply(lapply(fitlist_glm, function(x) summary(x)$coefficients[,1:2]), 
                   function(z) do.call(cbind,as.data.frame(z)))
reg_GLMtable <- cbind("Endpoint" = rep(colnames(tt5)[1:10], each = 2), 
                  "Modelclass" = rep("GLM", times = 20), 
                  "Coefficients" = rep(c("intercept", "dose slope"), times = 10), 
                  as.data.frame(do.call(rbind, reg_glm_list)))

reg_bglm_list <- lapply(lapply(fitlist_bglm, function(x) summary(x)$coefficients[,1:2]), 
                    function(z) do.call(cbind,as.data.frame(z)))
reg_BGLMtable <- cbind("Endpoint" = rep(colnames(tt5)[1:10], each = 2), 
                   "Modelclass" = rep("BGLM", times = 20), 
                   "Coefficients" = rep(c("intercept", "dose slope"), times = 10), 
                   as.data.frame(do.call(rbind, reg_bglm_list)))

reg_coeftable <- rbind(reg_GLMtable,reg_BGLMtable)

write.csv(reg_coeftable, 
          ".\\results\\Example_Analysis_NTP_Regression_Effect_sizes_and_SEs.csv", 
          row.names = FALSE)

# Group comparison models--------------------------------------------------

# change dose varirable to factor
tt5$dosef <- factor(tt5$dose)

## Generalized linear models ----------------------------------------------

## generalized linear models are fitted to each of the tumor class 
t24_glm_gc <- glm(cbind(t24,1-t24)~dosef, data = tt5, family = binomial())
t26_glm_gc <- glm(cbind(t26,1-t26)~dosef, data = tt5, family = binomial())
t27_glm_gc <- glm(cbind(t27,1-t27)~dosef, data = tt5, family = binomial())
t28_glm_gc <- glm(cbind(t28,1-t28)~dosef, data = tt5, family = binomial())
t29_glm_gc <- glm(cbind(t29,1-t29)~dosef, data = tt5, family = binomial())
t36_glm_gc <- glm(cbind(t36,1-t36)~dosef, data = tt5, family = binomial())
t41_glm_gc <- glm(cbind(t41,1-t41)~dosef, data = tt5, family = binomial())
t42_glm_gc <- glm(cbind(t42,1-t42)~dosef, data = tt5, family = binomial())
t59_glm_gc <- glm(cbind(t59,1-t59)~dosef, data = tt5, family = binomial())
t71_glm_gc <- glm(cbind(t71,1-t71)~dosef, data = tt5, family = binomial())

## glm_gcs are combined into list for further analysis
fitlist_glm_gc <- list(t24_glm_gc, t26_glm_gc, t27_glm_gc, t28_glm_gc, t29_glm_gc, 
                    t36_glm_gc, t41_glm_gc, t42_glm_gc, t59_glm_gc, t71_glm_gc)


## Bayesian generalized linear models -------------------------------------

## read "arm" library for bayesian glm_gc model fitting
#install.packages("arm")
library(arm)

## Bayesian generalized linear models are fitted to each of the tumor class 
t24_bglm_gc <- bayesglm(cbind(t24,1-t24)~dosef, data = tt5, family = binomial())
t26_bglm_gc <- bayesglm(cbind(t26,1-t26)~dosef, data = tt5, family = binomial())
t27_bglm_gc <- bayesglm(cbind(t27,1-t27)~dosef, data = tt5, family = binomial())
t28_bglm_gc <- bayesglm(cbind(t28,1-t28)~dosef, data = tt5, family = binomial())
t29_bglm_gc <- bayesglm(cbind(t29,1-t29)~dosef, data = tt5, family = binomial())
t36_bglm_gc <- bayesglm(cbind(t36,1-t36)~dosef, data = tt5, family = binomial())
t41_bglm_gc <- bayesglm(cbind(t41,1-t41)~dosef, data = tt5, family = binomial())
t42_bglm_gc <- bayesglm(cbind(t42,1-t42)~dosef, data = tt5, family = binomial())
t59_bglm_gc <- bayesglm(cbind(t59,1-t59)~dosef, data = tt5, family = binomial())
t71_bglm_gc <- bayesglm(cbind(t71,1-t71)~dosef, data = tt5, family = binomial())

## Bglm_gcs are combined into list for further analysis
fitlist_bglm_gc <- list(t24_bglm_gc, t26_bglm_gc, t27_bglm_gc, t28_bglm_gc, t29_bglm_gc, 
                     t36_bglm_gc, t41_bglm_gc, t42_bglm_gc, t59_bglm_gc, t71_bglm_gc)


## Bonferroni adjustment (GLM) --------------------------------------------

## p-values for each model are extracted from summary 
p_vec_glm_gc <- lapply(fitlist_glm_gc, function(x) summary(x)$coefficients[c(2:4),4])
p_vec_glm_gc <- as.vector(unlist(p_vec_glm_gc))

## adjusted p-values are obtained with p.adjust function 'method = "bonferroni"'
p_adj_glm_gc_bonf <- p.adjust(p_vec_glm_gc, method = "bonferroni")

## The sorting of the p-values corresponds to the order of the models in the list.
p_adj_glm_gc_bonf


## Bonferroni adjustment (Bglm_gc) -------------------------------------------

## p-values for each model are extracted from summary 
p_vec_bglm_gc <- lapply(fitlist_bglm_gc, function(x) summary(x)$coefficients[c(2:4),4])
p_vec_bglm_gc <- as.vector(unlist(p_vec_bglm_gc))

## adjusted p-values are obtained with p.adjust function 'method = "bonferroni"'
p_adj_bglm_gc_bonf <- p.adjust(p_vec_bglm_gc, method = "bonferroni")

## The sorting of the p-values corresponds to the order of the models in the list.
p_adj_bglm_gc_bonf


## MMM adjustment (GLM) ---------------------------------------------------

# Note: The p-value adjustment for the GLMS using the mmm function in the 
# group comparisons does not work because there are groups with zero 
# observations, leading to negative entries in the variance-covariance matrix.

## read "multcomp" library for application of 'Multiple marginal models'
#install.packages("multcomp")
# library(multcomp)
# 
# ## combine respective models into list of class 'mmm'
# mmm_glm_gc <- mmm(t24_glm_gc, t26_glm_gc, t27_glm_gc, t28_glm_gc, t29_glm_gc, 
#                t36_glm_gc, t41_glm_gc, t42_glm_gc, t59_glm_gc, t71_glm_gc)
# 
# ## use function 'glht' with the Dunnett matrix and 'summary' to obtain the adjusted p-values
# mmm_glht_glm_gc <- glht(mmm_glm_gc, 
#                         linfct = mlf(mcp(dosef = "Dunnett")), 
#                         alternative = "two.sided")
# p_adj_glm_gc_mmm <- summary(mmm_glht_glm_gc)
# 
# p_adj_glm_gc_mmm


## MMM adjustment (BGLM) --------------------------------------------------

## read "multcomp" library for application of 'Multiple marginal models'
#install.packages("multcomp")
library(multcomp)

## combine respective models into list of class 'mmm'
mmm_bglm_gc <- mmm(t24_bglm_gc, t26_bglm_gc, t27_bglm_gc, t28_bglm_gc, t29_bglm_gc, 
                t36_bglm_gc, t41_bglm_gc, t42_bglm_gc, t59_bglm_gc, t71_bglm_gc)

## use function 'glht' with the Dunnet matrix and 'summary' to obtain the adjusted p-values
mmm_glht_bglm_gc <- glht(mmm_bglm_gc, 
                         linfct = mlf(mcp(dosef = "Dunnett")), 
                         alternative = "two.sided")
p_adj_bglm_gc_mmm <- summary(mmm_glht_bglm_gc)

p_adj_bglm_gc_mmm


## Bootstrap adjustment (GLM) ---------------------------------------------

## extract the p-values of interest with the summary function and store into vector
l_pvalues_glm_gc <- lapply(fitlist_glm_gc, function(x) summary(x)$coefficients[c(2:4),4])
v_pvalues_glm_gc <- as.vector(unlist(l_pvalues_glm_gc))

## J is the number of endpoints
J <- 10

## specify the matrix where bootstrap results are stored
BPM <- matrix(nrow = 0, ncol = length(v_pvalues_glm_gc))

# set number of bootstrap iterations
nboot <- 10000

for (b in 1:nboot) {
  
  ## resample data for each of the tumor vectors
  RSDAT <- apply(tt5[seq(from = 1, to = J, by = 1)], 2, sample, replace = TRUE)
  RFDAT <- 1-RSDAT
  
  ## create new data.frame with the resampled data and success/occurence and failure/non-occurence
  colnames(RSDAT) <- paste("s",1:J, sep="")
  colnames(RFDAT) <- paste("f",1:J, sep="")
  RSDF <- cbind(as.data.frame(RSDAT), as.data.frame(RFDAT), "dosef"=tt5$dosef)
  
  ## create appropriate lists, so glms can be fitted automatically to all endpoints
  NAMS <- paste("s",1:J, sep="")
  NAMF <- paste("f",1:J, sep="")
  FORMVEC <- paste("cbind(", NAMS, ",",  NAMF, ") ~ dosef", sep="")
  
  FORMLIST <- as.list(FORMVEC)
  names(FORMLIST) <- NAMS
  RSFITLIST <- FORMLIST
  
  ## fit glms for the respective endpoints
  for(k in 1:J){RSFITLIST[[k]] <- glm(FORMLIST[[k]], data=RSDF, family=binomial)}
  
  ## extract p-values of interest and store into vector
  RSSUMLIST <- lapply(RSFITLIST, function(x) summary(x)$coefficients[c(2:4),4])
  RSPVEC <- as.vector(unlist(RSSUMLIST))
  
  ## extract the minimum p-value and check whether it is smaller 
  ## than the respective original unadjusted p-values.  
  ## Note the result for the respective endpoint
  BPM <- rbind(BPM,as.integer(min(RSPVEC) <= v_pvalues_glm))
  
}

## calculate how many times the minimum p value from the respective bootstrap run was smaller 
## than the respective unadjusted p-value

p_adj_glm_gc_bt <- colSums(BPM)/nboot

## The sorting of the p-values corresponds to the order of the models in the list.
p_adj_glm_gc_bt


## Bootstrap adjustment (BGLM) --------------------------------------------

## extract the p-values of interest with the summary function and store into vector
l_pvalues_bglm_gc <- lapply(fitlist_bglm_gc, function(x) summary(x)$coefficients[c(2:4),4])
v_pvalues_bglm_gc <- as.vector(unlist(l_pvalues_bglm_gc))

## J is the number of endpoints
J <- 10

## specify the matrix where bootstrap results are stored
BPM <- matrix(nrow = 0, ncol = length(v_pvalues_bglm))

# set number of bootstrap iterations
nboot <- 10000

for (b in 1:nboot) {
  
  ## resample data for each of the tumor vectors
  RSDAT <- apply(tt5[seq(from = 1, to = J, by = 1)], 2, sample, replace = TRUE)
  RFDAT <- 1-RSDAT
  
  ## create new data.frame with the resampled data and success/occurence and failure/non-occurence
  colnames(RSDAT) <- paste("s",1:J, sep="")
  colnames(RFDAT) <- paste("f",1:J, sep="")
  RSDF <- cbind(as.data.frame(RSDAT), as.data.frame(RFDAT), "dosef"=tt5$dosef)
  
  ## create appropriate lists, so glms can be fitted automatically to all endpoints
  NAMS <- paste("s",1:J, sep="")
  NAMF <- paste("f",1:J, sep="")
  FORMVEC <- paste("cbind(", NAMS, ",",  NAMF, ") ~ dosef", sep="")
  
  FORMLIST <- as.list(FORMVEC)
  names(FORMLIST) <- NAMS
  RSFITLIST <- FORMLIST
  
  ## fit glms for the respective endpoints
  for(k in 1:J){RSFITLIST[[k]] <- bayesglm(FORMLIST[[k]], data=RSDF, family=binomial)}
  
  ## extract p-values of interest and store into vector
  RSSUMLIST <- lapply(RSFITLIST, function(x) summary(x)$coefficients[c(2:4),4])
  RSPVEC <- as.vector(unlist(RSSUMLIST))
  
  ## extract the minimum p-value and check whether it is smaller 
  ## than the respective original unadjusted p-values.  
  ## Note the result for the respective endpoint
  BPM <- rbind(BPM,as.integer(min(RSPVEC) <= v_pvalues_bglm))
  
}

## calculate how many times the minimum p value from the respective bootstrap run was smaller 
## than the respective unadjusted p-value

p_adj_bglm_gc_bt <- colSums(BPM)/nboot

## The sorting of the p-values corresponds to the order of the models in the list.
p_adj_bglm_gc_bt


## Dataframe - Adjusted p-values ------------------------------------------

# Generation of the data frame containing the adjusted p values of all methods 
# (except GLM with MMM) for the group comparison models.

df_gc_pvalues <- 
  data.frame(Endpoint = rep(colnames(tt5)[1:10], each  = 3),
             Comparison = rep(c("37 vs 0", "75 vs 0", "150 vs 0"), times =10),
             Bonferroni = round(p_adj_glm_gc_bonf, 4),
             B_Bonferroni = round(p_adj_bglm_gc_bonf, 4),
             B_MMM = round(p_adj_bglm_gc_mmm$test$pvalues, 4),
             Bootstrap = p_adj_glm_gc_bt,
             B_Bootstrap = p_adj_bglm_gc_bt)

write.csv(df_gc_pvalues, 
          ".\\results\\Example_Analysis_NTP_GroupComparisons_pvalues.csv", 
          row.names = FALSE)

## Dataframe - Effect Sizes and SEs ---------------------------------------

# Generate data frame containing the effect sizes and standard errors 
# for all group comparison models

gc_glm_list <- lapply(lapply(fitlist_glm_gc, function(x) summary(x)$coefficients[2:4,1:2]), 
                       function(z) do.call(cbind,as.data.frame(z)))
gc_GLMtable <- cbind("Endpoint" = rep(colnames(tt5)[1:10], each = 3), 
                      "Modelclass" = rep("GLM", times = 30), 
                      "Coefficients" = rep(c("37 vs 0", "75 vs 0", "150 vs 0"), times =10), 
                      as.data.frame(do.call(rbind, gc_glm_list)))

gc_bglm_list <- lapply(lapply(fitlist_bglm_gc, function(x) summary(x)$coefficients[2:4,1:2]), 
                        function(z) do.call(cbind,as.data.frame(z)))
gc_BGLMtable <- cbind("Endpoint" = rep(colnames(tt5)[1:10], each = 3), 
                       "Modelclass" = rep("BGLM", times = 30), 
                      "Coefficients" = rep(c("37 vs 0", "75 vs 0", "150 vs 0"), times =10), 
                       as.data.frame(do.call(rbind, gc_bglm_list)))

gc_coeftable <- rbind(gc_GLMtable,gc_BGLMtable)

write.csv(gc_coeftable, 
          ".\\results\\Example_Analysis_NTP_GroupComparisons_Effect_sizes_and_SEs.csv", 
          row.names = FALSE)

