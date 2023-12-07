
# Setup -------------------------------------------------------------------

# set the appropriate wokring directory
setwd(".\\Code_and_Data\\")

# Packages ----------------------------------------------------------------

library(ggplot2)
library(tidyverse)

# Helper Functions --------------------------------------------------------

import_mult_csv <- function(name, to) {
  dat <- NULL
  
  for (i in 1:to) {
    tryCatch({
      dat <- rbind(dat, read.csv(paste0(name, i, ".csv")))
    },
    error = function(e) {
      
    })
  }
  return(dat)
}

# Data preprocessing for FWER results -------------------------------------

## Bonferroni Glm ---------------------------------------------------------

fwer_bonf_glm <-
  read.csv(".\\Simulation_results\\Regression_model\\bin_FWER_prv_bonf_glm_1.csv")

fwer_bonf_glm <-
  fwer_bonf_glm %>% group_by(method, model, pis, nresp, ntrt, corr, nsim) %>%
  summarise(
    percent = nrej / nsim,
    percenti = nreji / nsim,
    pererror = nerror / nsim,
    perfail = nfail / nsim
  )

## Bonferroni BGlm ---------------------------------------------------------

fwer_bonf_bglm <-
  read.csv(".\\Simulation_results\\Regression_model\\bin_FWER_prv_bonf_bglm_1.csv")

fwer_bonf_bglm <-
  fwer_bonf_bglm %>% group_by(method, model, pis, nresp, ntrt, corr, nsim) %>%
  summarise(
    percent = nrej / nsim,
    percenti = nreji / nsim,
    pererror = nerror / nsim,
    perfail = nfail / nsim
  )
fwer_bonf_bglm <-
  fwer_bonf_bglm %>% mutate(model = str_replace(model, "glm", "bglm"))

## MMM Glm ----------------------------------------------------------------

fwer_mmm_glm <-
  read.csv(".\\Simulation_results\\Regression_model\\bin_FWER_prv_mmm_glm_1.csv")

fwer_mmm_glm <-
  fwer_mmm_glm %>% group_by(method, model, pis, nresp, ntrt, corr, nsim) %>%
  summarise(
    percent = nrej / nsim,
    percenti = nreji / nsim,
    pererror = nerror / nsim,
    perfail = nfail / nsim
  )

## MMM BGlm ----------------------------------------------------------------

fwer_mmm_bglm <-
  read.csv(".\\Simulation_results\\Regression_model\\bin_FWER_prv_mmm_bglm_1.csv")

fwer_mmm_bglm <-
  fwer_mmm_bglm %>% group_by(method, model, pis, nresp, ntrt, corr, nsim) %>%
  summarise(
    percent = nrej / nsim,
    percenti = nreji / nsim,
    pererror = nerror / nsim,
    perfail = nfail / nsim
  )
fwer_mmm_bglm <-
  fwer_mmm_bglm %>% mutate(model = str_replace(model, "glm", "bglm"))

## Bootstrap Glm ----------------------------------------------------------

# The simulation runs were split for the bootstrap,
# as everything would otherwise have taken several weeks in one simulation.
# The import_mult_csv() function is therefore used to load the split csv files.
# For the sake of clarity, we have already packed all the individual csv files
# of the simulation into one file.
# A warning occurs (and can be ignored), when there is no file with the respective number.

fwer_bt_glm <-
  import_mult_csv(".\\Simulation_results\\Regression_model\\bin_FWER_prv_bt_glm_",
                  80)

fwer_bt_glm <-
  fwer_bt_glm %>% group_by(method, model, pis, nresp, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nerror = sum(nerror),
    nfail = sum(nfail),
    .groups = "keep"
  ) %>%
  group_by(method, model, pis, nresp, ntrt, corr, nsim) %>%
  summarise(
    percent = nrej / nsim,
    percenti = nreji / nsim,
    pererror = nerror / nsim,
    perfail = nfail / nsim
  )

## Bootstrap Bglm ----------------------------------------------------------

fwer_bt_bglm <-
  import_mult_csv(".\\Simulation_results\\Regression_model\\bin_FWER_prv_bt_bglm_",
                  80)

fwer_bt_bglm <-
  fwer_bt_bglm %>% group_by(method, model, pis, nresp, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nerror = sum(nerror),
    nfail = sum(nfail),
    .groups = "keep"
  ) %>%
  group_by(method, model, pis, nresp, ntrt, corr, nsim) %>%
  summarise(
    percent = nrej / nsim,
    percenti = nreji / nsim,
    pererror = nerror / nsim,
    perfail = nfail / nsim
  )

## All Combined FWER -------------------------------------------------------

bin_reg_FWER <-
  rbind(
    fwer_bonf_glm,
    fwer_mmm_glm,
    fwer_bt_glm,
    fwer_bt_bglm,
    fwer_bonf_bglm,
    fwer_mmm_bglm
  )

bin_reg_FWER <- bin_reg_FWER %>% mutate(N = ntrt * 4, J = nresp, cor = corr)

write.csv(
  bin_reg_FWER,
  ".\\Simulation_results\\Regression_model\\bin_reg_FWER_all.csv",
  row.names = FALSE
)


# Data preprocessing for Power results ------------------------------------

## Bonferroni Glm ---------------------------------------------------------

power_bonf_glm <-
  read.csv(".\\Simulation_results\\Regression_model\\bin_POWER_prv_bonf_glm_1.csv")

power_bonf_glm$nep_h0 <- as.factor(power_bonf_glm$nep_h0)
power_bonf_glm$nep_ha <- as.factor(power_bonf_glm$nep_ha)
power_bonf_glm$n_h0_ha <-
  power_bonf_glm$nep_h0:power_bonf_glm$nep_ha

power_bonf_glm <-
  power_bonf_glm %>%  group_by(method,
                               model,
                               nsim,
                               ntrt,
                               pi0,
                               nep_h0,
                               nep_ha,
                               slope,
                               corr,
                               n_h0_ha) %>%
  summarise(
    rej_p = nrej / nsim,
    reji_p = nreji / nsim,
    rejcor_p = nrejcor / nsim,
    pow_p = npow / nsim,
    powi_p = npowi / nsim,
    powcor_p = npowcor / nsim,
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  )

## Bonferroni BGlm ---------------------------------------------------------

power_bonf_bglm <-
  read.csv(".\\Simulation_results\\Regression_model\\bin_POWER_prv_bonf_bglm_1.csv")

power_bonf_bglm <-
  power_bonf_bglm %>% mutate(model = str_replace(model, "glm", "bglm"))

power_bonf_bglm$nep_h0 <- as.factor(power_bonf_bglm$nep_h0)
power_bonf_bglm$nep_ha <- as.factor(power_bonf_bglm$nep_ha)
power_bonf_bglm$n_h0_ha <-
  power_bonf_bglm$nep_h0:power_bonf_bglm$nep_ha

power_bonf_bglm <-
  power_bonf_bglm %>%  group_by(method,
                                model,
                                nsim,
                                ntrt,
                                pi0,
                                nep_h0,
                                nep_ha,
                                slope,
                                corr,
                                n_h0_ha) %>%
  summarise(
    rej_p = nrej / nsim,
    reji_p = nreji / nsim,
    rejcor_p = nrejcor / nsim,
    pow_p = npow / nsim,
    powi_p = npowi / nsim,
    powcor_p = npowcor / nsim,
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  )

## MMM Glm ----------------------------------------------------------------

power_mmm_glm <-
  read.csv(".\\Simulation_results\\Regression_model\\bin_POWER_prv_mmm_glm_1.csv")

power_mmm_glm$nep_h0 <- as.factor(power_mmm_glm$nep_h0)
power_mmm_glm$nep_ha <- as.factor(power_mmm_glm$nep_ha)
power_mmm_glm$n_h0_ha <- power_mmm_glm$nep_h0:power_mmm_glm$nep_ha

power_mmm_glm <-
  power_mmm_glm %>%  group_by(method,
                              model,
                              nsim,
                              ntrt,
                              pi0,
                              nep_h0,
                              nep_ha,
                              slope,
                              corr,
                              n_h0_ha) %>%
  summarise(
    rej_p = nrej / nsim,
    reji_p = nreji / nsim,
    rejcor_p = nrejcor / nsim,
    pow_p = npow / nsim,
    powi_p = npowi / nsim,
    powcor_p = npowcor / nsim,
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  )

## MMM BGlm ----------------------------------------------------------------

power_mmm_bglm <-
  read.csv(".\\Simulation_results\\Regression_model\\bin_POWER_prv_mmm_bglm_1.csv")

power_mmm_bglm <-
  power_mmm_bglm %>% group_by(method, model, ntrt, pi0, nep_h0, nep_ha, slope, corr) %>%
  summarise(
    nsim = sum(nsim),
    nfail = sum(nfail),
    nerror = sum(nerror),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    rej_p = nrej / nsim,
    reji_p = nreji / nsim,
    rejcor_p = nrejcor / nsim,
    pow_p = npow / nsim,
    powi_p = npowi / nsim,
    powcor_p = npowcor / nsim
  )

power_mmm_bglm$nep_h0 <- as.factor(power_mmm_bglm$nep_h0)
power_mmm_bglm$nep_ha <- as.factor(power_mmm_bglm$nep_ha)
power_mmm_bglm$n_h0_ha <-
  power_mmm_bglm$nep_h0:power_mmm_bglm$nep_ha
power_mmm_bglm <-
  power_mmm_bglm %>% mutate(model = str_replace(model, "glm", "bglm"))

## Bootstrap Glm ----------------------------------------------------------

# The simulation runs were split for the bootstrap,
# as everything would otherwise have taken several weeks in one simulation.
# The import_mult_csv() function is therefore used to load the split csv files.
# For the sake of clarity, we have already packed all the individual csv files
# of the simulation into one file.
# A warning occurs (and can be ignored), when there is no file with the respective number.

power_bt_glm <-
  import_mult_csv(".\\Simulation_results\\Regression_model\\bin_POWER_prv_bt_glm_",
                  60)

power_bt_glm <-
  power_bt_glm %>% group_by(method, model, ntrt, pi0, nep_h0, nep_ha, slope, corr) %>%
  summarise(
    nsim = sum(nsim),
    nfail = sum(nfail),
    nerror = sum(nerror),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    .groups = "keep"
  ) %>%
  summarise(
    nsim = nsim,
    rej_p = nrej / nsim,
    reji_p = nreji / nsim,
    rejcor_p = nrejcor / nsim,
    pow_p = npow / nsim,
    powi_p = npowi / nsim,
    powcor_p = npowcor / nsim,
    .groups = "keep"
  )

power_bt_glm$nep_h0 <- as.factor(power_bt_glm$nep_h0)
power_bt_glm$nep_ha <- as.factor(power_bt_glm$nep_ha)
power_bt_glm$n_h0_ha <- power_bt_glm$nep_h0:power_bt_glm$nep_ha

## Bootstrap BGlm ---------------------------------------------------------

power_bt_bglm <-
  import_mult_csv(".\\Simulation_results\\Regression_model\\bin_POWER_prv_bt_bglm_",
                  60)

power_bt_bglm <-
  power_bt_bglm %>% group_by(method, model, ntrt, pi0, nep_h0, nep_ha, slope, corr) %>%
  summarise(
    nsim = sum(nsim),
    nfail = sum(nfail),
    nerror = sum(nerror),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    .groups = "keep"
  ) %>%
  summarise(
    nsim = nsim,
    rej_p = nrej / nsim,
    reji_p = nreji / nsim,
    rejcor_p = nrejcor / nsim,
    pow_p = npow / nsim,
    powi_p = npowi / nsim,
    powcor_p = npowcor / nsim,
    .groups = "keep"
  )

power_bt_bglm$nep_h0 <- as.factor(power_bt_bglm$nep_h0)
power_bt_bglm$nep_ha <- as.factor(power_bt_bglm$nep_ha)
power_bt_bglm$n_h0_ha <- power_bt_bglm$nep_h0:power_bt_bglm$nep_ha

## All Combined Power -----------------------------------------------------

bin_reg_POWER <-
  rbind(
    power_bonf_glm,
    power_mmm_glm,
    power_bt_glm,
    power_bonf_bglm,
    power_mmm_bglm,
    power_bt_bglm
  )
bin_reg_POWER <-
  bin_reg_POWER %>% mutate(N = ntrt * 4, J0_JA = n_h0_ha, cor = corr)

write.csv(
  bin_reg_POWER,
  ".\\Simulation_results\\Regression_model\\bin_reg_POWER_all.csv",
  row.names = FALSE
)

