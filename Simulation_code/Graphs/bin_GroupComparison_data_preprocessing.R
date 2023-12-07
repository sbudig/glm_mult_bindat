
# Setup -------------------------------------------------------------------

# set the appropriate wokring directory
setwd(".\\Code_and_Data\\")

# Packages ----------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(patchwork)

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

# 2x001_2x0125_2x025_2x0375_2x05 ------------------------------------------

ss2x001_2x0125_2x025_2x0375_2x05 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_2x001_2x0125_2x025_2x0375_2x05_",
    20
  )

ss2x001_2x0125_2x025_2x0375_2x05 <-
  ss2x001_2x0125_2x025_2x0375_2x05 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h010_5pis_0125")



# 4x001_3x025_3x05 --------------------------------------------------------

ss4x001_3x025_3x05 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_4x001_3x025_3x05_",
    20
  )

ss4x001_3x025_3x05 <- ss4x001_3x025_3x05 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h010_4x001_3x025_3x05")


# 4x01_3x03_3x05 ----------------------------------------------------------

ss4x01_3x03_3x05 <-
  import_mult_csv(".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_4x01_3x03_3x05_",
                  20)

ss4x01_3x03_3x05 <- ss4x01_3x03_3x05 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h010_4x01_3x03_3x05")


# h02_01_ha8_04 -----------------------------------------------------------

ssh02_01_ha8_04 <-
  import_mult_csv(".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_01_ha8_04_",
                  20)

ssh02_01_ha8_04 <- ssh02_01_ha8_04 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_01_ha8_04_s03")


# h02_002_ha8_04 ----------------------------------------------------------

ssh02_002_ha8_04 <-
  import_mult_csv(".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_002_ha8_04_",
                  20)

ssh02_002_ha8_04 <- ssh02_002_ha8_04 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_002_ha8_04_s03")


# h02_04_ha8_01 -----------------------------------------------------------

ssh02_04_ha8_01 <-
  import_mult_csv(".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_04_ha8_01_",
                  20)

ssh02_04_ha8_01 <- ssh02_04_ha8_01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_04_ha8_01_s03")


# h02_04_ha8_002 ----------------------------------------------------------

ssh02_04_ha8_002 <-
  import_mult_csv(".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_04_ha8_002_",
                  20)

ssh02_04_ha8_002 <- ssh02_04_ha8_002 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_04_ha8_002_s03")


# h08_01_ha2_04 -----------------------------------------------------------

ssh08_01_ha2_04 <-
  import_mult_csv(".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_01_ha2_04_",
                  20)

ssh08_01_ha2_04 <- ssh08_01_ha2_04 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_01_ha2_04_s03")


# h08_002_ha2_04 ----------------------------------------------------------

ssh08_002_ha2_04 <-
  import_mult_csv(".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_002_ha2_04_",
                  20)

ssh08_002_ha2_04 <- ssh08_002_ha2_04 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_002_ha2_04_s03")


# h08_04_ha2_01 -----------------------------------------------------------

ssh08_04_ha2_01 <-
  import_mult_csv(".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_04_ha2_01_",
                  20)

ssh08_04_ha2_01 <- ssh08_04_ha2_01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_04_ha2_01_s03")


# h08_04_ha2_002 ----------------------------------------------------------

ssh08_04_ha2_002 <-
  import_mult_csv(".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_04_ha2_002_",
                  20)

ssh08_04_ha2_002 <- ssh08_04_ha2_002 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_04_ha2_002_s03")



# h02_01_ha8_04_d01 -------------------------------------------------------

ssh02_01_ha8_04_d01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_01_ha8_04_d01_",
    20
  )

ssh02_01_ha8_04_d01 <- ssh02_01_ha8_04_d01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_01_ha8_04_d01")


# bin_gc_ss_h02_01_ha8_04_d02 ---------------------------------------------

ss_h02_01_ha8_04_d02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_01_ha8_04_d02_",
    20
  )

ss_h02_01_ha8_04_d02 <- ss_h02_01_ha8_04_d02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_01_ha8_04_d02")


# bin_gc_ss_h02_01_ha8_04_s01 ---------------------------------------------

ss_h02_01_ha8_04_s01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_01_ha8_04_s01_",
    20
  )

ss_h02_01_ha8_04_s01 <- ss_h02_01_ha8_04_s01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_01_ha8_04_s01")


# bin_gc_ss_h02_01_ha8_04_s02 ---------------------------------------------

ss_h02_01_ha8_04_s02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_01_ha8_04_s02_",
    20
  )

ss_h02_01_ha8_04_s02 <- ss_h02_01_ha8_04_s02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_01_ha8_04_s02")


# bin_gc_ss_h02_002_ha8_04_d01 --------------------------------------------

ss_h02_002_ha8_04_d01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_002_ha8_04_d01_",
    20
  )

ss_h02_002_ha8_04_d01 <- ss_h02_002_ha8_04_d01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_002_ha8_04_d01")


# bin_gc_ss_h02_002_ha8_04_d02 --------------------------------------------

ss_h02_002_ha8_04_d02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_002_ha8_04_d02_",
    20
  )

ss_h02_002_ha8_04_d02 <- ss_h02_002_ha8_04_d02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_002_ha8_04_d02")


# bin_gc_ss_h02_002_ha8_04_s01 --------------------------------------------

ss_h02_002_ha8_04_s01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_002_ha8_04_s01_",
    20
  )

ss_h02_002_ha8_04_s01 <- ss_h02_002_ha8_04_s01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_002_ha8_04_s01")


# bin_gc_ss_h02_002_ha8_04_s02 --------------------------------------------

ss_h02_002_ha8_04_s02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_002_ha8_04_s02_",
    20
  )

ss_h02_002_ha8_04_s02 <- ss_h02_002_ha8_04_s02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_002_ha8_04_s02")


# bin_gc_ss_h02_04_ha8_01_d01 ---------------------------------------------

ss_h02_04_ha8_01_d01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_04_ha8_01_d01_",
    20
  )

ss_h02_04_ha8_01_d01 <- ss_h02_04_ha8_01_d01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_04_ha8_01_d01")


# bin_gc_ss_h02_04_ha8_01_d02 ---------------------------------------------

ss_h02_04_ha8_01_d02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_04_ha8_01_d02_",
    20
  )

ss_h02_04_ha8_01_d02 <- ss_h02_04_ha8_01_d02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_04_ha8_01_d02")


# bin_gc_ss_h02_04_ha8_01_s01 ---------------------------------------------

ss_h02_04_ha8_01_s01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_04_ha8_01_s01_",
    20
  )

ss_h02_04_ha8_01_s01 <- ss_h02_04_ha8_01_s01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_04_ha8_01_s01")


# bin_gc_ss_h02_04_ha8_01_s02 ---------------------------------------------

ss_h02_04_ha8_01_s02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_04_ha8_01_s02_",
    20
  )

ss_h02_04_ha8_01_s02 <- ss_h02_04_ha8_01_s02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_04_ha8_01_s02")


# bin_gc_ss_h02_04_ha8_002_d01 --------------------------------------------

ss_h02_04_ha8_002_d01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_04_ha8_002_d01_",
    20
  )

ss_h02_04_ha8_002_d01 <- ss_h02_04_ha8_002_d01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_04_ha8_002_d01")


# bin_gc_ss_h02_04_ha8_002_d02 --------------------------------------------

ss_h02_04_ha8_002_d02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_04_ha8_002_d02_",
    20
  )

ss_h02_04_ha8_002_d02 <- ss_h02_04_ha8_002_d02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_04_ha8_002_d02")


# bin_gc_ss_h02_04_ha8_002_s01 --------------------------------------------

ss_h02_04_ha8_002_s01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_04_ha8_002_s01_",
    20
  )

ss_h02_04_ha8_002_s01 <- ss_h02_04_ha8_002_s01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_04_ha8_002_s01")


# bin_gc_ss_h02_04_ha8_002_s02 --------------------------------------------

ss_h02_04_ha8_002_s02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h02_04_ha8_002_s02_",
    20
  )

ss_h02_04_ha8_002_s02 <- ss_h02_04_ha8_002_s02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h02_04_ha8_002_s02")


# bin_gc_ss_h08_01_ha2_04_d01 ---------------------------------------------

ss_h08_01_ha2_04_d01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_01_ha2_04_d01_",
    20
  )

ss_h08_01_ha2_04_d01 <- ss_h08_01_ha2_04_d01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_01_ha2_04_d01")


# bin_gc_ss_h08_01_ha2_04_d02 ---------------------------------------------

ss_h08_01_ha2_04_d02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_01_ha2_04_d02_",
    20
  )

ss_h08_01_ha2_04_d02 <- ss_h08_01_ha2_04_d02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_01_ha2_04_d02")


# bin_gc_ss_h08_01_ha2_04_s01 ---------------------------------------------

ss_h08_01_ha2_04_s01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_01_ha2_04_s01_",
    20
  )

ss_h08_01_ha2_04_s01 <- ss_h08_01_ha2_04_s01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_01_ha2_04_s01")


# bin_gc_ss_h08_01_ha2_04_s02 ---------------------------------------------

ss_h08_01_ha2_04_s02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_01_ha2_04_s02_",
    20
  )

ss_h08_01_ha2_04_s02 <- ss_h08_01_ha2_04_s02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_01_ha2_04_s02")


# bin_gc_ss_h08_002_ha2_04_d01 --------------------------------------------

ss_h08_002_ha2_04_d01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_002_ha2_04_d01_",
    20
  )

ss_h08_002_ha2_04_d01 <- ss_h08_002_ha2_04_d01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_002_ha2_04_d01")


# bin_gc_ss_h08_002_ha2_04_d02 --------------------------------------------

ss_h08_002_ha2_04_d02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_002_ha2_04_d02_",
    20
  )

ss_h08_002_ha2_04_d02 <- ss_h08_002_ha2_04_d02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_002_ha2_04_d02")


# bin_gc_ss_h08_002_ha2_04_s01 --------------------------------------------

ss_h08_002_ha2_04_s01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_002_ha2_04_s01_",
    20
  )

ss_h08_002_ha2_04_s01 <- ss_h08_002_ha2_04_s01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_002_ha2_04_s01")


# bin_gc_ss_h08_002_ha2_04_s02 --------------------------------------------

ss_h08_002_ha2_04_s02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_002_ha2_04_s02_",
    20
  )

ss_h08_002_ha2_04_s02 <- ss_h08_002_ha2_04_s02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_002_ha2_04_s02")


# bin_gc_ss_h08_04_ha2_01_d01 ---------------------------------------------

ss_h08_04_ha2_01_d01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_04_ha2_01_d01_",
    20
  )

ss_h08_04_ha2_01_d01 <- ss_h08_04_ha2_01_d01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_04_ha2_01_d01")


# bin_gc_ss_h08_04_ha2_01_d02 ---------------------------------------------

ss_h08_04_ha2_01_d02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_04_ha2_01_d02_",
    20
  )

ss_h08_04_ha2_01_d02 <- ss_h08_04_ha2_01_d02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_04_ha2_01_d02")


# bin_gc_ss_h08_04_ha2_01_s01 ---------------------------------------------

ss_h08_04_ha2_01_s01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_04_ha2_01_s01_",
    20
  )

ss_h08_04_ha2_01_s01 <- ss_h08_04_ha2_01_s01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_04_ha2_01_s01")


# bin_gc_ss_h08_04_ha2_01_s02 ---------------------------------------------

ss_h08_04_ha2_01_s02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_04_ha2_01_s02_",
    20
  )

ss_h08_04_ha2_01_s02 <- ss_h08_04_ha2_01_s02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_04_ha2_01_s02")


# bin_gc_ss_h08_04_ha2_002_d01 --------------------------------------------

ss_h08_04_ha2_002_d01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_04_ha2_002_d01_",
    20
  )

ss_h08_04_ha2_002_d01 <- ss_h08_04_ha2_002_d01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_04_ha2_002_d01")


# bin_gc_ss_h08_04_ha2_002_d02 --------------------------------------------

ss_h08_04_ha2_002_d02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_04_ha2_002_d02_",
    20
  )

ss_h08_04_ha2_002_d02 <- ss_h08_04_ha2_002_d02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_04_ha2_002_d02")


# bin_gc_ss_h08_04_ha2_002_s01 --------------------------------------------

ss_h08_04_ha2_002_s01 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_04_ha2_002_s01_",
    20
  )

ss_h08_04_ha2_002_s01 <- ss_h08_04_ha2_002_s01 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_04_ha2_002_s01")


# bin_gc_ss_h08_04_ha2_002_s02 --------------------------------------------

ss_h08_04_ha2_002_s02 <-
  import_mult_csv(
    ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_h08_04_ha2_002_s02_",
    20
  )

ss_h08_04_ha2_002_s02 <- ss_h08_04_ha2_002_s02 %>%
  group_by(method, model, ntrt, corr) %>%
  summarise(
    nsim = sum(nsim),
    nrej = sum(nrej),
    nreji = sum(nreji),
    nrejcor = sum(nrejcor),
    npow = sum(npow),
    npowi = sum(npowi),
    npowcor = sum(npowcor),
    nerror = sum(nerror),
    nfail = sum(nfail),
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
    error_p = nerror / nsim,
    fail_p = nfail / nsim
  ) %>%
  mutate(setting = "h08_04_ha2_002_s02")


# All Combined ------------------------------------------------------------

all_gc_ss <- rbind(
  ss2x001_2x0125_2x025_2x0375_2x05,
  ss4x001_3x025_3x05,
  ss4x01_3x03_3x05,
  ssh02_01_ha8_04,
  ssh02_002_ha8_04,
  ssh02_04_ha8_01,
  ssh02_04_ha8_002,
  ssh08_01_ha2_04,
  ssh08_002_ha2_04,
  ssh08_04_ha2_01,
  ssh08_04_ha2_002,
  ssh02_01_ha8_04_d01,
  ss_h02_01_ha8_04_d02,
  ss_h02_01_ha8_04_s01,
  ss_h02_01_ha8_04_s02,
  ss_h02_002_ha8_04_d01,
  ss_h02_002_ha8_04_d02,
  ss_h02_002_ha8_04_s01,
  ss_h02_002_ha8_04_s02,
  ss_h02_04_ha8_01_d01,
  ss_h02_04_ha8_01_d02,
  ss_h02_04_ha8_01_s01,
  ss_h02_04_ha8_01_s02,
  ss_h02_04_ha8_002_d01,
  ss_h02_04_ha8_002_d02,
  ss_h02_04_ha8_002_s01,
  ss_h02_04_ha8_002_s02,
  ss_h08_01_ha2_04_d01,
  ss_h08_01_ha2_04_d02,
  ss_h08_01_ha2_04_s01,
  ss_h08_01_ha2_04_s02,
  ss_h08_002_ha2_04_d01,
  ss_h08_002_ha2_04_d02,
  ss_h08_002_ha2_04_s01,
  ss_h08_002_ha2_04_s02,
  ss_h08_04_ha2_01_d01,
  ss_h08_04_ha2_01_d02,
  ss_h08_04_ha2_01_s01,
  ss_h08_04_ha2_01_s02,
  ss_h08_04_ha2_002_d01,
  ss_h08_04_ha2_002_d02,
  ss_h08_04_ha2_002_s01,
  ss_h08_04_ha2_002_s02
)

all_gc_ss <- all_gc_ss %>% mutate(N = ntrt * 4, cor = corr)

all_gc_ss$model <- factor(all_gc_ss$model)
all_gc_ss$model <- relevel(all_gc_ss$model, ref = "glm")

write.csv(all_gc_ss, ".\\Simulation_results\\Group_comparison_model\\bin_gc_ss_all.csv", row.names = FALSE)
