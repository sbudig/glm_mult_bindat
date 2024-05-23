
# Plot helper functions ---------------------------------------------------

# Removes decimals from zero
prettyZero <- function(l) {
  max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T) - 1
  lnew = formatC(
    l,
    replace.zero = T,
    zero.print = "0",
    digits = max.decimals,
    format = "f",
    preserve.width = T
  )
  return(lnew)
}

# FWER Regression model functions -----------------------------------------

# read and append all intermediate results from FWER simulations of regression models
r_a_a_bin_reg_FWER_csvs <- function(directory) {
  # List all files in the directory
  files <- list.files(directory, pattern = "^bin_reg_FWER.*\\.csv$", full.names = TRUE)
  
  # Read each file and store in a list
  data_list <- lapply(files, read.csv)
  
  # Combine all data frames in the list into one data frame
  combined_data <- do.call(rbind, data_list)
  
  return(combined_data)
}

# summarise and calculate percentages for FWER intermediate results
s_a_cp_bin_reg_FWER <- function(dataframe) {
  
  df_sum_cp <-
    dataframe %>% group_by(method, model, pis, nresp, ntrt, corr) %>%
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
    )  %>% 
    mutate(N = ntrt * 4, J = nresp, cor = corr)
  
  return(df_sum_cp)
}

# Power regression model functions ----------------------------------------

# read and append all intermediate results from Power simulations of regression models
r_a_a_bin_reg_Power_csvs <- function(directory) {
  # List all files in the directory
  files <- list.files(directory, pattern = "^bin_reg_POWER.*\\.csv$", full.names = TRUE)
  
  # Read each file and store in a list
  data_list <- lapply(files, read.csv)
  
  # Combine all data frames in the list into one data frame
  combined_data <- do.call(rbind, data_list)
  
  return(combined_data)
}

# summarise and calculate percentages for Power intermediate results
s_a_cp_bin_reg_Power <- function(dataframe) {
  
  df_sum_cp <-
    dataframe %>% group_by(method, model, ntrt, pi0, nep_h0, nep_ha, slope, corr) %>%
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
    ) %>% 
    mutate(N = ntrt * 4, 
           cor = corr,
           nep_h0 = factor(nep_h0),
           nep_ha = factor(nep_ha),
           n_h0_ha = nep_h0:nep_ha,
           J0_JA = nep_h0:nep_ha)

  return(df_sum_cp)
}

# Group Comparisons model functions ----------------------------------------

# read and append all intermediate results from GC models
r_a_a_bin_gc_ss_csvs <- function(directory) {
  # Get a list of all CSV files in the specified folder
  files <- list.files(directory, pattern = "bin_gc_ss.*\\.csv$", full.names = TRUE)
  
  # Initialize an empty list to store data frames
  df_list <- list()
  
  # Read each CSV file and add the extracted name column
  for (file in files) {
    df <- read.csv(file)
    df$setting <- gsub("bin_gc_ss_(.*?)_(\\d+)\\.csv", "\\1", basename(file))
    df_list[[length(df_list) + 1]] <- df
  }
  
  # Combine all data frames into a single data frame
  combined_data <- do.call(rbind, df_list)
  
  return(combined_data)
}

# summarise and calculate percentages for Power intermediate results
s_a_cp_bin_gc_ss <- function(dataframe) {
  
  df_sum_cp <- dataframe %>%
    group_by(method, model, ntrt, corr, setting) %>%
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
    mutate(setting = case_when(setting == '2x001_2x0125_2x025_2x0375_2x05' ~ 'h010_5pis_0125',
                               setting == '4x001_3x025_3x05' ~ 'h010_4x001_3x025_3x05',
                               setting == '4x01_3x03_3x05' ~ 'h010_4x01_3x03_3x05',
                               setting == 'h02_01_ha8_04' ~ 'h02_01_ha8_04_s03',
                               setting == 'h02_002_ha8_04' ~ 'h02_002_ha8_04_s03',
                               setting == 'h02_04_ha8_01' ~ 'h02_04_ha8_01_s03',
                               setting == 'h02_04_ha8_002' ~ 'h02_04_ha8_002_s03',
                               setting == 'h08_01_ha2_04' ~ 'h08_01_ha2_04_s03',
                               setting == 'h08_002_ha2_04' ~ 'h08_002_ha2_04_s03',
                               setting == 'h08_04_ha2_01' ~ 'h08_04_ha2_01_s03',
                               setting == 'h08_04_ha2_002' ~ 'h08_04_ha2_002_s03',
                               setting == 'h02_01_ha8_04_d01' ~ 'h02_01_ha8_04_d01',
                               setting == 'h02_01_ha8_04_d02' ~ 'h02_01_ha8_04_d02',
                               setting == 'h02_01_ha8_04_s01' ~ 'h02_01_ha8_04_s01',
                               setting == 'h02_01_ha8_04_s02' ~ 'h02_01_ha8_04_s02',
                               setting == 'h02_002_ha8_04_d01' ~ 'h02_002_ha8_04_d01',
                               setting == 'h02_002_ha8_04_d02' ~ 'h02_002_ha8_04_d02',
                               setting == 'h02_002_ha8_04_s01' ~ 'h02_002_ha8_04_s01',
                               setting == 'h02_002_ha8_04_s02' ~ 'h02_002_ha8_04_s02',
                               setting == 'h02_04_ha8_01_d01' ~ 'h02_04_ha8_01_d01',
                               setting == 'h02_04_ha8_01_d02' ~ 'h02_04_ha8_01_d02',
                               .default = as.character(setting))) %>%
    mutate(N = ntrt * 4, cor = corr) 
  
  df_sum_cp$model <- factor(df_sum_cp$model)
  df_sum_cp$model <- relevel(df_sum_cp$model, ref = "glm")
  
  return(df_sum_cp)
}

# long to wide format

l_t_w_bin_gc_ss <- function(longdataframe){
  
  # if not all settings have the same number of simulations,
  # these must be set equal or the column nsim must be removed
  # so that the change to the wide format works.
  # This doesnt influence any of the results.
  
  longdataframe <- longdataframe %>% mutate(nsim = longdataframe$nsim[1])
  
  widedataframe <- longdataframe %>%
    group_by(method) %>%
    pivot_wider(
      names_from = method,
      values_from = c(rej_p, reji_p, rejcor_p,
                      pow_p, powi_p, powcor_p,
                      error_p, fail_p)
    ) %>%
    type.convert(as.is = TRUE)
  
  return(widedataframe)
}
