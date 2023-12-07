@echo off
G:
REM Choose path for R
PATH C:\Program Files\R\R-4.1.1\bin\x64;%path%

REM Choose path where Simulation Code is saved
cd ...\Simulation_code\Regression_model

REM Set path where R library is saved
set "R_LIBS=%C:/Program Files/R/R-4.1.1/library"

REM make sure R is also set system environment variables

start Rscript bin_reg_FWER_bonf_bglm.R
start Rscript bin_reg_FWER_bonf_glm.R
start Rscript bin_reg_FWER_bt_bglm.R
start Rscript bin_reg_FWER_bt_glm.R
start Rscript bin_reg_FWER_mmm_bglm.R
start Rscript bin_reg_FWER_mmm_glm.R
start Rscript bin_reg_POWER_bonf_bglm.R
start Rscript bin_reg_POWER_bonf_glm.R
start Rscript bin_reg_POWER_bt_bglm.R
start Rscript bin_reg_POWER_bt_glm.R
start Rscript bin_reg_POWER_mmm_bglm.R
start Rscript bin_reg_POWER_mmm_glm.R