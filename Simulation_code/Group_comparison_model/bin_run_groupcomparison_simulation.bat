@echo off
G:
REM Choose path for R
PATH C:\Program Files\R\R-4.1.1\bin\x64;%path%

REM Choose path where Simulation Code is saved
cd ...\Simulation_code\Group_comparison_model

REM Set path where R library is saved
set "R_LIBS=%C:/Program Files/R/R-4.1.1/library"

REM make sure R is also set system environment variables

start bin_gc_ss_2x001_2x0125_2x025_2x0375_2x05.R 1
start bin_gc_ss_4x001_3x025_3x05.R 1
start bin_gc_ss_4x01_3x03_3x05.R 1
start bin_gc_ss_h02_01_ha8_04.R 1
start bin_gc_ss_h02_01_ha8_04_d01.R 1
start bin_gc_ss_h02_01_ha8_04_d02.R 1
start bin_gc_ss_h02_01_ha8_04_s01.R 1
start bin_gc_ss_h02_01_ha8_04_s02.R 1
start bin_gc_ss_h02_002_ha8_04.R 1
start bin_gc_ss_h02_002_ha8_04_d01.R 1
start bin_gc_ss_h02_002_ha8_04_d02.R 1
start bin_gc_ss_h02_002_ha8_04_s01.R 1
start bin_gc_ss_h02_002_ha8_04_s02.R 1
start bin_gc_ss_h02_04_ha8_01_s01.R 1
start bin_gc_ss_h02_04_ha8_01_s02.R 1
start bin_gc_ss_h02_04_ha8_002.R 1
start bin_gc_ss_h02_04_ha8_002_d01.R 1
start bin_gc_ss_h02_04_ha8_002_d02.R 1
start bin_gc_ss_h02_04_ha8_002_s01.R 1
start bin_gc_ss_h02_04_ha8_002_s02.R 1
start bin_gc_ss_h08_01_ha2_04.R 1
start bin_gc_ss_h08_01_ha2_04_d01.R 1
start bin_gc_ss_h08_01_ha2_04_d02.R 1
start bin_gc_ss_h08_01_ha2_04_s01.R 1
start bin_gc_ss_h08_01_ha2_04_s02.R 1
start bin_gc_ss_h08_002_ha2_04.R 1
start bin_gc_ss_h08_002_ha2_04_d01.R 1
start bin_gc_ss_h08_002_ha2_04_d02.R 1
start bin_gc_ss_h08_002_ha2_04_s01.R 1
start bin_gc_ss_h08_002_ha2_04_s02.R 1
start bin_gc_ss_h08_04_ha2_01.R 1
start bin_gc_ss_h08_04_ha2_01_d01.R 1
start bin_gc_ss_h08_04_ha2_01_d02.R 1
start bin_gc_ss_h08_04_ha2_01_s01.R 1
start bin_gc_ss_h08_04_ha2_01_s012.R 1
start bin_gc_ss_h08_04_ha2_002.R 1
start bin_gc_ss_h08_04_ha2_002_d01.R 1
start bin_gc_ss_h08_04_ha2_002_d02.R 1
start bin_gc_ss_h08_04_ha2_002_s01.R 1
start bin_gc_ss_h08_04_ha2_002_s02.R 1