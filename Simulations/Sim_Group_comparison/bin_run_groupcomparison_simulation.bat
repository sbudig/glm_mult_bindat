@echo off

REM If a specific version of R is required, the corresponding path can be selected
REM must be removed so that the line is taken into account
REM PATH C:\Program Files\R\R-4.1.1\bin\x64;%path%

REM Path to a specific R library can be set
REM set "R_LIBS=%C:/Program Files/R/R-4.1.1/library"

REM Select the path where the simulation code is stored
REM As specified now, the path in which this batch file is located is automatically taken
cd /d "%~dp0"

REM make sure R is also set as system environment variable
REM The number behind the scripts indiciate with which numbers extensions the results are saved
REM As the simulation takes quite a long time, it is advisable to split the simulation runs and 
REM use several numbers under which the partial results are saved.

start Rscript bin_gc_ss_2x001_2x0125_2x025_2x0375_2x05.R 1
start Rscript bin_gc_ss_4x001_3x025_3x05.R 1
start Rscript bin_gc_ss_4x01_3x03_3x05.R 1
start Rscript bin_gc_ss_h02_01_ha8_04.R 1
start Rscript bin_gc_ss_h02_01_ha8_04_d01.R 1
start Rscript bin_gc_ss_h02_01_ha8_04_d02.R 1
start Rscript bin_gc_ss_h02_01_ha8_04_s01.R 1
start Rscript bin_gc_ss_h02_01_ha8_04_s02.R 1
start Rscript bin_gc_ss_h02_002_ha8_04.R 1
start Rscript bin_gc_ss_h02_002_ha8_04_d01.R 1
start Rscript bin_gc_ss_h02_002_ha8_04_d02.R 1
start Rscript bin_gc_ss_h02_002_ha8_04_s01.R 1
start Rscript bin_gc_ss_h02_002_ha8_04_s02.R 1
start Rscript bin_gc_ss_h02_04_ha8_01_s01.R 1
start Rscript bin_gc_ss_h02_04_ha8_01_s02.R 1
start Rscript bin_gc_ss_h02_04_ha8_002.R 1
start Rscript bin_gc_ss_h02_04_ha8_002_d01.R 1
start Rscript bin_gc_ss_h02_04_ha8_002_d02.R 1
start Rscript bin_gc_ss_h02_04_ha8_002_s01.R 1
start Rscript bin_gc_ss_h02_04_ha8_002_s02.R 1
start Rscript bin_gc_ss_h08_01_ha2_04.R 1
start Rscript bin_gc_ss_h08_01_ha2_04_d01.R 1
start Rscript bin_gc_ss_h08_01_ha2_04_d02.R 1
start Rscript bin_gc_ss_h08_01_ha2_04_s01.R 1
start Rscript bin_gc_ss_h08_01_ha2_04_s02.R 1
start Rscript bin_gc_ss_h08_002_ha2_04.R 1
start Rscript bin_gc_ss_h08_002_ha2_04_d01.R 1
start Rscript bin_gc_ss_h08_002_ha2_04_d02.R 1
start Rscript bin_gc_ss_h08_002_ha2_04_s01.R 1
start Rscript bin_gc_ss_h08_002_ha2_04_s02.R 1
start Rscript bin_gc_ss_h08_04_ha2_01.R 1
start Rscript bin_gc_ss_h08_04_ha2_01_d01.R 1
start Rscript bin_gc_ss_h08_04_ha2_01_d02.R 1
start Rscript bin_gc_ss_h08_04_ha2_01_s01.R 1
start Rscript bin_gc_ss_h08_04_ha2_01_s012.R 1
start Rscript bin_gc_ss_h08_04_ha2_002.R 1
start Rscript bin_gc_ss_h08_04_ha2_002_d01.R 1
start Rscript bin_gc_ss_h08_04_ha2_002_d02.R 1
start Rscript bin_gc_ss_h08_04_ha2_002_s01.R 1
start Rscript bin_gc_ss_h08_04_ha2_002_s02.R 1