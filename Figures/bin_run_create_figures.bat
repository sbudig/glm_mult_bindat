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

start Rscript bin_figures.R