Supplementary information / reproducible research files for the manuscript 
Title: "Simultaneous inference of multiple binary endpoints in biomedical
research: small sample properties of multiple marginal models
and a resampling approach"

Authors: Budig, S., Jung, K., Hasler, M., Schaarschmidt, F.
Code was written by Budig, S.
In case of questions or comments please contact budig@cell.uni-hannover.de!


The code was written/evaluated/run in R with the following software versions:
R version 4.1.1 (2021-08-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252   
[3] LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reshape2_1.4.4  arm_1.12-2      lme4_1.1-29     Matrix_1.4-1   
 [5] plyr_1.8.7      forcats_0.5.1   stringr_1.4.0   dplyr_1.0.9    
 [9] purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.6   
[13] ggplot2_3.3.6   tidyverse_1.3.1 multcomp_1.4-19 TH.data_1.1-1  
[17] MASS_7.3-57     survival_3.3-1  mvtnorm_1.1-3  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.8.3     lubridate_1.8.0  lattice_0.20-45  zoo_1.8-10      
 [5] assertthat_0.2.1 utf8_1.2.2       R6_2.5.1         cellranger_1.1.0
 [9] backports_1.4.1  reprex_2.0.1     coda_0.19-4      httr_1.4.2      
[13] pillar_1.7.0     rlang_1.0.2      readxl_1.4.0     rstudioapi_0.13 
[17] minqa_1.2.4      nloptr_2.0.0     splines_4.1.1    munsell_0.5.0   
[21] broom_0.8.0      compiler_4.1.1   modelr_0.1.8     pkgconfig_2.0.3 
[25] tidyselect_1.1.2 codetools_0.2-18 fansi_1.0.3      crayon_1.5.1    
[29] tzdb_0.3.0       dbplyr_2.1.1     withr_2.5.0      grid_4.1.1      
[33] nlme_3.1-157     jsonlite_1.8.0   gtable_0.3.0     lifecycle_1.0.1 
[37] DBI_1.1.2        magrittr_2.0.3   scales_1.2.0     cli_3.3.0       
[41] stringi_1.7.6    fs_1.5.2         xml2_1.3.3       ellipsis_0.3.2  
[45] generics_0.1.2   vctrs_0.4.1      boot_1.3-28      sandwich_3.0-1  
[49] tools_4.1.1      glue_1.6.2       hms_1.1.1        abind_1.4-5     
[53] colorspace_2.0-3 rvest_1.0.2      haven_2.5.0     

This folder contains the following data and files that can be used to reproduce all analysis and figures of the manuscript.
It contains four subfolders containing the following files:

./Example_Analysis_NTP_bioassay/:
This folder contains all files for analysing the data example (NTP bioassay) from section 4 in the manuscript. 
The data set "miceF.rda" can be imported and analysed using the R file. 
The data set can be analysed using both a regression approach and group comparison models. 
The folder also contains the estimated values and p-values for the regression and group comparison models and the three methods.

./Figures/:
The Figures folder contains all the graphics used in the manuscript,
which can be created using the simulation code and the following data processing code and graphics code.

./Simulation_code/:
This folder contains 4 further subfolders: 

The first subfolder "./Graphs/" contains three R script files. 
The file "bin_all_graphs.R" can be used to create all the figures in the manuscript with the results of the simulations. 
The script uses the processed files from the "Simulation_results" folder "./Simulation_results/Regression_model/bin_reg_FWER_all.csv" and 
"./Simulation_results/Regression_model/bin_reg_POWER_all. csv" for generating the graphs for the regression model and 
the "./Simulation_results/Group_comparison_model/bin_gc_ss_all.csv" file for generating the graphs for the group comparison models. 

The file "./Simulation_results/Park_Generation_of_correlated_binary_data/bin_cor_endpoints.csv" is used to generate graph 7 in the appendix
The "bin_GroupComparison_data_preprocessing.R" and the "bin_Regression_data_preprocessing.R" scripts 
process the raw results of the regression and group comparison model, 
so that they can then be used with the "bin_all_graphs.R" code to generate the respective figures.

The second subfolder "./Group_comparison_model/" contains all the rscripts for the simulation of the group comparison model. 
These are a total of 43 R scripts, with each script containing one specific setting. 
Using "bin_run_groupcomparison_simulation.bat", all scripts can be called simultaneously under Windows. 
However, the correct path for the scripts folder, for the R programme and the R library must be inserted in the .bat file.
R must also be located in the system environment variables.

The third subfolder "./Park_Generation_of_correlated_binary_data/" contains the simulation 
for graphic 7 in Appendix A.1. This script can be easily executed.

The fourth subfolder "./Regression_model/" contains the rscripts for the simulation of the regression models. 
There are 6 scripts for the simulation of the FWER and 6 scripts for the simulation of the Power and FWERS. 
The scripts are divided according to model and method. 
Using the .bat file "bin_run_regression_simulation.bat" all Rscripts can be executed simultaneously under windows. 
However, the correct path for the scripts folder, for the R programme and the R library must be inserted in the .bat file.
R must also be located in the system environment variables.

./Simulation_results/:
This folder contains 3 further subfolders:

The subfolder "./Group_comparison_model/" contains all the results for the group comparison models. 
It contains both the raw results of the simulation code and the processed results which can be used to generate the graphics.

The subfolder "./Regression_model/" contains all the results for the regression models. 
It contains both the raw results of the simulation code and the processed results which can be used to generate the graphics.

The subfolder "./Park_Generation_of_correlated_binary_data/" contains the results that are used for figure 7.


Note: the simulations on a core would generally take several weeks. 
For this reason, the simulations, especially those with the bootstrap, were split across several cores. 
In the provided Rscripts, the number of simulations is set to 10000 everywhere. . 
It would therefore be advisable to reduce the number of simulations per Rscript and run several of these scripts in parallel.

