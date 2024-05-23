Supplementary information / reproducible research files for the manuscript 
Title: "Simultaneous inference of multiple binary endpoints in biomedical
research: small sample properties of multiple marginal models
and a resampling approach"

Authors: Budig, S., Jung, K., Hasler, M., Schaarschmidt, F.
Code was written by Budig, S.
In case of questions or comments please contact budig@cell.uni-hannover.de

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

This folder contains the following data and files that can be used to reproduce 
all analyses and figures of the manuscript.
It contains four subfolders with the following files:

./Example/:
This folder contains the R script "Example_Analysis_NTP_CAS NO. 93-15-2.R" for 
analysis of the example data (NTP bioassay) from Section 4 of the manuscript.
The data set can be analysed using both a regression approach and group comparison models. 
	./data/
	This folder contains the data set "miceF.rda", which can be imported and analysed 
	using the R script "Example_Analysis_NTP_CAS NO. 93-15-2.R"
	./results/
	A subfolder containing tables of the estimated effect sizes, standard errors and p-values 
	for the regression and group comparison models and the three methods. 
	The file "Table_1_Example_Analysis_NTP_Regression_pvalues.csv" is Table 1 from the manuscript.


./Figures/:
The Figures folder contains the R script ("bin_figures.R") and the corresponding 
source file ("bin_figures_source.R"), which can be used to create all the figures 
from the manuscript. For the figures of the simulation studies, the script requires 
the intermediate results, which are located in the corresponding folders of the simulations. 
There is also a batch file ("bin_run_create_figures.bat") that can be executed 
to generate the figures directly using the R script (Windows users). It must be ensured that 
all necessary packages are installed and that R is entered in the environment variables 
of the operating system. The files "Budig_Figure_1.eps" to "Budig_Figure_7.eps" are
the respective figures of the manuscript.

./Simulations/:
This folder contains 3 subfolders: 

	./Sim_Group_comparison/:
	This folder contains all the files for simulating the group comparison models. 
	The 43 R script files beginning with "bin_gc_ss" each contain a simulation scenario. 
	The R file "bin_gc_source.R" contains the functions and settings for the 
	43 simulation files. All scenarios can be simulated at once using the batch 
	file "bin_run_groupcomparison_simulation.bat". 
	The number of simulation runs and bootstrap runs can be set in the source file. 
	It is advisable to adjust the batch file to run one script per CPU thread. 
	As the simulation can take a long time depending on the CPU power, 
	it is also advisable to divide the simulation into smaller sub-simulations. 
	In the batch file, the number at the end of the line 
	"start Rscript bin_gc_ss....R number" can be used to set the suffix with which 
	the result is to be saved. The results of the simulations are saved in the 
	"intermediate_results" folder.
	
		./intermediate_results/
		This folder contains the intermediate results of the simulations of the 
		group comparison models

	./Sim_Park_cor_bin_dat/
	This folder contain the files for simulating the correlated binary data (Appendix A1).
	The R script "bin_checkPark_correlation_endpoints.R" contains the code for 
	the simulation. This script contains all functions and settings. 
	The script can be executed using the batch file "bin_run_checkPark_simulation.bat" 
	to produce the results, which are stored in the "intermediate_results" folder.
	
		./intermediate_results/
		This folder contains the result of the correlated binary data simulation.

	./Sim_Regression/
	The 12 R script files beginning with "bin_reg_FWER....R" or "bin_reg_Power....R" contain 
	the code to simulate either FWER or Power for the respective method. 
	The "bin_reg_source.R" file contains the functions and settings used 
	by the simulation files. The number of simulations and the number of bootstrap runs 
	can be set in this file. The simulations can be run simultaneously using the 
	batch file "bin_run_regression_simulation.bat". Ideally, the number of scripts 
	running simultaneously should match the number of cpu threads. 
	The simulation results are stored in the intermediate_results folder.

		./intermediate_results/
		This folder contains the intermediate results of the simulations of the 
		regression models


Note: Simulations on a core would normally take several weeks. 
For this reason, the simulations, especially those with the bootstrap, 
have been split across several cores. 
In the source files provided, the number of simulations is set to 10000. 
It would therefore be advisable to reduce the number of simulations 
per Rscript and to run several of these scripts in parallel.


