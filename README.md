# Anthropogenic habitat modification alters calling phenology of frogs

This repository contains data and code to reproduce Liu et al. (2022) Anthropogenic habitat modification alters calling phenology of frogs. Global Change Biology. The complete raw dataset is not fully available due to sensitivities in relation to locations of rare or threatened species and citizen scientist information (Rowley and Callaghan, 2020), but the FrogID data, with sensitive species’ localities removed or buffered, are made available annually (Rowley and Callaghan, 2020; https://doi.org/10.3897/zookeys.912.38253; data available at GBIF: https://doi.org/10.15468/wazqft and FrogID: https://www.frogid.net.au/explore). Here we provide processed breeding phenology estimates (duration, onset, median, and offset of the breeding season) and call characteristic data, and the code to reproduce the main analyses and figures in the manuscript. The following R information was used at the time of analysis:

sessionInfo()

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22000)

Matrix products: default

locale:
[1] LC_COLLATE=English_Australia.1252  LC_CTYPE=English_Australia.1252   
[3] LC_MONETARY=English_Australia.1252 LC_NUMERIC=C                      
[5] LC_TIME=English_Australia.1252    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gridExtra_2.3   cowplot_1.1.1   sf_1.0-7        mgcv_1.8-38     nlme_3.1-153   
 [6] ggeffects_1.1.2 lmerTest_3.1-3  lme4_1.1-29     Matrix_1.3-4    forcats_0.5.1  
[11] stringr_1.4.0   dplyr_1.0.8     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0    
[16] tibble_3.1.6    ggplot2_3.3.5   tidyverse_1.3.1

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.8          lubridate_1.8.0     lattice_0.20-45     class_7.3-19       
 [5] assertthat_0.2.1    utf8_1.2.2          R6_2.5.1            cellranger_1.1.0   
 [9] backports_1.4.1     reprex_2.0.1        e1071_1.7-9         httr_1.4.2         
[13] pillar_1.7.0        rlang_1.0.1         readxl_1.3.1        rstudioapi_0.13    
[17] minqa_1.2.4         nloptr_2.0.2        splines_4.1.2       munsell_0.5.0      
[21] proxy_0.4-26        broom_0.7.12        compiler_4.1.2      numDeriv_2016.8-1.1
[25] modelr_0.1.8        pkgconfig_2.0.3     tidyselect_1.1.2    fansi_1.0.2        
[29] crayon_1.5.0        tzdb_0.2.0          dbplyr_2.1.1        withr_2.5.0        
[33] MASS_7.3-54         jsonlite_1.8.0      gtable_0.3.0        lifecycle_1.0.1    
[37] DBI_1.1.2           magrittr_2.0.2      units_0.8-0         scales_1.1.1       
[41] KernSmooth_2.23-20  cli_3.2.0           stringi_1.7.6       fs_1.5.2           
[45] xml2_1.3.3          ellipsis_0.3.2      generics_0.1.2      vctrs_0.3.8        
[49] boot_1.3-28         tools_4.1.2         glue_1.6.2          hms_1.1.1          
[53] colorspace_2.0-3    sessioninfo_1.2.2   classInt_0.4-3      rvest_1.0.2        
[57] haven_2.4.3
