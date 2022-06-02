# AMP_NeutTiterBiomarker
Code implementing the methods to generate the figures and tables of the AMP Neutralization Titer Biomarker manuscript 

1. System requirements

  Unless otherwise specified software was tested on Ubuntu 18.04.5 LTS
  running R version 3.6.3 (2020-02-29)
  
  Specific R libraries and R files required by each program are listed below.

  gen_fig1_extFig2.R: 
    magrittr version 1.5
    dplyr version 0.8.5
    CCC_time.R

  gen_fig2_extFig4.R
    sievePH version 1.0.1
    plot.summary.sievePH.R
    
  gen_fig3.4.5.6_extFig5.6.7.8.9.10.R
    dplyr version 1.0.4
    ggplot2 version 3.3.3
    tidyr version 1.1.2
    gridExtra version 2.3
    triple-bnab-BH-titer.R
    tested on Windows 10.0.19042 Build 19042 running R version 4.0.2 (2020-06-22)

  gen_suppTable1.R
    survival version 3.1.12
    cmprsk version 2.2.9
    dplyr version 0.8.5

  gen_suppTable2.R
    dplyr version 0.8.5
    haven version 2.2.0

2. Installation guide
  
  Install required version of R. 
  Install required R packages.
  Create a local directory for code and data
  Unzip Supplemental Software file to a local directory (should take less than 2 minutes).
  Clone this repository (should take less than 2 minutes). 
  
3. Demo

  Each code file should run in under 2 minutes with the exception of gen_fig3.4.5.6_extFig5.6.7.8.9.10.R which 
  takes approximately 2 days to run.  Output corresponds to figure(s) (as pdf(s) file in the figures directory) or 
  table (as output in gen_suppTable1.Rout). 
  
  
  From the command line starting in the code directory run the following commands:

    R CMD BATCH gen_fig1_extFig2.R &
    R CMD BATCH gen_fig2_extFig4.R &
    R CMD BATCH gen_fig3.4.5.6_extFig5.6.7.8.9.10.R &
    R CMD BATCH gen_suppTable1.R &

4. Instructions for use

  Same as demo instructions above.
  
