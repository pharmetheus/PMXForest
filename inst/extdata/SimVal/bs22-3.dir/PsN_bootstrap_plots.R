#START OF AUTO-GENERATED PREAMBLE, WILL BE OVERWRITTEN WHEN THIS FILE IS USED AS A TEMPLATE
#Created 2022-04-13 at 13:19

xpose.runno <- '22-3bs.mod'
toolname <- 'bootstrap'
pdf.filename <- paste0('PsN_',toolname,'_plots.pdf')
working.directory<-'/share/Projects/Pharmetheus/PMX-REP-PMX-2/Analysis/Model/SimVal/bs22-3.dir/'
results.directory <- '..'
subset.variable<-NULL
tab.suffix <- '' 
rscripts.directory <- '/opt/noarch/psn/psn-5.2.6/PsN_5_2_6/R-scripts' # This is not used
tool.results.file <- 'bootstrap_results.csv'
raw.results.file <- 'raw_results_run22-3bs.csv'
theta.labels <- c('1. TVFREL','2. WT coef on TVCL','3. WT coef on TVV','4. TVCL','5. TVV','6. TVMAT','7 TVD1','GENO2=1 on CL','GENO3=1 on CL','GENO4=1 on CL','CLFOOD1','FRELFORM1','FRELGENO41','TV_AGE','TV_AST','TV_BILI','TV_CRCL','TV_BMI','TV_HT','TV_ALT','TV_SEX','TV_RACEL2','TV_ETHNIC','TV_NCIL')
theta.fixed <- c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
omega.labels <- c('1. IIV on RUV','2. D1','2. IIV on FREL','OMEGA(4,3)','3. IIV on CL','OMEGA(5,3)','OMEGA(5,4)','4. IIV on V','OMEGA(6,3)','OMEGA(6,4)','OMEGA(6,5)','5. IIV on MAT','OMEGA(7,3)','OMEGA(7,4)','OMEGA(7,5)','OMEGA(7,6)','BSV_AGE','OMEGA(8,3)','OMEGA(8,4)','OMEGA(8,5)','OMEGA(8,6)','OMEGA(8,7)','BSV_AST','OMEGA(9,3)','OMEGA(9,4)','OMEGA(9,5)','OMEGA(9,6)','OMEGA(9,7)','OMEGA(9,8)','BSV_BILI','OMEGA(10,3)','OMEGA(10,4)','OMEGA(10,5)','OMEGA(10,6)','OMEGA(10,7)','OMEGA(10,8)','OMEGA(10,9)','BSV_CRCL','OMEGA(11,3)','OMEGA(11,4)','OMEGA(11,5)','OMEGA(11,6)','OMEGA(11,7)','OMEGA(11,8)','OMEGA(11,9)','OMEGA(11,10)','BSV_BMI','OMEGA(12,3)','OMEGA(12,4)','OMEGA(12,5)','OMEGA(12,6)','OMEGA(12,7)','OMEGA(12,8)','OMEGA(12,9)','OMEGA(12,10)','OMEGA(12,11)','BSV_HT','OMEGA(13,3)','OMEGA(13,4)','OMEGA(13,5)','OMEGA(13,6)','OMEGA(13,7)','OMEGA(13,8)','OMEGA(13,9)','OMEGA(13,10)','OMEGA(13,11)','OMEGA(13,12)','BSV_ALT','OMEGA(14,3)','OMEGA(14,4)','OMEGA(14,5)','OMEGA(14,6)','OMEGA(14,7)','OMEGA(14,8)','OMEGA(14,9)','OMEGA(14,10)','OMEGA(14,11)','OMEGA(14,12)','OMEGA(14,13)','BSV_SEX','OMEGA(15,3)','OMEGA(15,4)','OMEGA(15,5)','OMEGA(15,6)','OMEGA(15,7)','OMEGA(15,8)','OMEGA(15,9)','OMEGA(15,10)','OMEGA(15,11)','OMEGA(15,12)','OMEGA(15,13)','OMEGA(15,14)','BSV_RACEL2','OMEGA(16,3)','OMEGA(16,4)','OMEGA(16,5)','OMEGA(16,6)','OMEGA(16,7)','OMEGA(16,8)','OMEGA(16,9)','OMEGA(16,10)','OMEGA(16,11)','OMEGA(16,12)','OMEGA(16,13)','OMEGA(16,14)','OMEGA(16,15)','BSV_ETHNIC','OMEGA(17,3)','OMEGA(17,4)','OMEGA(17,5)','OMEGA(17,6)','OMEGA(17,7)','OMEGA(17,8)','OMEGA(17,9)','OMEGA(17,10)','OMEGA(17,11)','OMEGA(17,12)','OMEGA(17,13)','OMEGA(17,14)','OMEGA(17,15)','OMEGA(17,16)','BSV_NCIL')
omega.fixed <- c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
sigma.labels <- c('1. RUV','EPSCOV')
sigma.fixed <- c(FALSE,TRUE)
n.eta <- 17
n.eps <- 2

#bootstrap-specific preamble
ESTIMATED.PARAMS <- c('4. TVCL','5. TVV','6. TVMAT','7 TVD1','GENO2=1 on CL','GENO3=1 on CL','GENO4=1 on CL','CLFOOD1','FRELFORM1','FRELGENO41','TV_AGE','TV_AST','TV_BILI','TV_CRCL','TV_BMI','TV_HT','TV_ALT','TV_SEX','TV_RACEL2','TV_ETHNIC','TV_NCIL','1. IIV on RUV','2. IIV on FREL','OMEGA(4,3)','3. IIV on CL','OMEGA(5,3)','OMEGA(5,4)','4. IIV on V','OMEGA(6,3)','OMEGA(6,4)','OMEGA(6,5)','5. IIV on MAT','OMEGA(7,3)','OMEGA(7,4)','OMEGA(7,5)','OMEGA(7,6)','BSV_AGE','OMEGA(8,3)','OMEGA(8,4)','OMEGA(8,5)','OMEGA(8,6)','OMEGA(8,7)','BSV_AST','OMEGA(9,3)','OMEGA(9,4)','OMEGA(9,5)','OMEGA(9,6)','OMEGA(9,7)','OMEGA(9,8)','BSV_BILI','OMEGA(10,3)','OMEGA(10,4)','OMEGA(10,5)','OMEGA(10,6)','OMEGA(10,7)','OMEGA(10,8)','OMEGA(10,9)','BSV_CRCL','OMEGA(11,3)','OMEGA(11,4)','OMEGA(11,5)','OMEGA(11,6)','OMEGA(11,7)','OMEGA(11,8)','OMEGA(11,9)','OMEGA(11,10)','BSV_BMI','OMEGA(12,3)','OMEGA(12,4)','OMEGA(12,5)','OMEGA(12,6)','OMEGA(12,7)','OMEGA(12,8)','OMEGA(12,9)','OMEGA(12,10)','OMEGA(12,11)','BSV_HT','OMEGA(13,3)','OMEGA(13,4)','OMEGA(13,5)','OMEGA(13,6)','OMEGA(13,7)','OMEGA(13,8)','OMEGA(13,9)','OMEGA(13,10)','OMEGA(13,11)','OMEGA(13,12)','BSV_ALT','OMEGA(14,3)','OMEGA(14,4)','OMEGA(14,5)','OMEGA(14,6)','OMEGA(14,7)','OMEGA(14,8)','OMEGA(14,9)','OMEGA(14,10)','OMEGA(14,11)','OMEGA(14,12)','OMEGA(14,13)','BSV_SEX','OMEGA(15,3)','OMEGA(15,4)','OMEGA(15,5)','OMEGA(15,6)','OMEGA(15,7)','OMEGA(15,8)','OMEGA(15,9)','OMEGA(15,10)','OMEGA(15,11)','OMEGA(15,12)','OMEGA(15,13)','OMEGA(15,14)','BSV_RACEL2','OMEGA(16,3)','OMEGA(16,4)','OMEGA(16,5)','OMEGA(16,6)','OMEGA(16,7)','OMEGA(16,8)','OMEGA(16,9)','OMEGA(16,10)','OMEGA(16,11)','OMEGA(16,12)','OMEGA(16,13)','OMEGA(16,14)','OMEGA(16,15)','BSV_ETHNIC','OMEGA(17,3)','OMEGA(17,4)','OMEGA(17,5)','OMEGA(17,6)','OMEGA(17,7)','OMEGA(17,8)','OMEGA(17,9)','OMEGA(17,10)','OMEGA(17,11)','OMEGA(17,12)','OMEGA(17,13)','OMEGA(17,14)','OMEGA(17,15)','OMEGA(17,16)','BSV_NCIL','1. RUV')
included.ids.file <- 'included_individuals1.csv'
skip.minimization.terminated=TRUE
skip.covariance.step.terminated=FALSE
skip.with.covstep.warnings=FALSE
skip.estimate.near.boundary=TRUE

setwd(working.directory)

############################################################################
#END OF AUTO-GENERATED PREAMBLE
#WHEN THIS FILE IS USED AS A TEMPLATE THIS LINE MUST LOOK EXACTLY LIKE THIS


pdf(file=pdf.filename,width=10,height=7)
library(PsNR)
library(magrittr)
library(methods)
library(xpose4)
library(dplyr)
library(PerformanceAnalytics)
#add R_info to the meta file
R_info(directory=working.directory)

meta <- PsNR::metadata(working.directory)
bootplots <- xpose4::boot.hist(results.file=raw.results.file,incl.ids.file=included.ids.file,
                       min.failed=skip.minimization.terminated,
                       cov.failed=skip.covariance.step.terminated,
                       cov.warnings=skip.with.covstep.warnings,
                       boundary=skip.estimate.near.boundary)
print(bootplots[1]) #parameters
if (PsNR::rplots_level(meta) > 1){
    print(bootplots[2:4]) #SEs ofv eigenvalues
}
if(file.exists("raw_results_dofv.csv")) {
  df <- read.csv("raw_results_dofv.csv",stringsAsFactors = F)
  needed_column <- fix_column_names(c("deltaofv", ESTIMATED.PARAMS))
  if(nrow(df)>1) {
    df <- df %>% 
      dplyr::select(!!needed_column) %>%
      dplyr::slice(-1)
    suppressWarnings(PerformanceAnalytics::chart.Correlation(df, histogram = TRUE, method = c("spearman")))
  }
}
dev.off()

