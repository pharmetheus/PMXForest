#START OF AUTO-GENERATED PREAMBLE, WILL BE OVERWRITTEN WHEN THIS FILE IS USED AS A TEMPLATE
#Created 2022-04-13 at 15:15

xpose.runno <- '7bs.mod'
toolname <- 'bootstrap'
pdf.filename <- paste0('PsN_',toolname,'_plots.pdf')
working.directory<-'/share/Projects/Pharmetheus/PMX-REP-PMX-2/Analysis/Model/SimVal/bs7n30.dir/'
results.directory <- '..'
subset.variable<-NULL
tab.suffix <- '' 
rscripts.directory <- '/opt/noarch/psn/psn-5.2.6/PsN_5_2_6/R-scripts' # This is not used
tool.results.file <- 'bootstrap_results.csv'
raw.results.file <- 'raw_results_run7bs.csv'
theta.labels <- c('1. TVFREL','2. WT coef on TVCL','3. WT coef on TVV','4. TVCL','5. TVV','6. TVMAT','7 TVD1','GENO2=1 on CL','GENO3=1 on CL','GENO4=1 on CL','CLFOOD1','FRELFORM1','FRELGENO41','FRELSEX1')
theta.fixed <- c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
omega.labels <- c('1. IIV on RUV','2. IIV on FREL','3. IIV on CL','4. IIV on V','5. IIV on MAT')
omega.fixed <- c(FALSE,FALSE,FALSE,FALSE,FALSE)
sigma.labels <- c('1. RUV')
sigma.fixed <- c(FALSE)
n.eta <- 5
n.eps <- 1

#bootstrap-specific preamble
ESTIMATED.PARAMS <- c('4. TVCL','5. TVV','6. TVMAT','7 TVD1','GENO2=1 on CL','GENO3=1 on CL','GENO4=1 on CL','CLFOOD1','FRELFORM1','FRELGENO41','FRELSEX1','1. IIV on RUV','2. IIV on FREL','3. IIV on CL','4. IIV on V','5. IIV on MAT','1. RUV')
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

