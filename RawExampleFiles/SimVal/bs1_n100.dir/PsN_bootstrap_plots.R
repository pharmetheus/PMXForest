#START OF AUTO-GENERATED PREAMBLE, WILL BE OVERWRITTEN WHEN THIS FILE IS USED AS A TEMPLATE
#Created 2020-06-05 at 13:05

rplots.level <- 1
xpose.runno <- '1'
toolname <- 'bootstrap'
pdf.filename <- paste0('PsN_',toolname,'_plots.pdf')
pdf.title <- 'bootstrap diagnostic plots run 1'
working.directory<-'/home/shared/Projects/Pharmetheus/PMX-MI-PMX-1/Analysis/Model/FRX/bs1_n100.dir/'
model.directory<-'/home/shared/Projects/Pharmetheus/PMX-MI-PMX-1/Analysis/Model/FRX/'
results.directory <- '/home/shared/Projects/Pharmetheus/PMX-MI-PMX-1/Analysis/Model/FRX/'
model.filename<-'run1.mod'
subset.variable<-NULL
mod.suffix <- '.mod'
mod.prefix <- 'run'
tab.suffix <- ''
rscripts.directory <- '/opt/psn/psn-4.8.1/PsN_4_8_1/R-scripts'
tool.results.file <- 'bootstrap_results.csv'
raw.results.file <- 'raw_results_run1.csv'
theta.labels <- c('1. TVFREL','2. WT coef on TVCL','3. WT coef on TVV','4. TVCL','5. TVV','6. TVMAT','7 TVD1','GENO=1 on CL','GENO=3 on CL','GENO=4 on CL','GENO=1 on CL_','GENO=3 on CL_','GENO=4 on CL_','CLFOOD1','FRELFORM1','FRELNCIL1')
theta.fixed <- c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
omega.labels <- c('1. IIV on RUV','2. IIV on FREL','3. IIV on CL','4. IIV on V','5. IIV on MAT')
omega.fixed <- c(FALSE,FALSE,FALSE,FALSE,FALSE)
sigma.labels <- c('1. RUV')
sigma.fixed <- c(FALSE)
n.eta <- 5
n.eps <- 1

#bootstrap-specific preamble
ESTIMATED.PARAMS <- c('4. TVCL','5. TVV','6. TVMAT','7 TVD1','GENO=1 on CL','GENO=3 on CL','GENO=4 on CL','GENO=1 on CL_','GENO=3 on CL_','GENO=4 on CL_','CLFOOD1','FRELFORM1','FRELNCIL1','1. IIV on RUV','2. IIV on FREL','3. IIV on CL','4. IIV on V','5. IIV on MAT','1. RUV')
dofv.is.run <- 0

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
# get libPaths
source(file.path(rscripts.directory,"common/R_info.R"))
R_info(directory=working.directory,only_libPaths=T)
source(paste0(rscripts.directory, "/bootstrap/cook.cov.calcul.R"))
source(paste0(rscripts.directory, "/bootstrap/plot.cook.cov.R"))
source(paste0(rscripts.directory, "/bootstrap/format.dofv.data.R"))
source(paste0(rscripts.directory, "/bootstrap/plot.dofv.R"))
source(paste0(rscripts.directory, "/common/fix.column.names.R"))

library(xpose4)
library(ggplot2)
library(grid)
require(plyr)
require(dplyr)
#add R_info to the meta file
R_info(directory=working.directory)
est.param.names <- fix_column_names(col_names=ESTIMATED.PARAMS)
if(packageVersion("xpose4")<"4.5.0"){
		warning("xpose4 version must be 4.5.0 or later for bootstrap plot")	
}							 
bootplots <- boot.hist(results.file=raw.results.file,incl.ids.file=included.ids.file,
                       min.failed=skip.minimization.terminated,
                       cov.failed=skip.covariance.step.terminated,
                       cov.warnings=skip.with.covstep.warnings,
                       boundary=skip.estimate.near.boundary)
print(bootplots[1]) #parameters

add_dOFV_plots <- FALSE
if (dofv.is.run){
  add_dOFV_plots <- TRUE
  
  # Read in and format data
  list_dofv <- format_dofv_data(dofv.raw.results.file,raw.results.file,est.param.names)
  all <- list_dofv$all
  df_est <- list_dofv$df_est

  # Plot dOFV distributions
  qdOFV_all <- plot_dofv(all,df_est,est.param.names)
  print(qdOFV_all)
}

if (rplots.level > 1){
    print(bootplots[2:4]) #SEs ofv eigenvalues
}
#calculate cook scores and cov ratios
list_cook.cov <- cook_cov_calcul(raw.results.file,included.ids.file,est.param.names)
#unlist
data_plots <- list_cook.cov$data_plots
failed_cov_ID <- list_cook.cov$failed_cov_ID
estimation_failures <- list_cook.cov$estimation_failures
samples <- list_cook.cov$samples

#plot cook scores and cov ratios to find influential individuals
if (rplots.level > 1) {
    gt <- plot_cook_cov(data_plots,failed_cov_ID,samples,estimation_failures)
    #plot
    grid.draw(gt)
}
dev.off()

