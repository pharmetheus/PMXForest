
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PMXForest

The goal of PMXForest is to make it easy to create Forest plots for
pharmacometric models, in particular NONMEM models. The user only have
to specify information about the covariate values to visualize, how the
covariates relates to the parameters of interest and where the
uncertainty information is. The latter can be a NONMEM .cov file, a PsN
bootstrap or SIR raw_results file. The Forest plots will include the
crucial information needed for an accurate interpretation of the plot
and the user have full control of the formatting and annotations.

## Installation

You can install the development version of PMXForest from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pharmetheus/PMXForest")
```

or download and install the source package.

## Example

This is a basic example which shows you how to create a Forest plot for
a NONMEM model.

First attach the PMXForest library:

``` r
library(PMXForest)
```

and set the ggplot theme

``` r
ggplot2::theme_set(ggplot2::theme_bw())
```

### Specify the covariate information

The following code sets up a data frame with information about the
coavriate combinations to visualise. It also specifies how the
covariates values should appear and how the covariate groups should be
described in the plot.

``` r
dfCovs <- createInputForestData(
  list("FOOD" = c(0,1),
       "GENO"=list("GENO1" = c(1,0,0,0),
                   "GENO3" = c(0,0,1,0),
                   "GENO4" = c(0,0,0,1)),
       "WT"   = c(65,115),
       "SEX"  = c(1,2))
  )

covnames <- c("Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM",
                         "WT 65 kg","WT 115 kg","Male","Female")

covariateGroupNames <- c("Food status","2D6 Genotype","Weight","Sex")
```

### Describe how the parameters are related to the covariates

The paramFunction specifies the relationship between the covariates and
the parameters of interest. It is not typically not necessary to
translate the whole NONMEM model to R code or specify any differential
equations. (PMXForest can be used together with mrgsolve, rxode2, etc
for more complex problems, but this is usually not necessary.) The code
below also specifies how the parameter names should apepar in the Forest
plot.

``` r
paramFunction <- function(thetas, df, ...) {
  
  if(df$WT !=-99) {
    TVCL <- thetas[4] * (df$WT / 75)**thetas[2]
  } else {
    TVCL <- thetas[4]
  }
  
  CLGENO <- 1
  if (df$GENO1 !=-99 && df$GENO1 == 1) CLGENO <- (1 + thetas[8])
  if (df$GENO3 !=-99 && df$GENO3 == 1) CLGENO <- (1 + thetas[9])
  if (df$GENO4 !=-99 && df$GENO4 == 1) CLGENO <- (1 + thetas[10])
  
  CLFOOD <- 1
  if (any(names(df)=="FOOD") && df$FOOD !=-99 && df$FOOD == 0) CLFOOD <- (1 + thetas[11])
  
  CL  <- CLFOOD*CLGENO*TVCL
  AUC <- 80/CL
  
  return(list(CL,AUC))
}

functionListName <- c("CL (L/h)","AUC (mgh/L)")
```

### Generate the data for the plot

``` r
covFile      <- "inst/extdata/SimVal/run7.cov"
extFile      <- "inst/extdata/SimVal/run7.ext"
dfSamples    <- getSamples(covFile,extFile,n=175)


dfres <- getForestDFSCM(dfCovs           = dfCovs,
                        cdfCovsNames     = covnames,
                        functionList     = list(paramFunction),
                        functionListName = functionListName,
                        noBaseThetas     = 14,
                        dfParameters     = dfSamples
)
```

### Create the Forest plot

``` r
forestPlot(dfres,groupNameLabels = covariateGroupNames,size=10)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />
