% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createInputForestData.R
\name{createInputForestData}
\alias{createInputForestData}
\title{createInputForestData}
\usage{
createInputForestData(listCovs, iMiss = -99)
}
\arguments{
\item{listCovs}{a list with covariate names and values of the form; list("AGE"=c(1,2)) etc.
For multiple covariates (i.e. not univariate) a list of lists should be used, e.g. list("AGEBW"=list("AGE"=c(1,2),"BW"=c(3,4)))}

\item{iMiss}{The value that should be filled for missing, default = -99}
}
\value{
a data frame with a covariate value (univariate) per row or several covariate values per row (multivariate). The data frame also contains a COVARIATEGROUPS column based on the names in the input list.
}
\description{
createInputForestData
}
\examples{
\dontrun{
## One covariate effect per 'group' on the Forest plot
dfCovs <- createInputForestData(
list("FORM" = c(0,1),
     "FOOD" = c(0,1),
     "GENO" = c(1,2,3,4),
     "RACEL"= c(1,2,3),
     "WT"   = c(65,115),
     "AGE"  = c(27,62),
     "CRCL" = c(83,150),
     "SEX"  = c(1,2)
     ),
     Miss=-99)
}

## In this case the GENO covariate levels (2,3 and 4) will be grouped in the Forest plot.
dfCovs <- createInputForestData(
list("FORM" = c(0,1),
     "FOOD" = c(0,1),
     "GENO"=list("GENO_2" = c(0,1,0,0),
                 "GENO_3" = c(0,0,1,0),
                 "GENO_4" = c(0,0,0,1))
          ),
          iMiss=-99
)

}
