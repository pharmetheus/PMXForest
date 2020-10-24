#' mvrnorm_vector
#'
#' @description Samples vectors of values from a vector of means and a varcov matrix
#' @param mu Vector of means
#' @param sigma Covariance matrix
#' @param fixed_mu Indicate that some parameters are fixed, i.e. guarantee that the exact
#' mu_value is returned for fixed parameters, especially important for zero parameters
#' which might otherwise return samples in the magnitude of <1E-16 instead of exactly 0
#' @param dSeed The seed number
#' @param iSampleIndex Number of samples to draw
#' @return A vector if iSampleIndex==1 otherwise a matrix with iSampleIndex rows
#' @export
#'
#' @examples
#' \dontrun{
#' ## Read the ext-file
#' extFile <- getExt(system.file("extdata","run1.ext",package="PMXForest"))
#'
#' ## Read the cov-file and extract the covmatrix
#' dfcov     <- read.table(system.file("extdata","run1.cov",package="PMXForest"),
#'                         fill = TRUE, header = TRUE, sep = "",
#'                                                 skip = 1, stringsAsFactors = FALSE)
#' sigma     <- data.matrix(dfcov[, 2:ncol(dfcov)])
#'
#' ## Get the final parameter estimates
#' finPar  <- subset(extFile, ITERATION == "-1000000000")
#' mu      <- as.numeric(finPar[, -(c(1, ncol(finPar)))])
#'
#' dfParameters <- as.data.frame(mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = 10))
#' }
mvrnorm_vector <- function(mu,sigma,fixed_mu=NULL,dSeed=NULL,iSampleIndex=1)
{
  if (iSampleIndex==0) return (mu)

  if (is.null(fixed_mu)) fixed_mu <- rep(0,length(mu))

  tmp_mu <- mu[which(fixed_mu==0)] #Get the non-fixed mu

  if (!is.null(dSeed)) set.seed(dSeed)

  ## Get the sigma
  if(!any(fixed_mu!=0)) {
    mySigma <- sigma
  } else {
    mySigma <- sigma[-which(fixed_mu==1),-which(fixed_mu==1)]
  }

  samples <- MASS::mvrnorm(n=iSampleIndex,tmp_mu,Sigma = mySigma)

  #samples <- mvrnorm(n=iSampleIndex,tmp_mu,sigma[-which(fixed_mu==1),-which(fixed_mu==1)])

  if (!is.matrix(samples)) {
    mu[which(fixed_mu==0)]<-samples
  } else {
    mum<-matrix(0,ncol=length(mu),nrow=iSampleIndex) #Create a matrix for the samples
    for (i in 1:iSampleIndex){
      mum[i,]<-mu
      mum[i,which(fixed_mu==0)]<-samples[i,]
    }
    return(mum)
  }
  return (mu)
}

