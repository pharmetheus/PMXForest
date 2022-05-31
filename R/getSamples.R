#' getSamples
#'
#' @description Function that gets samples, either from a regular bootstrap (Psn output) a SIR (PsN output), a covariance matrix (NM output)
#' or samples from a small bootstrap using the estimated covariance matrix
#'
#' @importFrom dplyr filter
#' @param input The input file name (character string).
#' bootstrap or SIR. If the file name extension is ".cov" it is assumed that the file
#' is a a NONMEM .cov file and if the extension is '.csv' is is assumed that it is a PsN raw results file from either a bootstrap or a SIR.
#' @param extFile The name of the NONMEM ext file. Needed to transform output to the ext-format.  Should not be NULL.
#' @param n The number of samples to generate. NULL is permissable if the input is a SIR raw_results file.
#' @param indexvec A vector that can be used to subset the columns in the raw results file so that the columns comes in the order THETA, SIGMA, OMEGA.
#' not required if a PsN raw_results_structure is present in the same directory as the csv file.
#' @param zerosindex A vector with indicies that indicates which parameters that are fixed to zero.
#'
#' @return A data frame with parameter vectors.
#'
#' In case a .cov file is provided as input the parameter vectors is the result of multivariate sampling from the variance-covariance matrix.
#'
#' In case a .csv file is provided as input and n=NULL, the input is assumed to be a SIR or bootstrap file. In case of a SIR file, all
#' samples present in the .csv file is returned. In case it is a bootstrap file, the function returns all samples with a nonzero value in the
#' ofv column.
#'
#' In case a .csv file is provided as input and n is not NULL, it is assumed that the raw results come from a boostrap analysis. A variance-
#' covariance matrix will be created from the samples with a non-zero value in the ofv column and
#' n parameter vectors will be returned obtained by multivariate sampling from the constructed variance-covariance matrix.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Covariance matrix
#' getSamples(covFile,extFile,n=175)
#'
#' # Bootstrap
#' getSamples(BSFile,extFile)
#'
#' # SIRFile
#' getSamples(SIRFile,extFile)
#'
#' # Small bootstrap mvnorm sampling
#' getSamples(BSFile,extFile,n=175)
#'
#' # Data frame
#' dfParameters<-read.csv(rawresFile)
#'
#' getSamples(dfParameters)
#'
#' # Data frame mvnorm sampling
#' getSamples(dfParameters,n=175)
#' }
getSamples <- function(input,extFile=NULL,n=NULL,indexvec=NULL,zerosindex=NULL) {

  ## Check the input argument.
  if(!is.character(input) && !is.data.frame(input)) stop("input needs to be a character string or a data frame")
  if (is.character(input)) {
    if(tools::file_ext(input) == "") stop("The input file name needs to have an extension")
    if(tools::file_ext(input) != "cov" &
       tools::file_ext(input) != "csv") stop("The input file name needs to have either a .cov or .csv extension")
    if(!file.exists(input)) stop(paste("Can not find",input))
  }
  ## Check the extFile argument
  if(!is.null(extFile) && tools::file_ext(extFile) != "ext") stop(paste(extFile, "does not have a .ext extension"))
  if(!is.null(extFile) && !file.exists(extFile)) stop(paste("Can not find",extFile))

  # A data frame was provided
  if (is.data.frame(input)) {
    dfParameters<-input
    if (!is.null(indexvec)) dfParameters <- subset(dfParameters)[, indexvec]
    # if ext file is provided
    if (!is.null(extFile)) { #If a ext file is provided, use that to transform input
      dfExt     <- subset(getExt(extFile = extFile), ITERATION == "-1000000000") # Get the final parameter estimates
      if (is.null(n)) {
        dfParameters <- addMissingColumns(dfParameters, dfExt, zerosindex)
        dfParameters <- cbind(dfParameters, OBJ = 0)
        dfParameters        <- rbind(dfExt[,-1],dfParameters) # Add final parameters in the top row.
        return(dfParameters)
      } else {
        # Now construct sigma and mu
        sigma <- cov(dfParameters)
        mu    <- colMeans(dfParameters)

        #Get the parameters which are fixed based on the covariance
        fixedmu = rep(FALSE,1,ncol(sigma))
        for (j in 1:ncol(sigma)) {
          fixedmu[j]<-all(sigma[,j]==0)
        }

        dfParameters <- as.data.frame(mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = n,fixed_mu=fixedmu))
        dfParameters <- addMissingColumns(dfParameters, dfExt, zerosindex)
        dfParameters <- cbind(dfParameters, OBJ = 0)
        dfParameters        <- rbind(dfExt[,-1],dfParameters) # Add final parameters in the top row.
        return(dfParameters)
      }
    } else {
      if (is.null(n)) return(dfParameters)
      # Now construct sigma and mu
      sigma <- cov(dfParameters)
      mu    <- colMeans(dfParameters)

      #Get the parameters which are fixed based on the covariance
      fixedmu = rep(FALSE,1,ncol(sigma))
      for (j in 1:ncol(sigma)) {
        fixedmu[j]<-all(sigma[,j]==0)
      }
      tmpnames<-names(dfParameters)
      dfParameters <- as.data.frame(mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = n,fixed_mu=fixedmu))
      names(dfParameters)<-tmpnames
      dfParameters <- cbind(dfParameters, OBJ = 0)
      return(dfParameters)
  }
 }
  ## Sample from the varcov matrix if a .cov was provided
  if(tools::file_ext(input) == "cov") {

    ## Check that n is provided
    if(is.null(n)) stop("n - number of samples - needs to be provided")

    dfExt     <- subset(getExt(extFile = extFile), ITERATION == "-1000000000") # Get the final parameter estimates
    dfcov     <- read.table(input, fill = TRUE, header = TRUE, sep = "", skip = 1, stringsAsFactors = FALSE)
    sigma     <- data.matrix(dfcov[, 2:ncol(dfcov)])
    strNames  <- names(dfExt[, -(c(1, ncol(dfExt)))])
    mu        <- as.numeric(dfExt[, -(c(1, ncol(dfExt)))]) # The mean is the final estimates

    #Get the parameters which are fixed based on the covariance
    fixedmu = rep(FALSE,1,ncol(sigma))
    for (j in 1:ncol(sigma)) {
      fixedmu[j]<-all(sigma[,j]==0)
    }
    # Draw n samples from cov matrix
    dfParameters        <- as.data.frame(mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = n,fixed_mu = fixedmu))
    names(dfParameters) <- strNames
    dfParameters        <- cbind(dfParameters, OBJ = 0) # Add a dummy column with the OBJ to make it look more like a ext file
    dfParameters        <- rbind(dfExt[,-1],dfParameters) # Add final parameters in the top row.
    return(dfParameters)
  }





  ## Sample from a PsN raw results file if a .csv was provided
 # args <- list(...)

  if (tools::file_ext(input) == "csv") { # Submitted a bootstrap file, SSE or a SIR file
    ## Figure out the parameter positions in the input file
    if(!is.null(indexvec)) {
      indexvec <- as.numeric(indexvec)
    } else {
      rr_struct <- file.path(dirname(input), "raw_results_structure")

      ## Check that the raw_results structure file exists
      if(!file.exists(rr_struct)) stop(paste(rr_struct,"does not exist and indexvec is not provided"))

      res <- getParamPositions(rr_struct)
      if (!is.null(res)) {
        s1 <- s2 <- s3 <- NULL
        if (length(res$THETA) > 0 && res$THETA[2]>0) s1 <- seq(res$THETA[1], res$THETA[1] + res$THETA[2] - 1)
        if (length(res$SIGMA) > 0 && res$SIGMA[2]>0) s2 <- seq(res$SIGMA[1], res$SIGMA[1] + res$SIGMA[2] - 1)
        if (length(res$OMEGA) > 0 && res$OMEGA[2]>0) s3 <- seq(res$OMEGA[1], res$OMEGA[1] + res$OMEGA[2] - 1)
        indexvec <- c(s1, s2, s3) + 1
      }
    }

    ## Read the csv and ext files
    dfParameters <- read.csv(file = input, header = TRUE,stringsAsFactors = FALSE)

    ## Add SIGMA and OMEGA last in dfParameters if they are completely missing
    if (is.null(s2) || is.null(s3)) {
      if (is.null(s2)) {
        dfParameters["SIGMA.1.1."]<-0
        s2<-ncol(dfParameters)-1
      }
      if (is.null(s3)) {
        dfParameters["OMEGA.1.1."]<-0
        s3<-ncol(dfParameters)-1
      }
      indexvec <- c(s1,s2,s3) + 1
    }

    dfExt <- subset(getExt(extFile = extFile), ITERATION == "-1000000000")

    ## If its a SIR results file
    if ("resamples" %in% names(dfParameters) & "samples_order" %in% names(dfParameters)) {
      dfParameters <- subset(dfParameters, resamples == 1)[, indexvec]
      dfParameters <- addMissingColumns(dfParameters, dfExt, zerosindex)
      dfParameters <- cbind(dfParameters, OBJ = 0)
      return(dfParameters)
    }

    ## Regular raw_results from a bootsrap
    if (!is.null(n)) { # Also provided n, i.e. construct a cov matrix
      dfParameters <- dfParameters %>% filter(ofv != 0)

      dfParameters <- dfParameters[, indexvec]

      # Now construct sigma and mu
      sigma <- cov(dfParameters)
      mu    <- colMeans(dfParameters)

      #Get the parameters which are fixed based on the covariance
      fixedmu = rep(FALSE,1,ncol(sigma))
      for (j in 1:ncol(sigma)) {
        fixedmu[j]<-all(sigma[,j]==0)
      }

      dfParameters <- as.data.frame(mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = n,fixed_mu=fixedmu))
      dfParameters <- addMissingColumns(dfParameters, dfExt, zerosindex)
      dfParameters <- cbind(dfParameters, OBJ = 0)
      return(dfParameters)

    } else {
      ## Remove rows with ofv = 0
      dfParameters <- dfParameters %>% filter(ofv != 0)
      dfParameters <- dfParameters[, indexvec]
      dfParameters <- addMissingColumns(dfParameters, dfExt, zerosindex)
      dfParameters <- cbind(dfParameters, OBJ = 0)
      return(dfParameters)
    }
  }

}

