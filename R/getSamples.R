#' getSamples
#'
#' @description Function that gets samples, either from a regular bootstrap (PsN
#'   output), a SIR (PsN output), or a covariance matrix (NONMEM output).
#'   It includes protective coding to ensure the parameter estimates in the
#'   .ext data align with the dimensions of the .cov matrix.
#'
#' @importFrom dplyr filter
#' @param input The input file name (character string). If the file name
#'   extension is `.cov` it is assumed that the file is a NONMEM .cov file.
#'   If the extension is `.csv`, it is assumed to be a PsN raw results file
#'   from a bootstrap or SIR. Can also be a data.frame with the same
#'   structure as such files.
#' @param extFile The name of the NONMEM .ext file (character string) OR a
#'   data.frame containing the parameter iterations. Providing a data.frame
#'   allows for manual alignment or "padding" of parameters if NONMEM has
#'   omitted columns in the .cov file.
#' @param n The number of samples to generate when sampling from a
#'   variance-covariance matrix. This is required if `input` is a `.cov` file.
#'   NULL is permissible if the input is a SIR or bootstrap raw_results file.
#' @param indexvec A vector used to subset columns in the raw results file
#'   so parameters appear in the order: THETA, SIGMA, OMEGA. Not required
#'   if a PsN `raw_results_structure` file is in the same directory.
#' @param zerosindex A vector of indices indicating which parameters are
#'   fixed to zero.
#'
#' @return A data frame with $n + 1$ rows (if $n$ is specified) or the full set
#'   of resamples. The first row contains the final parameter estimates
#'   (ITERATION -1000000000), followed by the generated or extracted samples.
#'
#' @details
#'   **Protective Coding:** When sampling from a `.cov` file, the function
#'   verifies that the number of parameters in the `extFile` matches the
#'   dimensions of the covariance matrix. If a mismatch is detected (a known
#'   NONMEM edge case), the function will `stop()` and prompt the user to
#'   provide a manually aligned data.frame to `extFile`.
#'
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
getSamples <- function(input,
                       extFile    = NULL,
                       n          = NULL,
                       indexvec   = NULL,
                       zerosindex = NULL) {

  ## Check the input argument.
  if(!is.character(input) && !is.data.frame(input)) stop("input needs to be a character string or a data frame")
  if (is.character(input)) {
    if(tools::file_ext(input) == "") stop("The input file name needs to have an extension")
    if(tools::file_ext(input) != "cov" &
       tools::file_ext(input) != "csv") stop("The input file name needs to have either a .cov or .csv extension")
    if(!file.exists(input)) stop(paste("Can not find",input))
  }

  ## Check the extFile argument (Updated to allow data.frame)
  if (!is.null(extFile) && !is.data.frame(extFile)) {
    if(tools::file_ext(extFile) != "ext") stop(paste(extFile, "does not have a .ext extension"))
    if(!file.exists(extFile)) stop(paste("Can not find",extFile))
  }

  if(class(input) != "data.frame" && tools::file_ext(input) == "csv" && is.null(extFile)) stop("Need to provide an .ext file when input is a .csv file.")

  # Load the ext data
  # This section ensures dfExt is available for all subsequent blocks
  if (!is.null(extFile)) {
    if (is.data.frame(extFile)) {
      dfExt <- extFile
    } else {
      dfExt <- getExt(extFile = extFile)
    }
  }

  # A data frame was provided
  if (is.data.frame(input)) {
    dfParameters<-input
    if (!is.null(indexvec)) dfParameters <- subset(dfParameters)[, indexvec]
    # if ext file is provided
    if (!is.null(extFile)) { #If a ext file is provided, use that to transform input
      dfExtSub <- subset(dfExt, ITERATION == "-1000000000") # Get the final parameter estimates
      if (is.null(n)) {
        dfParameters <- addMissingColumns(dfParameters, dfExtSub, zerosindex)
        dfParameters <- cbind(dfParameters, OBJ = 0)
        dfParameters        <- rbind(dfExtSub[,-1],dfParameters) # Add final parameters in the top row.
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
        dfParameters <- addMissingColumns(dfParameters, dfExtSub, zerosindex)
        dfParameters <- cbind(dfParameters, OBJ = 0)
        dfParameters        <- rbind(dfExtSub[,-1],dfParameters) # Add final parameters in the top row.
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

    dfExtSub  <- subset(dfExt, ITERATION == "-1000000000") # Get the final parameter estimates
    dfcov     <- read.table(input, fill = TRUE, header = TRUE, sep = "", skip = 1, stringsAsFactors = FALSE)
    sigma     <- data.matrix(dfcov[, 2:ncol(dfcov)])

    # Identify parameter columns in ext (excluding ITERATION and OBJ)
    strNames  <- names(dfExtSub[, -(c(1, ncol(dfExtSub)))])
    mu        <- as.numeric(dfExtSub[, -(c(1, ncol(dfExtSub)))]) # The mean is the final estimates

    ## Dimensionality Check: Protective Coding
    n_ext <- length(strNames)
    n_cov <- nrow(sigma)

    if (n_ext != n_cov) {
      stop(paste0(
        "Critical Mismatch: The .ext data contains ", n_ext, " parameters, ",
        "but the .cov matrix has ", n_cov, " parameters.\n",
        "NONMEM may have omitted parameters from the .cov file. To fix this, ",
        "provide a modified data.frame to 'extFile' that matches the .cov dimensions."
      ))
    }

    #Get the parameters which are fixed based on the covariance
    fixedmu = rep(FALSE,1,ncol(sigma))
    for (j in 1:ncol(sigma)) {
      fixedmu[j]<-all(sigma[,j]==0)
    }
    # Draw n samples from cov matrix
    dfParameters        <- as.data.frame(mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = n,fixed_mu = fixedmu))
    names(dfParameters) <- strNames
    dfParameters        <- cbind(dfParameters, OBJ = 0) # Add a dummy column with the OBJ to make it look more like a ext file
    dfParameters        <- rbind(dfExtSub[,-1],dfParameters) # Add final parameters in the top row.
    return(dfParameters)
  }

  ## Sample from a PsN raw results file if a .csv was provided
  if (tools::file_ext(input) == "csv") { # Submitted a bootstrap file, SSE or a SIR file

    ## Read the csv data
    dfParameters <- read.csv(file = input, header = TRUE,stringsAsFactors = FALSE)
    dfExtSub <- subset(dfExt, ITERATION == "-1000000000")

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
    }

    ## If its a SIR results file
    if ("resamples" %in% names(dfParameters) & "samples_order" %in% names(dfParameters)) {
      dfParameters <- subset(dfParameters, resamples == 1)[, indexvec]
      dfParameters <- addMissingColumns(dfParameters, dfExtSub, zerosindex)
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
      dfParameters <- addMissingColumns(dfParameters, dfExtSub, zerosindex)
      dfParameters <- cbind(dfParameters, OBJ = 0)
      return(dfParameters)

    } else {
      ## Remove rows with ofv = 0
      dfParameters <- dfParameters %>% filter(ofv != 0)
      dfParameters <- dfParameters[, indexvec]
      dfParameters <- addMissingColumns(dfParameters, dfExtSub, zerosindex)
      dfParameters <- cbind(dfParameters, OBJ = 0)
      return(dfParameters)
    }
  }
}
