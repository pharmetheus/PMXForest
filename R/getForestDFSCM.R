#' getForestDFSCM
#'
#' @description Get a data frame with Forest border for each univariate or multivariate covariate (and value(s)) in the input data frame. If a list a data frame will be created from the list, see function dfCreateInputForestData
#'
#' @import doParallel
#' @import foreach
#' @import dplyr
#'
#' @param dfCovs A data frame with covariates to include, if a covariate value is set -99 or NA,
#' they are assumed missing and will not be included in any FFEM transformations. If a dfCovs is
#' a list an attempt will be made to create the appropriate data frame with the createInputForestData function.
#' @param cdfCovsNames A string vector with names of the rows in dfCovs, if not used, names will be
#' automatically assigned based on the covariate values and column names in dfCovs.
#' @param functionList A list of functions with input (basethetas, covthetas,dfrow and ...) from which the change
#' from the reference value will be calculated. If the function returns a vector of values, each value will be used
#' but functionListName must contain the names with a length of all return for all functions in the functionList
#' @param functionListName A vector of strings (names) of the parameters for each function in the functionList
#' @param noBaseThetas the number of parameters from `dfParameters` to pass to the functions in `functionList`.
#' @param dfParameters A data frame with parameter samples from the uncertainty distribution.
#' The vector of final parameter estimates is assumed to be in the first row.
#' The column order is assumed the same as in the NONMEM ext file except the ITERATION and OBJ columns whichshould not be included.
#' @param quiet If output should be allowed during the function call, default= TRUE. (This option is mainly for debugging purposes.)
#' @param probs A vector of probabilities that should be computed for each of the parameters from functionList. These will be used as the
#' as the uncertainties in the Forest plots. The probs vector position one and two will be used for plotting the uncertanties (i.e. columns q1 and q2). Default is c(0.05, 0.95).
#' @param pointFunction The function used to calculate the point for each covariate in the forest plot. default=median
#' This function is also used for the reference covariate combination
#' @param dfRefRow A data frame  (one row or equal number of rows as dfCovs) with the covariate values that will be used as the reference, if NULL the typical subject is used as reference.
#' @param cGrouping A vector of numbers defining how to group the y-axis of the Forest plot, the length of the vector should match the number of rows in dfCovs.
#' If NULL (default) an educated guess of the grouping will be set
#' @param ncores the number of cores to use for the calculations, default = 1 which means no parallellization
#' @param cstrPackages a character vector with package names needed to run the calculations in parallel, default = NULL
#' @param cstrExports a character vector with variables needed to run the calculations in parallel, default = NULL
#' @param iMiss The missing value number. -99 by default.
#' @param oneHot An optional one-hot encoding specification (see
#' \code{\link{oneHotEncode}}). When supplied, raw multi-level categorical columns
#' are one-hot encoded before the parameter functions are evaluated - in `dfCovs`
#' and `dfRefRow` for `getForestDFSCM()`, in `dfData` for `getForestDFemp()` - so
#' a single parameter function written against the dummy columns can serve both
#' the parametric and empirical workflows. Default `NULL` (no encoding). For
#' `getForestDFSCM()` the raw column is dropped after encoding.
#' @param oneHotSep The separator between the covariate name and the level in the
#' dummy column names created when `oneHot` is used. Defaults to `"_"`. Set to
#' `""` to match dummy columns named like `GENO1`, `GENO2`.
#' @param ... additional variables to be forwarded to the the functionList functions
#'
#'
#' @return A data frame with summary statistics for each parameters and covariate combinations:
#'
#' @export
#'
#' @examples
#' dfData  <- read.csv(
#'   system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
#' )
#' extFile <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")
#' covFile <- system.file("extdata", "SimVal/run7.cov", package = "PMXForest")
#'
#' # Covariate rows and uncertainty samples (run7 has 14 THETAs)
#' dfCovs    <- setupDfCovs(dfData, covariates = c("WT", "AGE", "CRCL", "FOOD"),
#'                          idVar = "ID")
#' dfSamples <- getSamples(covFile, extFile, n = 50)
#'
#' # Rebuild the part of the model needed for CL
#' paramFunction <- function(thetas, df, ...) {
#'   TVCL <- thetas[4]
#'   if (any(names(df) == "WT") && df$WT != -99) {
#'     TVCL <- thetas[4] * (df$WT / 75)^thetas[2]
#'   }
#'   if (any(names(df) == "FOOD") && df$FOOD != -99 && df$FOOD == 0) {
#'     TVCL <- TVCL * (1 + thetas[11])
#'   }
#'   list(CL = TVCL)
#' }
#'
#' dfres <- getForestDFSCM(
#'   dfCovs           = dfCovs,
#'   functionList     = list(paramFunction),
#'   functionListName = "CL",
#'   noBaseThetas     = 14,
#'   dfParameters     = dfSamples
#' )
#' head(dfres[, c("COVNAME", "GROUPNAME", "PARAMETER", "POINT", "Q1", "Q2")])
getForestDFSCM <- function(dfCovs,
                           cdfCovsNames = NULL,
                           functionList = list(
                             function(basethetas, dfrow, ...) {
                               return(basethetas[1])
                             }
                           ),
                           functionListName = "PAR1",
                           noBaseThetas,
                           dfParameters,
                           quiet = TRUE,
                           probs = c(0.05, 0.95),
                           pointFunction = median,
                           dfRefRow = NULL,
                           cGrouping = NULL,
                           ncores = 1,
                           cstrPackages = NULL,
                           cstrExports = NULL,
                           iMiss = -99,
                           oneHot = NULL,
                           oneHotSep = "_",
                           ...) {

  if (!is.null(dfRefRow) && nrow(dfRefRow)!=1 && nrow(dfRefRow)!=nrow(dfCovs)) {
    stop("The number of reference rows (dfRefRow) should be either NULL (missing used as reference), one (this row used as reference) or equal to dfCovs (change reference for each covariate combination)")
  }

  ## Remove samples with problems. Will use THETA1 == NA as an indicator for a problematic sample
  dfParameters <- dfParameters %>% filter(!is.na(THETA1))

  ## Try to make dfCovs into a data.frame if it isn't that already
  if (!is.data.frame(dfCovs)) {
    dfCovs <- createInputForestData(dfCovs)
  }
  ## Coerce to a plain data.frame: the internal `[` code assumes base-R
  ## drop-to-vector behaviour, which a tibble does not provide.
  dfCovs <- as.data.frame(dfCovs)
  if (!is.null(dfRefRow) && is.data.frame(dfRefRow)) dfRefRow <- as.data.frame(dfRefRow)
  dfCovs[is.na(dfCovs)] <- iMiss

  ## Optionally one-hot encode raw multi-level categorical columns so a single
  ## parameter function can be written against the dummy columns.
  if (!is.null(oneHot)) {
    dfCovs <- oneHotEncode(dfCovs, spec = oneHot, sep = oneHotSep,
                           missVal = iMiss, dropOriginal = TRUE)
    if (!is.null(dfRefRow) && is.data.frame(dfRefRow)) {
      ## A reference row built by setupDfRefRow() may already carry the dummy
      ## columns and not the raw one; encoding is then a no-op bar the warning.
      dfRefRow <- suppressWarnings(
        oneHotEncode(dfRefRow, spec = oneHot, sep = oneHotSep,
                     missVal = iMiss, dropOriginal = TRUE)
      )
    }
  }

  groupnames<-NULL #Store the temp groupnames
  if (any(names(dfCovs)=="COVARIATEGROUPS")) {
    groupnames<-dfCovs[,"COVARIATEGROUPS"]
    dfCovs[,"COVARIATEGROUPS"]<-NULL
  }

  ## Define the grouping
  getGroups <- function(df) {
    cGroups <- c()
    cUnique <- c()
    iGroup <- 0
    for (i in 1:nrow(df)) {
      tmp <- paste0(names(dfCovs[i, , drop = FALSE])[as.numeric(dfCovs[i, , drop = FALSE]) != iMiss], collapse = ",")
      if (tmp %in% cUnique) {
        tmpl <- which(tmp == cUnique)
        cGroups <- c(cGroups, tmpl)
      } else {
        iGroup <- iGroup + 1
        cGroups <- c(cGroups, iGroup)
        cUnique <- c(cUnique, tmp)
      }
    }
    return(cGroups)
  }
  if (is.null(cGrouping)) cGrouping <- getGroups(dfCovs)


  ## Register to allow for paralell computing. Tear the cluster down on exit
  ## (including on error) rather than only at the end of a successful run.
  if (ncores>1) {
    registerDoParallel(cores = ncores)
    on.exit(stopImplicitCluster(), add = TRUE)
  }

  ## Calculate the parameters
  internalCalc<-function(k) {
    thetas <- as.numeric(dfParameters[k, 1:noBaseThetas])
    dfrest <- data.frame()

    for (i in 1:nrow(dfCovs)) {
      n <- 1
      for (j in 1:length(functionList)) {
        val <- functionList[[j]](thetas = thetas, df = dfCovs[i, ,drop=FALSE], ...)
        if (!is.null(dfRefRow)) {
          indi<-min(i,nrow(dfRefRow))
          valbase <- functionList[[j]](thetas = thetas, df = dfRefRow[indi,,drop=FALSE], ...)
        }
        else {
          dfMissing <- as.data.frame(dfCovs[1,,drop=FALSE])
          dfMissing[,] <- iMiss
          valbase <- functionList[[j]](thetas = thetas, df = dfMissing, ...)
        }
        listcount <- length(val)
        for (l in 1:listcount) {
          dfrest <- bind_rows(dfrest, data.frame(
            ITER = k,
            COVS = i, NAME = functionListName[n], VALUE = val[[l]],
            VALUEBASE = valbase[[l]],
            stringsAsFactors = FALSE
          ))
          n <- n + 1
        }
      }
    }
    return(dfrest)
  }

  if (ncores>1) {
  dfres <- foreach(
    k = 1:nrow(dfParameters), .packages = cstrPackages,
    ## Bundle the whole local environment for PSOCK workers (Windows). foreach's
    ## static global detection does not reliably follow `internalCalc`'s free
    ## variables; `cstrExports` remains available for anything outside this frame.
    .export = c(ls(environment()), cstrExports),
    .verbose = !quiet, .combine = bind_rows
  ) %dopar% {
    internalCalc(k)
    }
  } else {
    dfres<-data.frame()
    for (k in 1:nrow(dfParameters)) dfres<-bind_rows(dfres,internalCalc(k))
  }

  getCovNameString <- function(dfrow) {
    strName <- ""
    colnames <- names(dfrow)
    for (i in 1:ncol(dfrow)) {
      if (dfrow[1, i] != iMiss) {
        if (strName == "") {
          strName <- paste0(colnames[i], "=", dfrow[1, i])
        } else {
          strName <- paste0(strName, ", ", colnames[i], "=", dfrow[1, i])
        }
      }
    }
    if (strName == "") {
      strName <- "Ref cov"
    }
    return(strName)
  }


  dfret <- data.frame()
  for (i in 1:nrow(dfCovs)) {
    if (is.null(cdfCovsNames)) {
      covname <- getCovNameString(dfCovs[i, , drop = FALSE])
    }
    else {
      covname <- cdfCovsNames[i]
    }
    group <- cGrouping[i]
    for (j in 1:length(functionListName)) {
      dft <- dfres[dfres$COVS == i & dfres$NAME == functionListName[j], ]

      quant <- quantile(dft$VALUE,probs = probs, names = FALSE,na.rm = T)
      #Calculate the point value of the forest plot
      FUNCVAL=pointFunction(dft$VALUE)

      #Define the relative without parameter uncertainty
      dft$RELINTERNAL <- dft$VALUE/dft$VALUEBASE
      #Get quantile and pointvalue when uncertainty is not taken into account
      quantrel <- quantile(dft$RELINTERNAL,probs = probs, names = FALSE,na.rm = T)
      FUNCNOVAR=pointFunction(dft$RELINTERNAL)

      #Calculate reference value based one the pointFunction
      func_base<-pointFunction(dft$VALUEBASE)
      true_base <- dft$VALUEBASE[dft$ITER == 1]

      #Quantiles of the value relative to the reference. Computed on the ratio
      #directly (rather than dividing `quant` by the scalar base) so the endpoints
      #stay ordered as lower/upper even when the reference value is negative.
      #Fall back to quant/base when the base is NA so an all-NA reference yields
      #NA columns rather than erroring in quantile().
      quant_reffunc  <- if (is.na(func_base))  quant/func_base  else quantile(dft$VALUE/func_base,  probs = probs, names = FALSE, na.rm = T)
      quant_reffinal <- if (is.na(true_base)) quant/true_base else quantile(dft$VALUE/true_base, probs = probs, names = FALSE, na.rm = T)
      groupname<-group
      if (!is.null(groupnames)) groupname<-groupnames[i]
      dfrow <- cbind(dfCovs[i, , drop = FALSE], data.frame(
        GROUP = group,
        GROUPNAME = groupname,
        COVNUM = i, COVNAME = covname, PARAMETER = functionListName[j],
        REFFUNC = func_base, REFFINAL = true_base, POINT=FUNCVAL, POINT_NOVAR_REL_REFFUNC=FUNCNOVAR,
        POINT_REL_REFFUNC = FUNCVAL/func_base, POINT_REL_REFFINAL = FUNCVAL/true_base,
        COVEFF = !all(dft$RELINTERNAL==1)))
      for (k in 1:length(probs)) {
        dfp <- data.frame(X1 = 1)
        dfp[[paste0("Q", k)]] <- quant[k]
        dfp[[paste0("Q",k,"_REL_REFFUNC")]] <-  quant_reffunc[k]
        dfp[[paste0("Q",k,"_REL_REFFINAL")]] <-  quant_reffinal[k]
        dfrow <- cbind(dfrow, dfp[, 2:4])
      }
      for (k in 1:length(probs)) {
        dfp <- data.frame(X1 = 1)
        dfp[[paste0("Q",k,"_NOVAR_REL_REFFUNC")]] <- quantrel[k]
        dfrow <- cbind(dfrow, dfp[, 2])
        names(dfrow)[ncol(dfrow)] <- paste0("Q", k,"_NOVAR_REL_REFFUNC")
      }
      dfret <- rbind(dfret, dfrow)
    }
  }

  ## (cluster teardown is handled by on.exit() registered above)

  ## Add a column with YES/NO depending on if refRow was provided or if it was set to the default NULL
  dfret <- dfret %>% mutate(REFROW = ifelse(is.null(dfRefRow),"NO","YES"))

  ## Make sure GROUPNAME, COVNAME, PARAMETER are factors
  dfret$GROUPNAME <- factor(dfret$GROUPNAME,levels=unique(dfret$GROUPNAME))
  dfret$COVNAME   <- factor(dfret$COVNAME,levels=unique(dfret$COVNAME))
  dfret$PARAMETER <- factor(dfret$PARAMETER,levels=unique(dfret$PARAMETER))

  return(dfret)
}



