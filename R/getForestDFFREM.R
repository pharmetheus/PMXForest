#' getForestDFSCM
#'
#' @description Get a data frame with Forest border for each univariate or multivariate covariate (and value(s)) in the input data frame. If a list a data frame will be created from the list, see function dfCreateInputForestData
#'
#' @param noSigmas  Number of sigma (epsilon) parameters in the model.
#' @param noParCov Number of parameters for which covariate relations are sought (often the same as noBaseThetas).
#' @param noSkipOm Nof diag omegas (variances) that should not be part of the FREM calculations. Such omegas has to come before the large FREM omega block.
#' @param parNames Names of the parameters
#' @param covNames Names of the covariates
#' @param availCov Names of the covariates to use in teh calculation of the FFEM model.
#' @inheritParams getForestDFSCM
#'
#'
#' @return A data frame with summary statistics for each parameters and covariate combinations:
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dfresCOV <- getForestDFSCM(
#'             dfCovs = dfCovs1,
#'             cdfCovsNames = covnames,
#'             functionList = list(paramFunction),
#'             functionListName = functionListName,
#'             noBaseThetas = 16,
#'             dfParameters = dfSamplesCOV,
#'             probs = c(0.05, 0.95),
#'             dfRefRow = dfRefRow,
#'             quiet = TRUE)
#' }

getForestDFFREM <- function(dfCovs,
                            cdfCovsNames = NULL,
                            functionList = list(function(basethetas,covthetas, dfrow, ...) {
                              return(basethetas[1] * exp(covthetas[1]))
                            }),
                            covNames,
                            functionListName = "PAR1",
                            noBaseThetas,
                            dfParameters,
                            noCovThetas,
                            noSigmas,
                            noParCov = noBaseThetas,
                            noSkipOm = 0,
                            parNames = paste("Par",1:noParCov,sep = ""),
                            availCov = covNames,
                            quiet = FALSE,
                            probs = c(0.025, 0.975),
                            pointFunction = median,
                            dfRefRow = NULL,
                            cGrouping = NULL,
                            ncores = 1,
                            cstrPackages = NULL,
                            cstrExports = NULL,
                            iMiss = -99,
                            ...) {


  if(!is.list(covNames)) stop("covNames must be a list")

  if(!all(covNames$covNames %in% names(dfCovs))) stop("All covariates in the frem model need to be present in dfCovs.")

  #resList <- list()


  ## Try to make dfCovs into a data.frame if it isn't that already
  if (!is.data.frame(dfCovs)) {
    dfCovs <- dfCreateInputForestData(dfCovs)
  }
  ## Replace potential NAs in dfCovs with the missing value token
  dfCovs[is.na(dfCovs)] <- iMiss

  ## Create function to sort out grouping
  getGroups <- function(df) {
    cGroups <- c()
    cUnique <- c()
    iGroup <- 0
    for (i in 1:nrow(df)) {
      tmp <- paste0(names(dfCovs[i, ])[as.numeric(dfCovs[i, ]) != iMiss], collapse = ",")
      if (tmp %in% cUnique) {
        tmpl <- which(tmp == cUnique)
        cGroups <- c(cGroups, tmpl)
      }else {
        iGroup <- iGroup + 1
        cGroups <- c(cGroups, iGroup)
        cUnique <- c(cUnique, tmp)
      }
    }
    return(cGroups)
  }

  if (is.null(cGrouping)) {
    cGrouping <- getGroups(dfCovs)
  }

  ## Register to allow for parallell computing
  registerDoParallel(cores = ncores)

  ## Calculate the parameters
  dfres <- foreach( k = 1:nrow(dfParameters), .packages = cstrPackages,
                    .export = cstrExports, .verbose = !quiet, .combine = bind_rows
  ) %dopar% {

    dfext  <- cbind(first = 0, dfParameters[k, ])
    thetas <- as.numeric(dfext[2:(noBaseThetas + 1)])
    dfrest <- data.frame()

    for (i in 1:nrow(dfCovs)) {
      currentNames <- names(dfCovs[i, ])[as.numeric(dfCovs[i, ]) != iMiss]

      if (any(!currentNames %in% covNames$covNames) && quiet) {
        warning(paste0("Can't find some of the covariates: ", currentNames, " in the FREM model, perhaps they are structural covariates!"))
      }

      ## Calculate the ffemObj in case reference covariates have been specified
      if (!is.null(dfRefRow)) {
        ffemObjRef <- PMXFrem::calcFFEM(
          noBaseThetas = noBaseThetas,
          noCovThetas  = noCovThetas,
          noSigmas     = noSigmas,
          dfext        = dfext,
          covNames     = covNames$covNames,
          availCov     = names(dfRefRow)[as.numeric(dfRefRow) != iMiss & names(dfRefRow)[as.numeric(dfRefRow) != iMiss] %in% covNames$covNames],
          quiet        = quiet,
          noSkipOm     = noSkipOm,
          noParCov     = noParCov
        )
      }

      ## Calculate the ffemObj for each set of parameters
      ffemObj <- PMXFrem::calcFFEM(
        noBaseThetas = noBaseThetas,
        noCovThetas  = noCovThetas,
        noSigmas     = noSigmas,
        dfext        = dfext,
        covNames     = covNames$covNames,
        availCov     = names(dfCovs[i, ])[as.numeric(dfCovs[i, ]) != iMiss & names(dfCovs[i, ])[as.numeric(dfCovs[i, ]) != iMiss] %in% covNames$covNames],
        quiet        = quiet,
        noSkipOm     = noSkipOm,
        noParCov     = noParCov
      )

      coveffects      <- rep(0, length(parNames))
      coveffects_base <- rep(0, length(parNames))
      data47_jxrtp    <- dfCovs[i, ]

      ## Compute the ffem expressions
      for (j in 1:length(parNames)) {
        ffem_expr <- str_replace_all(ffemObj$Expr[j], pattern = "data\\$", replacement = "data47_jxrtp$")

        if (length(names(dfCovs[i, ])[as.numeric(dfCovs[i, ]) != iMiss]) != 0) {
          coveffects[j] <- as.numeric(eval(parse(text = ffem_expr)))
        }

        if (!is.null(dfRefRow)) {
          ffem_expr_base <- str_replace_all(ffemObjRef$Expr[j],pattern = "data\\$", replacement = "dfRefRow$")
          coveffects_base[j] <- as.numeric(eval(parse(text = ffem_expr_base)))
        }
      }

      ## Send the expresssions to the functions in functionList
      n <- 1
      for (j in 1:length(functionList)) {
        val <- functionList[[j]](thetas, coveffects,dfrow = dfCovs[i, ], ...)

        ## Do it for the reference
        if (!is.null(dfRefRow)) {
          valbase <- functionList[[j]](basethetas = thetas,covthetas = coveffects_base, dfrow = dfRefRow, ...)
        } else {
          dfmissing <- dfCovs[1, ]
          dfmissing[, ] <- iMiss
          valbase <- functionList[[j]](basethetas = thetas,covthetas = rep(0, length(parNames)), dfrow = dfmissing, ...)
        }

        ## Do it for the parameters
        listcount <- length(val)
        for (l in 1:listcount) {
          dfrest <- bind_rows(dfrest, data.frame(
            ITER = k,
            COVS = i, NAME = functionListName[n], VALUE = val[[l]],
            VALUEBASE = valbase[[l]], stringsAsFactors = FALSE
          ))
          n <- n + 1
        }
      }
    }
    dfrest
  }

  ## Assemble the return data.frame
  getCovNameString <- function(dfrow) {
    strName <- ""
    colnames <- names(dfrow)

    for (i in 1:ncol(dfrow)) {
      if (dfrow[1, i] != iMiss) {
        if (strName == "") {
          strName <- paste0(colnames[i], "=", dfrow[1,i])
        } else {
          strName <- paste0(strName, ", ", colnames[i],"=", dfrow[1, i])
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
      covname <- getCovNameString(dfCovs[i, ])
    } else {
      covname <- cdfCovsNames[i]
    }

    group <- cGrouping[i]

    for (j in 1:length(functionListName)) {
      dft         <- dfres[dfres$COVS == i & dfres$NAME == functionListName[j], ]
      quant       <- quantile(dft$VALUE,probs = probs, names = FALSE, na.rm = T)
      #Calculate the point value of the forest plot
      FUNCVAL=pointFunction(dft$VALUE)

      #Define the relative without parameter uncertainty
      dft$RELINTERNAL <- dft$VALUE/dft$VALUEBASE
      #Get quantile and pointvalue when uncertainty is not taken into account
      quantrel <- quantile(dft$RELINTERNAL,probs = probs, names = FALSE,na.rm = T)
      FUNCNOVAR=pointFunction(dft$RELINTERNAL)

      #Calculate reference value based one the pointFunction
      func_base<-pointFunction(dft$VALUEBASE)
      true_base   <- dft$VALUEBASE[dft$ITER == 1]
      dfrow       <- cbind(dfCovs[i, ], data.frame(GROUP = group,
        COVNUM = i, COVNAME = covname, PARAMETER = functionListName[j],
        REFFUNC = func_base, REFTRUE = true_base, FUNC=FUNCVAL, FUNC_NOVAR=FUNCNOVAR,
        COVEFF = !all(dft$RELINTERNAL==1)))

      for (k in 1:length(probs)) {
        dfp <- data.frame(X1 = 1)
        dfp[[paste0("q", k)]] <- quant[k]
        dfrow <- cbind(dfrow, dfp[, 2])
        names(dfrow)[ncol(dfrow)] <- paste0("q", k)
      }
      for (k in 1:length(probs)) {
        dfp <- data.frame(X1 = 1)
        dfp[[paste0("q", k,"_NOVAR")]] <- quantrel[k]
        dfrow <- cbind(dfrow, dfp[, 2])
        names(dfrow)[ncol(dfrow)] <- paste0("q", k,"_NOVAR")
      }
      dfret <- rbind(dfret, dfrow)
    }
  }

  stopImplicitCluster()
  return(dfret)
}


