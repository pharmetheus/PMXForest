#' getForestDFFREM
#'
#' @description Get a data frame with Forest border for each univariate or multivariate covariate (and value(s)) in the input data frame. If a list a data frame will be created from the list, see function dfCreateInputForestData
#'
#' @import doParallel
#' @import foreach
#' @import dplyr
#' @import stringr
#'
#' @param noSigmas  Number of sigma (epsilon) parameters in the model.
#' @param noParCov Number of parameters for which covariate relations are sought (often the same as noBaseThetas).
#' @param noSkipOm Nof diag omegas (variances) that should not be part of the FREM calculations. Such omegas has to come before the large FREM omega block.
#' @param parNames Names of the parameters
#' @param covNames Names of the covariates
#' @param availCov Names of the covariates to use in the calculation of the FFEM model.
#' @inheritParams getForestDFSCM
#'
#'
#' @return A data frame with summary statistics for each parameters and covariate combinations:
#'
#' @export
#'
#' @examples
#' \dontrun{
#'dfresCOVfrem <-getForestDFFREM(dfCovs           = dfCovs,
#'                               cdfCovsNames     = covnames,
#'                               covNames         = getCovNames(modFile),
#'                               functionList     = list(paramFunction),
#'                               functionListName = functionListName,
#'                               noBaseThetas     = noBaseThetas,
#'                               noParCov         = 4,
#'                               noSigmas         = 2,
#'                               noSkipOm         = 2,
#'                               noCovThetas      = noCovThetas,
#'                               dfParameters     = dfSamplesCOVfrem,
#'                               probs            = c(0.05, 0.95),
#'                               dfRefRow         = NULL,
#'                               quiet            = TRUE,
#'                               ncores           = 1,
#'                               cstrPackages     = c("PMXFrem","dplyr"))
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
                            probs = c(0.05, 0.95),
                            pointFunction = median,
                            dfRefRow = NULL,
                            cGrouping = NULL,
                            ncores = 1,
                            cstrPackages = NULL,
                            cstrExports = NULL,
                            iMiss = -99,
                            ...) {


  if(!is.list(covNames)) stop("covNames must be a list")

  if (!is.null(dfRefRow) && nrow(dfRefRow)!=1 && nrow(dfRefRow)!=nrow(dfCovs)) {
    stop("The number of reference rows (dfRefRow) should be either NULL (missing used as reference), one (this row used as reference) or equal to dfCovs (change reference for each covariate combination)")
  }
  ## Try to make dfCovs into a data.frame if it isn't that already
  if (!is.data.frame(dfCovs)) {
    dfCovs <- createInputForestData(dfCovs)
  }
  ## Replace potential NAs in dfCovs with the missing value token
  dfCovs[is.na(dfCovs)] <- iMiss

  #If a needed covariate is not present in dfCovs, set it to missing
  if(!all(covNames$covNames %in% names(dfCovs))) {
    if (!quiet) warning("Not all covariates in frem model are present in dfCovs, setting them to missing")
    dfCovs[,covNames$covNames[!(covNames$covNames %in% names(dfCovs))]] <-iMiss
  }



  groupnames<-NULL #Store the temp groupnames
  if (any(names(dfCovs)=="COVARIATEGROUPS")) {
    groupnames<-dfCovs[,"COVARIATEGROUPS"]
    dfCovs[,"COVARIATEGROUPS"]<-NULL
  }

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

      if (any(!currentNames %in% covNames$covNames) && !quiet) {
        warning(paste0("Can't find some of the covariates: ", currentNames, " in the FREM model, perhaps they are structural covariates!"))
      }

      ## Calculate the ffemObj in case reference covariates have been specified
      if (!is.null(dfRefRow)) {

          indi <- min(i,nrow(dfRefRow))

          ffemObjRef <- PMXFrem::calcFFEM(
            noBaseThetas = noBaseThetas,
            noCovThetas  = noCovThetas,
            noSigmas     = noSigmas,
            dfext        = dfext,
            covNames     = covNames$covNames,
            availCov     = names(dfRefRow[indi,])[as.numeric(dfRefRow[indi,]) != iMiss][names(dfRefRow[indi,])[as.numeric(dfRefRow[indi,]) != iMiss] %in% covNames$covNames],
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
        availCov     = names(dfCovs[i, ])[as.numeric(dfCovs[i,]) != iMiss][names(dfCovs[i, ])[as.numeric(dfCovs[i,]) != iMiss] %in% covNames$covNames],
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
          data47_jxrtp_ref<-dfRefRow[indi,]
          ffem_expr_base <- str_replace_all(ffemObjRef$Expr[j],pattern = "data\\$", replacement = "data47_jxrtp_ref$")
          coveffects_base[j] <- as.numeric(eval(parse(text = ffem_expr_base)))
        }
      }

      ## Send the expresssions to the functions in functionList
      n <- 1
      for (j in 1:length(functionList)) {

        val <- functionList[[j]](thetas, coveffects,dfrow = dfCovs[i, ], ...)

        ## Do it for the reference
        if (!is.null(dfRefRow)) {
          valbase <- functionList[[j]](basethetas = thetas,covthetas = coveffects_base, dfrow = dfRefRow[indi,], ...)
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
    groupname<-group
    if (!is.null(groupnames)) groupname<-groupnames[i]

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
      dfrow       <- cbind(dfCovs[i, ],
        data.frame(GROUP = group,
        GROUPNAME = groupname,
        COVNUM = i, COVNAME = covname, PARAMETER = functionListName[j],
        REFFUNC = func_base, REFFINAL = true_base, POINT=FUNCVAL, POINT_NOVAR_REL_REFFUNC=FUNCNOVAR,
        POINT_REL_REFFUNC = FUNCVAL/func_base, POINT_REL_REFFINAL = FUNCVAL/true_base,
        COVEFF = !all(dft$RELINTERNAL==1)))

      for (k in 1:length(probs)) {
        dfp <- data.frame(X1 = 1)
        dfp[[paste0("Q", k)]] <- quant[k]
        dfp[[paste0("Q",k,"_REL_REFFUNC")]] <-  quant[k]/func_base
        dfp[[paste0("Q",k,"_REL_REFFINAL")]] <-  quant[k]/true_base
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

  stopImplicitCluster()

  ## Add a column with YES/NO depending on if refRow was provided or if it was set to the default NULL
  dfret <- dfret %>% mutate(REFROW = ifelse(is.null(dfRefRow),"NO","YES"))

  ## Make sure GROUPNAME, COVNAME, PARAMETER are factors
  dfret$GROUPNAME <- factor(dfret$GROUPNAME,levels=unique(dfret$GROUPNAME))
  dfret$COVNAME   <- factor(dfret$COVNAME,levels=unique(dfret$COVNAME))
  dfret$PARAMETER <- factor(dfret$PARAMETER,levels=unique(dfret$PARAMETER))

  return(dfret)
}


