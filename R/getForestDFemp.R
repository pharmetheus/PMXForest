#' getForestDFemp
#'
#' @param covExpressionList A list of expressions that define the covariate values to visualize in the Forest plot.
#' @param metricfunction The function to use to summarise the parameter values in each covariate category. Default is median.
#' @param cGrouping A vector of the same length as covExpressionList to indicate groupings of covariates. Default is no grouping. See example.
#' @param fixedSpacing ??
#' @inheritParams getForestDFSCM
#'
#' @return A data frame
#' @export
#'
#' @examples
#' \dontrun{
#' lsExpr<-rev(list(expression(SEX==2),
#' expression(SEX==1),
#' expression(CRCL>=146),
#' expression(CRCL<94),
#' expression(AGE>=57),
#' expression(AGE<35),
#' expression(WT>104),
#' expression(WT<70),
#' expression(RACEL==3),
#' expression(RACEL==2),
#' expression(RACEL==1),
#' expression(GENO==4),
#' expression(GENO==3),
#' expression(GENO==2),
#' expression(GENO==1),
#' expression(FOOD==1),
#' expression(FOOD==0),
#' expression(FORM==1),
#' expression(FORM==0)))
#'
#' cGrouping = c(1,1,2,2,3,3,3,3,4,4,4,5,5,6,6,7,7,8,8,9,9)
#'
#' covnamesEmp <- c("Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian","African American","Asian and other","WT<70 kg","WT>104 kg",
#'                "Age<35 y","Age>57 y","CRCL<94 mL/min","CRCL>146 mL/min","Male","Female")
#'
#'
#'
#' dfresCOVemp<-getForestDFemp(
#'   metricFunction     = median,
#'   cGrouping          = cGrouping,
#'   dfData             = dfData,
#'   covExpressionsList = lsExpr,
#'   cdfCovsNames       = covnamesEmp,
#'   functionList       = list(paramFunction),
#'   functionListName   = functionListName,
#'   noBaseThetas       = 16,
#'   dfParameters       = dfSamplesCOV,
#'   ncores             = 6)
#' }
getForestDFemp <- function(dfData,
                          covExpressionsList,
                          noBaseThetas,
                          dfParameters,
                          cdfCovsNames = NULL,
                          functionList = list(function(basethetas, dfrow, ...) {
                            return(basethetas[1])
                          }),
                          functionListName = "PAR1",
                          strID = "ID",
                          quiet = TRUE,
                          metricFunction = median,
                          probs = c(0.05, 0.5, 0.95),
                          dfRefRow = NULL,
                          cGrouping = NULL,
                          fixedSpacing = TRUE,
                          groupdist = 0.3,
                          withingroupdist = 0.2,
                          ncores = 1,
                          cstrPackages = NULL,
                          cstrExports = NULL,
                          ...) {

  ## Check input
  if(!is.null(cdfCovsNames) & (length(cdfCovsNames) != length(covExpressionsList))) stop("cdfCovsNames should have the same length as covExpressionsList.")

  ## Remove samples with problems. Will use THETA1 == NA as an indicator for a problematic samples
  dfParameters <- dfParameters %>% filter(!is.na(THETA1))

  ## Create a fake dfCovs
  dfCovs <- data.frame(COV=rep(-99,length(covExpressionsList)))


  resList <- list()
  dfData$TMPINDEX1 <- 1:nrow(dfData) # Add a temp index to row evaluate later on

  if (is.null(cGrouping)) cGrouping <- 1:length(covExpressionsList)

  if (ncores > 1) {
    registerDoParallel(cores = ncores)
  }
  dfres <- foreach(
    k = 1:nrow(dfParameters), .packages = cstrPackages,
    .export = cstrExports, .verbose = !quiet, .combine = bind_rows
  ) %dopar% {
    thetas <- as.numeric(dfParameters[k, 1:noBaseThetas])
    dftmp <- data.frame()

    ## Calculate TVpars
    # For all subjects, calculate the functionList plot
    for (m in 1:nrow(dfData)) {
      for (j in 1:length(functionList)) {
        val <- functionList[[j]](thetas = thetas, df = dfData[m, ], ...)
        listcount <- length(val)
        n <- 1

        for (l in 1:listcount) {
          dftmp <- bind_rows(dftmp, data.frame(
            ITER = k,
            SUBJ = m, NAME = functionListName[n], VALUE = val[[l]],
            stringsAsFactors = FALSE
          ))
          n <- n + 1
        }
      }
    }


    # The default ref=1
    valbase <- rep(1, length(functionListName))

    # Calculate the ref for this sample
    for (j in 1:length(functionList)) {
      if (!is.null(dfRefRow)) {
        valbase_tmp <- functionList[[j]](thetas = thetas, df = dfRefRow, ...)


        listcount <- length(valbase_tmp)
        n <- 1

        for (l in 1:listcount) {
          valbase[n] <- valbase_tmp[[l]]
          n <- n + 1
        }
      }
    }

    if (is.null(dfRefRow)) { # Calculate a reference based on observed data
      for (j in 1:length(functionListName)) {
        valbase[j] <- metricFunction(subset(dftmp, ITER == k & NAME == functionListName[j])$VALUE)
      }
    }


    ## Compute statistic for each "cov value"
    dfrest <- data.frame()
    for (i in 1:length(covExpressionsList)) { # For all rows in the forestPlot
      # Get the subjects included in this Expression
      subjs <- subset(dfData, eval(covExpressionsList[[i]]))$TMPINDEX1
      if (is.null(subjs) || length(subjs) == 0) {
        stop(paste0("Error: no available data for subset: ", as.character(covExpressionsList[[i]])))
      }
      NSubjs <- length(unique(subset(dfData, eval(covExpressionsList[[i]]))[strID]))
      # Get the values out
      dft <- subset(dftmp, ITER == k & SUBJ %in% subjs)
      for (m in 1:length(functionListName)) {
        val <- metricFunction(subset(dft, NAME == functionListName[m])$VALUE)
        # Store the val
        dfrest <- bind_rows(dfrest, data.frame(
          ITER = k,
          SUBJ = -1, NAME = functionListName[m], VALUE = val,
          COVS = i, NSUBJSCOVGROUP = NSubjs,
          VALUEBASE = valbase[m], stringsAsFactors = FALSE
        ))
      }
    }
    dfrest
  }

  dfret <- data.frame()
  for (i in 1:length(covExpressionsList)) {
    if (is.null(cdfCovsNames)) {
      covname <- as.character(covExpressionsList[[i]])
    }
    else {
      covname <- cdfCovsNames[i]
    }
    group <- cGrouping[i]
    for (j in 1:length(functionListName)) {
      dft <- dfres[dfres$COVS == i & dfres$NAME == functionListName[j] &
        dfres$SUBJ == -1, ]
      quant <- quantile(dft$VALUE,
        probs = probs, names = FALSE,
        na.rm = T
      )
      mean_base <- mean(dft$VALUEBASE)
      median_base <- median(dft$VALUEBASE)
      true_base <- dft$VALUEBASE[dft$ITER == 1]
      dfrow <- cbind(dfCovs[i, ], data.frame(
        GROUP = group,
        COVNUM = i, COVNAME = covname, PARAMETER = functionListName[j],
        REFMEAN = mean_base, REFTRUE = true_base, REFMEDIAN = median_base
      ))
      for (k in 1:length(probs)) {
        dfp <- data.frame(X1 = 1)
        dfp[[paste0("q", k)]] <- quant[k]
        dfrow <- cbind(dfrow, dfp[, 2])
        names(dfrow)[ncol(dfrow)] <- paste0("q", k)
      }
      dfret <- rbind(dfret, dfrow)
    }
  }
  if (!fixedSpacing) {
    dfret$Y <- NA
    for (i in 1:length(sort(unique(dfret$GROUP)))) {
      dft <- dfret[dfret$GROUP == sort(unique(dfret$GROUP))[i], ]
      num_in_group <- nrow(dft) / length(unique(dfret$PARAMETER))
      for (n in 1:length(unique(dft$COVNUM))) {
        dfret$Y[dfret$GROUP == sort(unique(dfret$GROUP))[i] &
          dfret$COVNUM == unique(dft$COVNUM)[n]] <- (i -
          1) + n / (num_in_group + 1)
      }
    }
  }
  else {
    dfret$Y <- NA
    a <- 0
    for (i in 1:length(sort(unique(dfret$GROUP)))) {
      dft <- dfret[dfret$GROUP == sort(unique(dfret$GROUP))[i], ]
      num_in_group <- nrow(dft) / length(unique(dfret$PARAMETER))
      for (n in 1:length(unique(dft$COVNUM))) {
        dfret$Y[dfret$GROUP == sort(unique(dfret$GROUP))[i] &
          dfret$COVNUM == unique(dft$COVNUM)[n]] <- a
        a <- a + withingroupdist
      }
      a <- a + groupdist
    }
  }
  stopImplicitCluster()
  return(dfret)
}
