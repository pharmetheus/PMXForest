#' getForestDFemp
#'
#'
#' @import doParallel
#' @import foreach
#' @import dplyr
#'
#' @param covExpressionList A list of expressions that define the covariate values to visualize in the Forest plot.
#' @param metricfunction The function to use to summarise the parameter values in each covariate category. Default is median.
#' @param cGrouping A vector of the same length as covExpressionList to indicate groupings of covariates. Default is no grouping. See example.
#' @inheritParams getForestDFSCM
#'
#' @return A data frame
#' @export
#'
#' @examples
#' \dontrun{
#' lsExpr<-rev(list("SEX"=expression(SEX==2),
#' "SEX"=expression(SEX==1),
#' "CRCL"=expression(CRCL>=146),
#' "CRCL"=expression(CRCL<94),
#' "AGE"=expression(AGE>=57),
#' "AGE"=expression(AGE<35),
#' "WT"=expression(WT>104),
#' "WT"=expression(WT<70),
#' "RACEL"=expression(RACEL==3),
#' "RACEL"expression(RACEL==2),
#' "RACEL"expression(RACEL==1),
#' "GENO"=expression(GENO==4),
#' "GENO"=expression(GENO==3),
#' "GENO"=expression(GENO==2),
#' "GENO"=expression(GENO==1),
#' "FOOD"=expression(FOOD==1),
#' "FOOD"=expression(FOOD==0),
#' "FORM"=expression(FORM==1),
#' "FORM"=expression(FORM==0)))
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
                          probs = c(0.025, 0.975),
                          pointFunction = median,
                          dfRefRow = NULL,
                          cGrouping = NULL,
                          ncores = 1,
                          cstrPackages = NULL,
                          cstrExports = NULL,
                          iMiss=-99,
                          ...) {

  ## Check input
  if(!is.null(cdfCovsNames) & (length(cdfCovsNames) != length(covExpressionsList))) stop("cdfCovsNames should have the same length as covExpressionsList.")

  ## Remove samples with problems. Will use THETA1 == NA as an indicator for a problematic samples
  dfParameters <- dfParameters %>% filter(!is.na(THETA1))

  ## Create a fake dfCovs
  dfCovs <- data.frame(COV=rep(iMiss,length(covExpressionsList)))

  if (!is.null(dfRefRow) && ((is.data.frame(dfRefRow) && nrow(dfRefRow)!=1 && nrow(dfRefRow)!=nrow(dfCovs)) ||
                             (is.expression(dfRefRow[[1]]) && length(dfRefRow)!=1 && length(dfRefRow)!=nrow(dfCovs)))) {
    stop("The number of reference rows/expressions (dfRefRow) should be either NULL (missing used as reference), one (this row used as reference) or equal to covExpressionsList (change reference for each covariate expression)")
  }

  resList <- list()
  dfData$TMPINDEX1 <- 1:nrow(dfData) # Add a temp index to row evaluate later on

  if (is.null(cGrouping)) cGrouping <- 1:length(covExpressionsList)

  groupnames<-names(covExpressionsList)

  ## Register to allow for parallell computing
  registerDoParallel(cores = ncores)


  dfres <- foreach(
    k = 1:nrow(dfParameters), .packages = cstrPackages,
    .export = cstrExports, .verbose = !quiet, .combine = bind_rows
  ) %dopar% {
    thetas <- as.numeric(dfParameters[k, 1:noBaseThetas])
    dftmp <- data.frame()

    ## Calculate TVpars
    # For all subjects, calculate the functionList plot
    for (m in 1:nrow(dfData)) {
      n <- 1
      for (j in 1:length(functionList)) {
        val <- functionList[[j]](thetas = thetas, df = dfData[m, ], ...)
        listcount <- length(val)

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

    dfvalbase<-data.frame() #Prepare for a different reference per expressionList index
    # The default ref=1
    valbase <- rep(1, length(functionListName))
    for (m in 1:length(covExpressionsList)) { #Allow different ref rows per covExpressionList
      n <- 1
      # Calculate the ref for this sample
      for (j in 1:length(functionList)) {
        if (!is.null(dfRefRow) && is.data.frame(dfRefRow)) {
          if (m==1 || nrow(dfRefRow)>1) {
            valbase_tmp <- functionList[[j]](thetas = thetas, df = dfRefRow[min(m,nrow(dfRefRow)),], ...)
            listcount <- length(valbase_tmp)

            for (l in 1:listcount) {
              valbase[n] <- valbase_tmp[[l]]
              n <- n + 1
            }
          }
        }
      }
      dfvalbase<-rbind(dfvalbase,data.frame(rbind(valbase)))
    }


    if (is.null(dfRefRow)) { # Calculate a reference based on observed data. Set the reference to the `metricFunction` of the typical individual predictions
      for (j in 1:length(functionListName)) {
        valbase[j] <- metricFunction(subset(dftmp, ITER == k & NAME == functionListName[j])$VALUE)
      }
      dfvalbase<-as.data.frame(data.frame(rbind(valbase))[rep(1,length(covExpressionsList)),])
    }

    ## Compute statistic for each "cov value"
    dfrest <- data.frame()
    for (i in 1:length(covExpressionsList)) { # For all rows in the forestPlot
      # Get the subjects included in this Expression
      subjs <- subset(dfData, eval(covExpressionsList[[i]]))$TMPINDEX1
      if (is.null(subjs) || length(subjs) == 0) {
        stop(paste0("Error: no available data for subset: ", as.character(covExpressionsList[[i]])))
      }
      NSubjs <- length(unique(subset(dfData, eval(covExpressionsList[[i]]))[[strID]]))
      # Get the values out
      dft <- subset(dftmp, ITER == k & SUBJ %in% subjs)

      #If we have an expression reference
      if (!is.null(dfRefRow) && is.expression(dfRefRow[[1]])) {
        if (i==1 || length(dfRefRow)>1)
          subjsref <- subset(dfData, eval(dfRefRow[[i]]))$TMPINDEX1
        if (is.null(subjsref) || length(subjsref) == 0) {
          stop(paste0("Error: no available data for reference subset: ", as.character(dfRefRow[[i]])))
        }
        # Get the reference values out
        dfr <- subset(dftmp, ITER == k & SUBJ %in% subjsref)
      }


      for (m in 1:length(functionListName)) {
        #Get the metric value
        val <- metricFunction(subset(dft, NAME == functionListName[m])$VALUE)
        #Get the metric ref value if expression ref
        if (!is.null(dfRefRow) && is.expression(dfRefRow[[1]])) {
          VB<-metricFunction(subset(dfr, NAME == functionListName[m])$VALUE)
        } else {
          VB=dfvalbase[i,m]
        }

        # Store the val
        dfrest <- bind_rows(dfrest, data.frame(
          ITER = k,
          SUBJ = -1, NAME = functionListName[m], VALUE = val,
          COVS = i, NSUBJSCOVGROUP = NSubjs,
          VALUEBASE = VB, stringsAsFactors = FALSE
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
    groupname<-groupnames[i]

    for (j in 1:length(functionListName)) {
      dft <- dfres[dfres$COVS == i & dfres$NAME == functionListName[j] &
        dfres$SUBJ == -1, ]
         quant <- quantile(dft$VALUE,
        probs = probs, names = FALSE,
        na.rm = T
      )
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
      dfrow <- cbind(data.frame(COV=as.character(covExpressionsList[[i]])), data.frame(
        GROUP = group,
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

  return(dfret)
}
