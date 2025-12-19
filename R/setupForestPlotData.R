#' Set up the data.frame for \code{forestPlot} to use
#'
#' @description Creates the data.frame that is actually used for generating the Forest plots. Used by \code{forestPlot}.
#'
#' @inheritParams getPlotVars
#' @param dfres A data.frame as output from getForestDFSCM, getForestDFFREM or getForestDFemp.
#' @param parameters A character vector with the parameters present in the dfres$PARAMETER column to be included in the output.
#' @param parameterLabels A vector of labels for the parameters. Should either have the same length as \code{parameters} or the same length as the number of rows in \code{dfres}.
#' If a vector of the same length as \code{parameters} should also be ordered to match \code{parameters}. If \code{parameters} is not specified, then is should match the order in \code{dfres}. Is by default the same as \code{parameters}. Used for faceting the errorbar panels in the Forest plot in the x-direction (columns).
#' @param parameterLabelsPrefix A string to be added before the \code{parameterLabels} in the facet strip for the parameter panels.
#' @param groupNameLabels A vector of labels for the covariate groups. Should either have the same length as the number unique values in \code{dfres$GROUPNAME} or the same length as the number of rows in \code{dfres}.
#'  Is by default the same as \code{dfres$GROUPNAME}.
#' @param statisticsLabels A character string that will precede the \code{parameterLabels} in the facet labels for the statistics panels. Default is `Statistics:`.
#' @param sigdigits An integer number specifying the number of significant digits to use in the statistics tables.
#' @param onlySignificant Logical. Should only the significant covariates be included (TRUE) or all covariates regardless of significance (FALSE).
#'
#' @return A processed data.frame to the used for creating the Forest plot. Only the columns used in the actual Forest plot is included:
#' \describe{
#' \item{GROUPNAME}{The names of the covariate group, as specified in the call to \code{\link{createInputForestData}}.}
#' \item{GROUPNAMELABEL}{The names of the covariate group as they will appear in the y-axis Facet. Is by default the same as \code{GROUPNAME}.}
#' \item{COVNAME}{The label for each row in each facet in the Forest plot. For exampe Males and Females for the GROUPNAME Sex.}
#' \item{PARAMETER}{The parameters to be plotted in the Forest plot. Comes from \code{functionListName} in the call to \code{\link{getForestDFSCM}}.}
#' \item{PARAMETERLABEL}{The labels for \code{PARAMETER} as they will appear in the Forest plot. Is by default the same as \code{PARAMETER}. Used for faceting the errorbar panels in the Forest plot in the x-direction (columns).}
#' \item{REF}{The parameter value to be used for the reference line, either REFTRUE or REFFUNC in \code{dfres} or 1 if \code{plotRelative=TRUE}.}
#' \item{point}{The value used to plot the point estimate for each \code{COVNAME}.}
#' \item{q1 and q2}{The values used to plot the lower (q1) and (q2) interval edges of the confidence intervals.}
#' \item{reference}{A character string indicating what was used as the reference, either \code{func} or \code{true}, for the parameter function or the first row in \code{dfParameters}}
#' \item{STATISTIC}{A character variable with the numerical summary for each \code{COVNAME} to be used in the statistics table.}
#' \item{STATISTICSLABEL}{Similar to \code{PARAMETERLABEL}. Is the label used for \code{PARAMETER} in the table plots. Used for faceting the table panels in the Forest plot in the x-direction (columns).}
#' \item{onlySignificant}{Logical. Should only GROUPNAME with with TRUE in the COVEFF column be included in the Forest plot?}
#' \item{setSignEff}{NULL of a list. If a list, it should be a list ov vectors with a string for PARAMETER as the forst element and a string for GROUPNAME as the second element. If they match PARAMETER and GROUPNAME in dfres then
#' COVEFF will be set to TRUE else be set to FALSE. This will override the existing values in COVEFF.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' plotData<- setupForestPlotData(dfres)
#' }
setupForestPlotData <- function(dfres,
                                parameters            = unique(dfres$PARAMETER),
                                parameterLabels       = NULL,
                                parameterLabelsPrefix = NULL,
                                groupNameLabels       = NULL,
                                statisticsLabels      = NULL,
                                plotRelative          = TRUE,
                                noVar                 = FALSE,
                                reference             = "func",
                                sigdigits             = 2,
                                onlySignificant       = FALSE,
                                setSignEff            = NULL) {

  ## Input checks
  if(!is.null(parameterLabels)) {
    if(!((length(parameterLabels) == length(parameters)) |
         (length(parameterLabels) ==nrow(dfres)))) {
      stop("The number of parameter labels must either be the same as the number of parameters or have the same length as the number of rows in dfres.")
    }
  }

  if(!is.null(groupNameLabels)) {
    if(!((length(groupNameLabels) == length(unique(dfres$GROUPNAME))) |
         (length(groupNameLabels)  == nrow(dfres)))) {
      stop("The number of group name labels must either be the same as the number of unique values in GROUPNAME or have the same length as the number of rows in dfres.")
    }
  }

  ## Filter the parameters to use and redefine the levels definition (in case it is different from default)
  dfres <- dfres %>%
    filter(PARAMETER %in% parameters) %>%
    mutate(PARAMETER=factor(PARAMETER,levels=parameters)) %>%
    arrange(COVNUM,PARAMETER)


  ## Define the PARAMETERLABEL column
  if(!is.null(parameterLabels)) {
    dfres$PARAMETERLABEL <- parameterLabels
    dfres$PARAMETERLABEL <- factor(dfres$PARAMETERLABEL,levels=unique(dfres$PARAMETERLABEL))
  } else {
    dfres$PARAMETERLABEL <- dfres$PARAMETER
  }

  ## Name the STATISTICS column
  if(is.null(statisticsLabels)) {
    dfres$STATISTICSLABEL <- dfres$PARAMETERLABEL
  } else {
    dfres$STATISTICSLABEL <- factor(paste0(statisticsLabels,dfres$PARAMETERLABEL),levels = paste0(statisticsLabels,unique(dfres$PARAMETERLABEL)))
  }

  ## Add the parameter labels prefix if requested.
  if(!is.null(parameterLabelsPrefix)) {
    dfres$PARAMETERLABEL <- paste0(parameterLabelsPrefix, dfres$PARAMETERLABEL)
    dfres$PARAMETERLABEL <- factor(dfres$PARAMETERLABEL,levels=unique(dfres$PARAMETERLABEL))
  }

  ## Name the GROUPNAME column
  if(!is.null(groupNameLabels)) {

    if(length(groupNameLabels) == length(unique(dfres$GROUPNAME))) {
      names(groupNameLabels) <- unique(dfres$GROUPNAME)

      dfres$GROUPNAMELABEL <- dfres$GROUPNAME

      for (gName in unique(dfres$GROUPNAME)) {
        dfres$GROUPNAMELABEL <- ifelse(dfres$GROUPNAME == gName,groupNameLabels[gName],dfres$GROUPNAMELABEL)
      }

    } else {
      dfres$GROUPNAMELABEL <- groupNameLabels
    }

  } else {
    dfres$GROUPNAMELABEL <- dfres$GROUPNAME
  }

  # Retain order
  dfres$GROUPNAMELABEL <- factor(dfres$GROUPNAMELABEL,levels=unique(dfres$GROUPNAMELABEL[order(dfres$GROUPNAME)]))


  ## Determine which statistic to use
  vars <- getPlotVars(plotRelative,noVar,reference)

  ## Set COVEFF if setSignEff is non-null
  if(!is.null(setSignEff)) {
    dfres <- setCOVEFF(dfres,setSignEff)
  }

  ## Remove the non-significant covariates if onlySignificant
  if(onlySignificant) {
    dfres <- droplevels(dfres %>% filter(PARAMETER%in% parameters) %>% group_by(GROUPNAME) %>% filter(any(COVEFF==TRUE)))
  }

  ## Create the plot data
  plotData <- dfres %>%
    select(GROUPNAME,GROUPNAMELABEL,COVNUM,COVEFF,COVNAME,PARAMETER,PARAMETERLABEL,STATISTICSLABEL,!!vars) %>%
    # select(GROUPNAME,GROUPNAMELABEL,COVNUM,COVEFF,COVNAME,PARAMETER,PARAMETERLABEL,STATISTICSLABEL,REF=vars["ref"],point=vars["point"],q1=vars["q1"],q2=vars["q2"]) %>%
    mutate(reference=reference) %>%
    group_by(PARAMETER,COVNAME) %>%
    mutate(REF=ifelse(plotRelative,1,REF)) %>%
    ungroup %>%
    #mutate(COVNAME = factor(COVNUM,labels=unique(COVNAME))) %>%
    mutate(
      meanlabel  = table1::signif_pad(point, sigdigits),
      lowcilabel = table1::signif_pad(q1, sigdigits),
      upcilabel  = table1::signif_pad(q2,sigdigits),
      STATISTIC  = paste0(meanlabel, " [", lowcilabel, "-", upcilabel,"]"),
      STATISTIC  = stringr::str_pad(STATISTIC,max(stringr::str_length(STATISTIC)),side="right",' ')
    ) %>%
    select(-meanlabel,-lowcilabel,-upcilabel,-COVNUM)

  return(plotData)
}
