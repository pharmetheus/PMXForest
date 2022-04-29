#' Get the correct dfres variables to plot based on user specification
#'
#' @description Determine which variables from dfres to base the forest plot on
#' @keywords internal
#' @param plotRelative Should the plot be made on the relative scale (TRUE or FALSE)
#' @param noVar Should the uncertainty in the reference be included (FALSE) or not (TRUE)
#' @param reference Which reference should be used, the final parameter estimates (final) or the parameter function (func)
#'
#' @return A vector with the names from dfres to use for POINT, q1, q2 and ref
#' @export
#'
#' @examples
#' # Get the relative variables with uncertainty in the reference:
#' # "POINT_REL_REFFUNC"    "Q1_REL_REFFUNC"    "Q2_REL_REFFUNC"           "REFFUNC"
#'
#' getPlotVars()
#'
#' # Get the relative variables without uncertainty in the reference:
#' # "POINT_NOVAR_REL_REFFUNC"    "Q1_NOVAR_REL_REFFUNC"    "Q2_NOVAR_REL_REFFUNC"                 "REFFUNC"
#'
#' getPlotVars(noVar=TRUE)
#'
getPlotVars <- function(plotRelative=TRUE,noVar=FALSE,reference="func") {

  if(!(reference %in% c("func","final"))) {stop("reference needs to be either func or final.")}

  if(!plotRelative & noVar & reference == "func") {
    vars <- c(point="POINT",q1="Q1",q2="Q2",ref="REFFUNC")

  } else if(!plotRelative & noVar & reference == "final") {
    vars <- c(point="POINT",q1="Q1",q2="Q2",ref="REFFINAL")

  } else if(!plotRelative & !noVar & reference == "func") {
    stop("The combination of plotRelative=FALSE and noVar=FALSE is not possible\n.")

  } else if(!plotRelative & !noVar & reference == "final") {
    stop("The combination of plotRelative=FALSE and noVar=FALSE is not possible\n.")

  } else if(plotRelative & noVar & reference == "func") {
    vars <- c(point="POINT_NOVAR_REL_REFFUNC",q1="Q1_NOVAR_REL_REFFUNC",q2="Q2_NOVAR_REL_REFFUNC",ref="REFFUNC")

  } else if(plotRelative & noVar & reference == "final") {
    stop("The combination of plotRelative=TRUE, noVar=TRUE and reference=final is not possible\n.")

  } else if(plotRelative & !noVar & reference == "func") {
    vars <- c(point="POINT_REL_REFFUNC",q1="Q1_REL_REFFUNC",q2="Q2_REL_REFFUNC",ref="REFFUNC")

  } else if(plotRelative & !noVar & reference == "final") {
    vars <- c(point="POINT_REL_REFFINAL",q1="Q1_REL_REFFINAL",q2="Q2_REL_REFFINAL",ref="REFFINAL")

  }
  return(vars)
}


#' Set up the data.frame for \code{forestPlot} to use
#'
#' @description Creates the data.frame that is actually used for generating the Forest plots. Used by \code{forestPlot}.
#'
#' @inheritParams getPlotVars
#' @param dfres A data.frame as output from getForestDFSCM, getForestDFFREM or getForestDFemp.
#' @param parameters A character vector with the parameters present in the dfres$PARAMETER column to be included in the output.
#' @param parameterLabels A vector of labels for the parameters. Should either have the same length as \code{parameters} or the same length as the number of rows in \code{dfres}.
#' Is by default the same as \code{parameters}. Used for faceting the errorbar panels in the Forest plot in the x-direction (columns).
#' @param groupNameLabels A vector of labels for the covariate groups. Should either have the same length as the number unique values in \code{dfres$GROUPNAME} or the same length as the number of rows in \code{dfres}.
#'  Is by default the same as \code{dfres$GROUPNAME}.
#' @param statisticsLabels A vector of labels for the parameters in the table plots. Should either have the same length as \code{parameters} or the same length as the number of rows in \code{dfres}.
#'  Can be \code{plotmath} expressions. Is by default the same as \code{dfres$PARAMETER}.
#' @param sigdigits An integer number specifying the number of significant digits to use in the statistics tables.
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
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' plotData<- setupForestPlotData(dfres)
#' }
setupForestPlotData <- function(dfres,parameters=unique(dfres$PARAMETER),parameterLabels=NULL,groupNameLabels=NULL,statisticsLabels=NULL,
                                plotRelative=TRUE,noVar=FALSE,reference="func",sigdigits=2) {

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

  if(!is.null(statisticsLabels)) {
    if(!((length(statisticsLabels) == length(parameters)) |
         (length(statisticsLabels) ==nrow(dfres)))) {
      stop("The number of statistics labels must either be the same as the number of parameters or have the same length as the number of rows in dfres.")
    }
  }
 # if(!is.null(statisticsLabels) & length(statisticsLabels) != length(unique(parameters))) stop("The number of statistics labels must be the same as the number of parameters.")

  ## Filter the parameters to use
  dfres <- dfres %>% filter(PARAMETER %in% parameters)

  ## Name the PARAMETER column
  if(!is.null(parameterLabels)) {
    dfres$PARAMETERLABEL <- parameterLabels
  } else {
    dfres$PARAMETERLABEL <- dfres$PARAMETER
  }

  # Retain order
  dfres$PARAMETERLABEL <- factor(dfres$PARAMETERLABEL,levels=unique(dfres$PARAMETERLABEL[order(dfres$PARAMETER)]))

  ## Name the GROUPNAME column
  if(!is.null(groupNameLabels)) {

    if(length(groupNameLabels) == length(unique(dfres$GROUPNAME))) {
      names(groupNameLabels) <- unique(dfres$GROUPNAME)

      dfres <- dfres %>%
        mutate(GROUPNAME2 = GROUPNAME) %>%
        group_by(GROUPNAME2) %>%
        mutate(GROUPNAMELABEL = groupNameLabels[GROUPNAME[1]]) %>%
        ungroup %>%
        select(-GROUPNAME2)
    } else {
      dfres$GROUPNAMELABEL <- groupNameLabels
    }

  } else {
    dfres$GROUPNAMELABEL <- dfres$GROUPNAME
  }

  # Retain order
  dfres$GROUPNAMELABEL <- factor(dfres$GROUPNAMELABEL,levels=unique(dfres$GROUPNAMELABEL[order(dfres$GROUPNAME)]))

  ## Name the STATISTICS column
  if(!is.null(statisticsLabels)) {
    dfres$STATISTICSLABEL <- statisticsLabels
  } else {
    dfres$STATISTICSLABEL <- dfres$PARAMETERLABEL
  }

  # Retain order
  dfres$STATISTICSLABEL <- factor(dfres$STATISTICSLABEL,levels=unique(dfres$STATISTICSLABEL[order(dfres$PARAMETER)]))

  ## Determine which statistic to use
  vars <- getPlotVars(plotRelative,noVar,reference)

  plotData <- dfres %>%
    select(GROUPNAME,GROUPNAMELABEL,COVNUM,COVEFF,COVNAME,PARAMETER,PARAMETERLABEL,STATISTICSLABEL,REF=vars["ref"],point=vars["point"],q1=vars["q1"],q2=vars["q2"]) %>%
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


#'Forest plots
#'
#' @import dplyr ggplot2 ggpubr
#'
#' @description Create Forest plots consisting of alternating columns of errorbar plots and tabulated numerical statistics.
#' @inheritParams setupForestPlotData
#' @param plotData A \code{data.frame} to be used for creating the Forest plot. If \code{plotData} is provided \code{dfres} is ignored. See Details.
#' @param referenceParameters Character string indicating how the reference parameter estimates were obtained.
#' "final" means the final parameter estimates from the model. "func" means the mean parameter estimates across the posterior distribution.
#' @param referenceInfo A character string or NULL. If set to NULL, no information about how the reference line was derived will be displayed. If set to "auto",
#' generic information based on \code{referenceParameters} and the \code{REFROW} column in \code{dfres} will be included at the bottom of the plot. If not NULL or "auto",
#' then the character string will be displayed below the plot. Uses \code{ggpubr::annotate_figure} and \code{ggpubr::text_grob} to place and format the text.
#' Formatting instructions to \code{ggpubr::text_grob} are passed on through \code{...}. Default is "auto".
#' @param labelfun A label function compatible with \code{labeller}. Used to format \code{parameterLabels} used as column facet labels for the errorbar and table plots.
#' Default is \code{label_value}.
#' @param groupname_labelfun A label function compatible with \code{labeller}. Used to format \code{groupNameLabels} used as row facet labels for the errorbar and table plots.
#' Default is \code{label_value}.
#' @param ref_area Numerical vector indicating the horizontal size of the reference area. The default is \code{c(0.8,1.25)}. Used to multiply the reference value to derive the actual
#' xmin and xmax values for the reference area.
#' @param ref_fill_col Reference area color.
#' @param ref_fill_alpha Reference area alpha value.
#' @param ref_line_size Reference line size.
#' @param ref_line_type Reference line linetype.
#' @param ref_line_col Reference line color.
#' @param ci_line_type Errorbar linetype.
#' @param ci_line_col Errorbar color.
#' @param ci_line_size Errorbar size.
#' @param point_shape Point estimate shape.
#' @param point_size Size of the point estimate symbol.
#' @param tabTextSize The size (pt) of the text in the table plots.
#' @param point_color Point estimate color.
#' @param strip_right_size Size of the facet text on the right (in pt). Default is NULL, which will fall back on the theme default.
#' @param strip_top_size Size of the facet text on the top (in pt). Default is NULL, which will fall back on the theme default.
#' @param ref_subj_label A character string indicating the label to be used for the reference subject in the figure legend. Default is "Reference subject".
#' @param ref_area_label A character string indicating the label to be used for the reference area in the figure legend. Default is "Reference area".
#' @param point_label A character string indicating the label to be used for the point estimates in the figure legend. Default is "Point estimate".
#' @param ci_label A character string indicating the label to be used for the confidence interval in the figure legend. Default is "Confidence interval".
#' @param statisticsLabel A vector of character string of the same length as \code{parameters} to be used as facet labels for the table plots.
#' @param xlb x-axis label for the errorbar plots
#' @param return Either "plot" (default), "plotList" or "data". "plot" returns the \code{ggpubr::ggarrange} object. "plotList" returns a list of the different plots that make up
#' the Forest plot (see Details). "data" returns a \code{data.frame} with the actual data used to create the Forest plot.
#' @param table Logical. Should the table plots be included in the Forest plot or not.
#' @param rightStrip Should the facet title in the rightmost panel in the Forest plot be displayed? If TRUE (the default) \code{groupNameLabels} will be used as the facets labels.
#' @param keepYlabs  Logical. Should the labels on the y-axis be kept for all errorbar panels.
#' @param keepRightStrip Logical. Should the right facet titles be kept for all table plots. Only when \code{rightStrip} is \code{TRUE}.
#' @param stackedPlots Should the plots for the parameters be stacked instead of being vertical. Useful if there are many parameters to visualize. Works best with \code{keepYlab} and \code{keepRightStrip} set to \code{TRUE}.
#' @param errbartabwidth A numerical vector of length 2 times the number of parameters included in the Forest plot, specifying the relative width of the errorbar and table plots. Have no effect when \code{table} is FALSE.
#' @param errbarplotscale Scaling factor for the width of the leftmost errorbar plot to compensate for y-axis labels.
#' @param tabplotscale Scaling factor for the width of the rightmost column (usually a table plot) to adjust for the size of the right strip.
#' @param ... Arguments passed on to \code{ggpubr::text_grob}.
#'
#' @details
#' If \code{plotData} is NULL, \code{dfres} is passed to \code{setupForestPlotData} to create a \code{data.frame} that contains the exact data to be plotted, including the numerical statistics.
#' If \code{return} is set to "data" then the \code{data.frame} is returned. This makes it possible to to review the data and/or do modifications if necessary. \code{forestPlot} can
#' be called with \code{plotData} set to the name of the \code{data.frame}, in which case \code{dfres} is not required and the \code{plotData} \code{data.frame} will be used for
#' creating the Forest plot.
#'
#' The Forest plots are created as a combination of separate errorbar plots and plots with the table information. Each panel is a separate plot, which are combined using \code{ggpubr::ggarrange}.
#' In other words, even if the plot looks like a regular faceted plot it is not. For example, in a one parameter Forest plot with \code{table=TRUE}, the left panel is an errorbar plot with the right
#' facet labels suppressed and the right plot is a plot with text, and which have the y-axis labels suppressed. In a plot with three parameters the errobar plot for the middle parameter will have both the y-axis labels and
#' facet labels suppressed. the same is true for that parameter's table plot. The arguments \code{keepYlabs} and \code{keepRightStrip} controls if the y-axis labels and right facet labels for panels "in the middle" should
#' keep the axis and facet labels or not.
#'
#' \code{stackedPlots} switch from the default stacked horizontal orientation used if multiple parameters are included to vertical stacking. For this to provide a nice display it
#' is necessary to have \code{keepYlabs} set to \code{TRUE} and \code{keepRightStrip=TRUE}. \code{stackedPlots} can be combind with \code{table=FALSE} and  \code{rightStrip=FALSE}.
#'
#' The graphical settings use standard \code{ggplot} syntax.
#' @return A \code{ggpubr::ggarrange} object, a list of plots or a \code{data.frame}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a Forest plot, including the table part, for all parameters in dfres.
#'
#' forestPlot(dfres)
#'
#' # Exclude the table part of the plot
#'
#' forestPlot(dfres,table=FALSE)
#'
#' # Create a Forest plot for the specified parameter (which need to be included in dfres).
#'
#' forestPlot(dfres,parameters="CL")
#'
#' # Specify group name labels, the size of the table text and x-axis label
#'
#' forestPlot(dfres,parameters=c("CL"),groupNameLabels = c("Age (y)","Sex","Weight (kg)"),tabTextSize = 20,xlb="Relative parameter value")
#'
#' # Stack the plots instead of plot them horizontally
#'
#' forestPlot(dfres,parameters = c("CL","Frel"),stackedPlots = TRUE,keepYlabs = TRUE,keepRightStrip = TRUE)
#' }
forestPlot <- function(dfres,
                       plotData=NULL,
                       plotRelative=TRUE,
                       noVar = TRUE,
                       referenceParameters = "func",
                       sigdigits = 2,
                       parameters=unique(dfres$PARAMETER), #Labels
                       parameterLabels=parameters,
                       groupNameLabels = NULL,
                       referenceInfo = "auto",
                       labelfun=label_value,
                       groupname_labelfun=label_value,
                       ref_area=c(0.8,1.2),
                       ref_fill_col="gray",
                       ref_fill_alpha=0.5,
                       ref_line_size=1,
                       ref_line_type="dotted",
                       ref_line_col="black",
                       ci_line_type="solid",
                       ci_line_col="blue",
                       ci_line_size=1,
                       point_shape=16,
                       point_color="blue",
                       point_size=3,
                       tabTextSize=10,
                       keepYlabs = FALSE,
                       keepRightStrip = FALSE,
                       stackedPlots   =FALSE,
                       strip_right_size = NULL,
                       strip_top_size = NULL,
                       ref_subj_label = "Reference subject",
                       ref_area_label = "Reference area",
                       point_label    = "Point estimate",
                       ci_label       = "Confidence interval",
                       statisticsLabel = paste("Statistics:",parameterLabels),
                       xlb = ifelse(plotRelative,"Relative parameter value","Parameter value"),
                       return = "plot",
                       table=TRUE,
                       rightStrip=TRUE,
                       errbartabwidth = rep(c(2,1),length(parameters)),
                       errbarplotscale = 1.45,
                       tabplotscale    = 1.1,
                       onlySignificant = FALSE,
                       ...) {


  ## Check if only significant covariates is to be used
  if(onlySignificant) {
     dfres <- droplevels(dfres %>% filter(PARAMETER%in% parameters) %>% group_by(GROUPNAME) %>% filter(any(COVEFF==TRUE)))
  }

  ## Setup the plotting data frame
  if(is.null(plotData)) {
    plotData<- setupForestPlotData(dfres,reference=referenceParameters,plotRelative = plotRelative,parameters=parameters,parameterLabels=parameterLabels,groupNameLabels=groupNameLabels,
                                   statisticsLabel=statisticsLabel,
                                   noVar = noVar,sigdigits = sigdigits)
  } else {
    message("Will use the provided plotData and ignore dfres\n")
  }

  #######
  ## Create two functions, one for the error bar charts and one for the tables
  #######
  parPlot <- function(data,parameters,parameterLabels,label_fun=label_parsed,group_name_label_fun=label_value) {

    ## This function is designed to only deal with one parameter at a time.
    if(length(unique(data$PARAMETER)) !=1 ) stop("Can only deal with one parameter at a time.")


    ref_value <- data %>% distinct(GROUPNAMELABEL,PARAMETERLABEL,REF)

    rect_data <- data.frame(ref_value=ref_value$REF,
                            xmin=ref_value$REF * min(ref_area),
                            xmax = ref_value$REF * max(ref_area),
                            ymin = -Inf, ymax = Inf,
                            GROUPNAMELABEL=ref_value$GROUPNAMELABEL,
                            PARAMETERLABEL=ref_value$PARAMETERLABEL)

    p1a <- ggplot(data,aes(x=point,y=COVNAME,xmin=q1,xmax=q2)) +
      geom_blank() +

      geom_rect(data=rect_data,
                aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="Reference area"),alpha=ref_fill_alpha,
                inherit.aes = FALSE) +

      geom_vline(data=rect_data,
                 aes(xintercept = ref_value,
                     linetype = "Reference subject",
                     color    = "Reference subject"
                 ),
                 size     = ref_line_size,
                 key_glyph = "path") +

      geom_errorbarh(aes(color="CI",linetype="CI"),key_glyph = "path",height=0,size=ci_line_size) +
      geom_point(aes(shape="Point estimate"),color=point_color,size=point_size) +

      scale_fill_manual(name     = NULL, values = c("Reference area"    = ref_fill_col),
                        labels=c(ref_area_label)) +
      scale_color_manual(name    = NULL, values = c("Reference subject" = ref_line_col,"CI"=ci_line_col),
                         labels=c(ci_label,ref_subj_label)) +
      scale_linetype_manual(name = NULL, values = c("Reference subject" = ref_line_type,"CI"=ci_line_type),
                            labels=c(ci_label,ref_subj_label)) +
      scale_shape_manual(name    = NULL, values = c("Point estimate"    = point_shape),
                         labels=c(point_label)) +

      guides(linetype=guide_legend(override.aes=list(size=1))) +
      facet_grid(GROUPNAMELABEL~PARAMETERLABEL,scales = "free",
                 labeller = labeller(PARAMETERLABEL= label_fun,
                                     GROUPNAMELABEL = group_name_label_fun)) +
      ylab(NULL) +
      xlab(xlb) +
      theme(plot.margin = unit(c(5.5,0,5.5,5.5), "pt"))
  }


  tablePlot <- function(data,parameters,tabTextSize=10,label_fun=label_parsed,group_name_label_fun=label_value) {

    p2a <- ggplot(data,aes(x=1,y=COVNAME)) +
      geom_text(aes(label =STATISTIC),size=tabTextSize*0.36) +  # 0.36 will scale size down to regular ggplot size
      facet_grid(GROUPNAMELABEL~STATISTICSLABEL,scales = "free",
                 labeller = labeller(STATISTICSLABEL= label_fun,
                                     GROUPNAMELABEL = group_name_label_fun),space = "free") +
      theme(axis.text = element_blank()) +
      theme(axis.ticks = element_blank()) +
      theme(axis.title = element_blank()) +
      theme(panel.grid.major = element_blank()) +
      theme(panel.grid.minor = element_blank())

    return(p2a)
  }

  ######
  ## Create the plots and determine how each panel in the plot should be formatted in terms of axis text and strip
  ######
  plotList <- list()
  tabList <- list()

  ## The errorbar plots
  for(i in 1:length(parameters)) {

    plotList[[i]] <- parPlot(subset(plotData,PARAMETER==parameters[i]),parameters,parameterLabels,
                             label_fun = labelfun,
                             group_name_label_fun = groupname_labelfun) +
      theme(plot.margin = unit(c(5.5,0,5.5,5.5), "pt")) +
      theme(strip.text.x=element_text(size=strip_top_size))

    ## Deal with y-axis elements
    if(i==1 && i != length(parameters)) {
      plotList[[i]] <- plotList[[i]]
    }

    if(i <= length(parameters) && i!=1) {
      if(!keepYlabs) {
        plotList[[i]] <- plotList[[i]] +
          theme(axis.text.y = element_blank()) +
          theme(axis.ticks.y = element_blank()) +
          theme(axis.title.y = element_blank())
      }
    }

    ## Deal with the strip
    if(i != length(parameters)) {
      plotList[[i]] <- plotList[[i]] +
        theme(strip.text.y=element_blank())
    }

    if(i == length(parameters)) {

      if(table) {
        plotList[[i]] <- plotList[[i]] +
          theme(strip.text.y=element_blank())
      }

      if(!table && !rightStrip) {
        plotList[[i]] <- plotList[[i]] +
          theme(strip.text.y=element_blank())
      }

      if(!table && rightStrip) {
        plotList[[i]] <- plotList[[i]] +
          theme(strip.text.y=element_text(size=strip_right_size))
      }
    }

  }

  ## The table plots
  for(i in 1:length(parameters)) {

    tabList[[i]] <- tablePlot(subset(plotData,PARAMETER==parameters[i]),parameters,
                              label_fun = labelfun,
                              group_name_label_fun = groupname_labelfun,
                              tabTextSize = tabTextSize) +
      theme(plot.margin = unit(c(5.5,5.5,5.5,0), "pt")) +
      theme(strip.text.x=element_text(size=strip_top_size))

    if(i < length(parameters)) {
      if(!keepRightStrip || !rightStrip) {
      tabList[[i]] <- tabList[[i]] +
        theme(strip.text.y=element_blank())
      } else {
        tabList[[i]] <- tabList[[i]] +
          theme(strip.text.y=element_text(size=strip_right_size))
      }
    }

    if(i == length(parameters) & rightStrip) {  # Last panel with right strip
      tabList[[i]] <- tabList[[i]] +
        theme(strip.text.y=element_text(size=strip_right_size))
    }

    if(i == length(parameters) & !rightStrip) {  # Last panel with noright strip
      tabList[[i]] <- tabList[[i]] +
        theme(strip.text.y=element_blank())
    }

  }

  ## If the tables are to be included, assemble a combined list with alternating plots and tables,
  ## otherwise just return  plotList
  totList <- list()
  is.even <- function(x) x %% 2 == 0

  if(table) {
    for(j in 1:(2*length(plotList))) {

      if(!is.even(j)) {
        totList[[j]] <- plotList[[ceiling(j/2)]]
      } else {
        totList[[j]] <- tabList[[j/2]]
      }
    }

  } else {
    totList <- plotList
  }

  if(table) {

    ## Figure out the widths of the plot components
    errbartabwidth[1] <- errbartabwidth[1]*errbarplotscale
    errbartabwidth[length(errbartabwidth)] <- errbartabwidth[length(errbartabwidth)]*tabplotscale

    if(!stackedPlots) {
      myPlot <- ggpubr::ggarrange(plotlist=totList,nrow=1,widths = errbartabwidth,align="h",common.legend = T)
    } else {
      myPlot <- ggpubr::ggarrange(plotlist=totList,nrow=length(totList)/2,ncol=2,widths = errbartabwidth,align="h",common.legend = T)
    }
  } else{
    errbartabwidth[1] <- errbartabwidth[1]*errbarplotscale

    if(!stackedPlots) {
      myPlot <- ggpubr::ggarrange(plotlist=totList,nrow=1,widths = errbartabwidth,align="h",common.legend = T)
    } else {
      myPlot <- ggpubr::ggarrange(plotlist=totList,nrow=length(totList),ncol=1,widths = errbartabwidth,align="h",common.legend = T)
    }
  }

  ## Add information about the reference subject if needed
  if(is.null(referenceInfo)) {
    myPlot <- myPlot
  } else if(!is.null(referenceInfo) & referenceInfo != "auto") {
    myPlot <- ggpubr::annotate_figure(myPlot,bottom=ggpubr::text_grob(referenceInfo,...))
  } else if(referenceInfo=="auto") {

    if(dfres[1,"REFROW"]=="NO" && referenceParameters=="final") {
      refText <- "The reference line is based on the final parameter estimates and the reference covariate values in the model."
    } else if(dfres[1,"REFROW"]=="NO" && referenceParameters=="func") {
      refText <- "The reference line is based on the average parameter estimates over the posterior parameter distribution and the reference covariate values in the model."
    } else if(dfres[1,"REFROW"]=="YES" && referenceParameters=="final") {
      refText <- "The reference line is based on the final parameter estimates and selected covariate values."
    } else if(dfres[1,"REFROW"]=="YES" && referenceParameters=="func") {
      refText <- "The reference line is based on the average parameter estimates over the posterior parameter distribution and selected covariate values."
    }

    myPlot <- ggpubr::annotate_figure(myPlot,bottom=ggpubr::text_grob(refText,...))
  }

  ## return either the arranged plot or the list of plot objects
  if(return=="data") {
    return(plotData)
  } else if(return=="plotList") {
    return(totList)
  } else {
    return(myPlot)
  }
}
