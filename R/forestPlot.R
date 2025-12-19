#'Forest plots
#'
#' @import dplyr ggplot2 ggpubr table1
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
#' @param commonXlab Logical. Should a common x-axis title be used or should each panel have its own. Default is FALSE.
#' @param size.legend.text The font size used in the legend in points or relative size to the theme base font. Default is rel(0.8).
#' @param return Either "plot" (default), "plotList" or "data". "plot" returns the \code{ggpubr::ggarrange} object. "plotList" returns a list of the different plots that make up
#' the Forest plot (see Details). "data" returns a \code{data.frame} with the actual data used to create the Forest plot.
#' @param table Logical. Should the table plots be included in the Forest plot or not.
#' @param rightStrip Should the facet title in the rightmost panel in the Forest plot be displayed? If TRUE (the default) \code{groupNameLabels} will be used as the facets labels.
#' @param keepYlabs  Logical. Should the labels on the y-axis be kept for all errorbar panels.
#' @param keepRightStrip Logical. Should the right facet titles be kept for all table plots. Only when \code{rightStrip} is \code{TRUE}.
#' @param stackedPlots Should the plots for the parameters be stacked instead of being vertical. Useful if there are many parameters to visualize. Works best with \code{keepYlab} and \code{keepRightStrip} set to \code{TRUE}.
#' @param errbartabwidth A numerical vector of the same length as the number of panels in the plot (eror bar panels + table panels). Specifies the relative width of the panels.
#' @param errbarplotscale Scaling factor for the width of the leftmost errorbar plot to compensate for y-axis labels.
#' @param tabplotscale Scaling factor for the width of the rightmost column (usually a table plot) to adjust for the size of the right strip.
#' @param onlySignificantErrorBars Logical. Should error bars be hidden for non-significant covariates (TRUE) or be shown for all covariates regardless of significance (FALSE).
#' @param addcodeErr A string of code to be applied to each of the panels with error bars.
#' @param ... Arguments passed on to \code{ggpubr::text_grob}.
#'
#' @details
#' If \code{plotData} is NULL, \code{dfres} is passed to \code{setupForestPlotData} to create a \code{data.frame} that contains the exact data to be plotted, including the numerical statistics.
#' If \code{return} is set to "data" then the \code{data.frame} is returned. This makes it possible to to review the data and/or do modifications if necessary. \code{forestPlot} can
#' be called with \code{plotData} set to the name of the \code{data.frame}, in which case \code{dfres} is not required and the \code{plotData} \code{data.frame} will be used for
#' creating the Forest plot.
#'
#' If \code{referenceInfo} is set to \code{auto} (the default) then generic information based on \code{referenceParameters} and the \code{REFROW} column in \code{dfres} will be included at the bottom of the plot.
#' The \code{REFROW} column in \code{dfres} will be \code{NO} if \code{dfRefRow} is \code{NULL} in the call to \code{getForestDFSCM}, \code{getForestDFemp} and \code{getForestDFFREM} and \code{YES} if it is not.
#'
#' If \code{REFROW} is \code{NO}   and \code{referenceParameters} is \code{final} then the reference information text will be: "The reference line is based on the final parameter estimates and the reference covariate values in the model.".
#'
#' If \code{REFROW} is \code{NO}   and \code{referenceParameters} is \code{func} then the reference information text will be: "The reference line is based on the average parameter estimates over the posterior parameter distribution and the reference covariate values in the model.".
#'
#' If \code{REFROW} is \code{YES}  and \code{referenceParameters} is \code{final} then the reference information text will be: "The reference line is based on the final parameter estimates and selected covariate values.".
#'
#' If \code{REFROW} is \code{YES}  and \code{referenceParameters} is \code{func}  then the reference information text will be: "The reference line is based on the average parameter estimates over the posterior parameter distribution and selected covariate values.".
#'
#'If \code{referenceInfo} is not \code{NULL} and not \code{auto} then the reference information text will be \code{referenceInfo}.
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
#' forestPlot(dfres,parameters=c("CL"),groupNameLabels = c("Age (y)","Sex","Weight (kg)"),
#'            tabTextSize = 20,xlb="Relative parameter value")
#'
#' # Stack the plots instead of plot them horizontally
#'
#' forestPlot(dfres,parameters = c("CL","Frel"),stackedPlots = TRUE,
#'            keepYlabs = TRUE,keepRightStrip = TRUE)
#' }
forestPlot <- function(dfres,
                       plotData=NULL,
                       plotRelative=TRUE,
                       noVar = TRUE,
                       referenceParameters = "func",
                       sigdigits = 2,
                       parameters=unique(dfres$PARAMETER), #Labels
                       parameterLabels=parameters,
                       parameterLabelsPrefix = NULL,
                       groupNameLabels = NULL,
                       referenceInfo = "auto",
                       labelfun=label_value,
                       groupname_labelfun=label_value,
                       ref_area=c(0.8,1.25),
                       ref_fill_col  = "lightgrey",
                       ref_fill_alpha=0.4,
                       ref_line_size=1,
                       ref_line_type="dotted",
                       ref_line_col="black",
                       ci_line_type="solid",
                       ci_line_col="black",
                       ci_line_size=0.7,
                       point_shape=16,
                       point_color="black",
                       point_size=2.5,
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
                       statisticsLabel = "Statistics:", #paste("Statistics:",parameterLabels),
                       xlb = ifelse(plotRelative,"Relative parameter value","Parameter value"),
                       commonXlab = FALSE,
                       size.legend.text = rel(0.8),
                       return = "plot",
                       table=TRUE,
                       rightStrip=TRUE,
                       errbartabwidth = rep(c(2,1),length(parameters)),
                       errbarplotscale = 1.45,
                       tabplotscale    = 1.1,
                       onlySignificant = FALSE,
                       onlySignificantErrorBars = FALSE,
                       setSignEff = NULL,
                       size=theme_get()$text$size*0.8,
                       addcodeErr="NULL",
                       xlim=c(NA,NA),
                       ...) {


  ## Setup the plotting data frame
  if(is.null(plotData)) {
    plotData<- setupForestPlotData(dfres,
                                   reference             = referenceParameters,
                                   plotRelative          = plotRelative,
                                   parameters            = parameters,
                                   parameterLabels       = parameterLabels,
                                   parameterLabelsPrefix = parameterLabelsPrefix,
                                   groupNameLabels       = groupNameLabels,
                                   statisticsLabel       = statisticsLabel,
                                   noVar                 = noVar,
                                   sigdigits             = sigdigits,
                                   onlySignificant       = onlySignificant,
                                   setSignEff            = setSignEff)
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

    ## Set the limits of the error bars to the point value in case we don't want to have error bars for non-significant covariates.
    if(onlySignificantErrorBars) {
      data <- data %>% mutate(q1 = ifelse(COVEFF,q1,point),
                              q2 = ifelse(COVEFF,q2,point))

    }

    p1a <- ggplot(data,aes(x=point,y=COVNAME,xmin=q1,xmax=q2)) +
      geom_blank() +

      geom_rect(data=rect_data,
                aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="Reference area"),alpha=ref_fill_alpha,
                inherit.aes = FALSE) +

      geom_vline(data           =rect_data,
                 aes(xintercept = ref_value,
                     linetype   = "Reference subject",
                     color      = "Reference subject"
                 ),
                 size           = ref_line_size,
                 key_glyph      = "path") +

      geom_errorbarh(aes(color="CI",linetype="CI"),key_glyph = "path",height=0,size=ci_line_size) +
      geom_point(aes(shape="Point estimate"),color=point_color,size=point_size) +

      #browser()
      scale_fill_manual(name     = NULL, values = c("Reference area" = ref_fill_col),
                        labels=c("Reference area" = ref_area_label)) +
      scale_color_manual(name    = NULL, values = c("Reference subject" = ref_line_col,"CI"=ci_line_col),
                         labels=c("Reference subject" =ref_subj_label,"CI"=ci_label)) +
      scale_linetype_manual(name = NULL, values = c("Reference subject" = ref_line_type,"CI"=ci_line_type),
                            labels=c("Reference subject" =ref_subj_label,"CI"=ci_label)) +
      scale_shape_manual(name    = NULL, values = c("Point estimate"    = point_shape),
                         labels=c("Point estimate" = point_label)) +

      guides(linetype=guide_legend(override.aes=list(linewidth=1))) +
      facet_grid(GROUPNAMELABEL~PARAMETERLABEL,scales = "free",
                 labeller = labeller(PARAMETERLABEL= label_fun,
                                     GROUPNAMELABEL = group_name_label_fun)) +
      ylab(NULL) +
      xlab(xlb) +
      coord_cartesian(xlim=xlim) +
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

    plotList[[i]] <- parPlot(subset(plotData,PARAMETER==as.character(parameters[i])),parameters,parameterLabels,
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

    plotList[[i]] <- eval(parse(text=paste("plotList[[i]]+",addcodeErr)))

  }

  ## The table plots
  for(i in 1:length(parameters)) {

    tabList[[i]] <- tablePlot(subset(plotData,PARAMETER==as.character(parameters[i])),parameters,
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

  ## Remove all xlabs. Will add a common one later
  if(commonXlab) {
    for(i in 1:length(totList)) {
      totList[[i]] <- totList[[i]] + rremove("xlab")
    }
  }

  ## Set the font size of the legend
  totList[[1]] <- totList[[1]] + theme(legend.text = element_text(size=size.legend.text))

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

  ## Add common x-axis label
  if(commonXlab) {
    myPlot <- ggpubr::annotate_figure(myPlot,bottom=ggpubr::text_grob(xlb,...))
  }

  ## Add information about the reference subject if needed
  if(is.null(referenceInfo)) {
    myPlot <- myPlot
  } else if(!is.null(referenceInfo) & referenceInfo != "auto") {
    myPlot <- ggpubr::annotate_figure(myPlot,bottom=ggpubr::text_grob(referenceInfo,size=size,...))
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

    myPlot <- ggpubr::annotate_figure(myPlot,bottom=ggpubr::text_grob(refText,size=size,...))
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
