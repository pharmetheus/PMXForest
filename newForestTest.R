# library(PMXForest)
#library(ggplot2)
# library(table1) # To get access to signif_pad
# library(ggpubr)

set.seed(865765)

dfCovs <- createInputForestData(
  list( "FORM" = c(0,1),
        "FOOD" = c(0,1),
        "GENO" = c(1,2,3,4),
        "RACEL"= c(1,2,3),
        "WT"   = c(65,115),
        "AGE"  = c(27,62),
        "CRCL" = c(83,150),
        "SEX"  = c(1,2)),
  iMiss=-99)

cGrouping <- c(1,1,2,2,3,3,3,3,4,4,4,5,5,6,6,7,7,8,8)
covnames  <- c("Oral tablets","FDC","Fasted","Fed","2D6 UM","2D6 EM","2D6 IM","2D6 PM","Caucasian",
               "African American","Asian and other","WT 65 kg","WT 115 kg",
               "Age 27 y","Age 62 y","CRCL 83 mL/min","CRCL 150 mL/min","Male","Female")

## The print names of the covariates
covariate <- c("Formulation","Food status","2D6 genotype","Race","Weight","Age","Createnine clearance","Sex")

paramFunction <- function(thetas, df, ...) {

  FRELNCIL <- 1
  if (any(names(df) == "NCIL") && df$NCIL == 1) FRELNCIL <- 1 + thetas[16]

  FRELFORM <- 1
  if (any(names(df) == "FORM") && df$FORM == 0) FRELFORM <- 1 + thetas[15]

  FRELCOV <- FRELFORM * FRELNCIL

  CLFOOD <- 1
  if (df$FOOD == 0) CLFOOD <- 1 + thetas[14]

  CLCOV <- CLFOOD

  TVFREL <- thetas[1]
  if(any(names(df) == "GENO")) {
    if(df$GENO == 1) TVFREL <- TVFREL * (1 + thetas[11])
    if(df$GENO == 3) TVFREL <- TVFREL * (1 + thetas[12])
    if(df$GENO == 4) TVFREL <- TVFREL * (1 + thetas[13])
  }

  TVFREL <- FRELCOV * TVFREL

  if (!any(names(df) == "WT") || df$WT == -99) {
    TVCL <- thetas[4]
  } else {
    TVCL <- thetas[4] * (df$WT / 75)**thetas[2]
  }

  if(any(names(df) == "GENO")) {
    if (df$GENO == 1) TVCL <- TVCL * (1 + thetas[8])
    if (df$GENO == 3) TVCL <- TVCL * (1 + thetas[9])
    if (df$GENO == 4) TVCL <- TVCL * (1 + thetas[10])
  }

  TVCL <- CLCOV * TVCL

  if (!any(names(df) == "WT") || df$WT == -99) {
    TVV <- thetas[5]
  } else {
    TVV <- thetas[5] * (df$WT / 75)**thetas[3]
  }
  TVMAT <- thetas[6]
  TVD1  <- thetas[7]

  FREL <- TVFREL
  CL <- TVCL
  V <- TVV

  return(list(CL,FREL,V))
}

functionListName <- c("CL","FREL","V")


runno   <- 1
extFile <- system.file("extdata",paste0("run",runno,".ext"),package="PMXForest")
covFile <- system.file("extdata",paste0("run",runno,".cov"),package="PMXForest")
dfSamplesCOV <- getSamples(covFile,extFile=extFile,n=200)


dfres <- getForestDFSCM(dfCovs           = dfCovs,
                        cdfCovsNames     = covnames,
                        functionList     = list(paramFunction),
                        functionListName = functionListName,
                        noBaseThetas     = noBaseThetas,
                        dfParameters     = dfSamplesCOV,
                        dfRefRow         = NULL
)


refRow <- dfCovs %>% mutate(
  FORM = -99,
  FOOD = -99,
  GENO = -99,
  RACEL= -99,
  WT   =- 99,
  AGE  = -99,
  CRCL = -99,
  SEX = -99) %>%
  group_by(COVARIATEGROUPS) %>%
  mutate(WT=sample(c(50,75,100),1))

dfresRefRow <- getForestDFSCM(dfCovs           = dfCovs,
                              cdfCovsNames     = covnames,
                              functionList     = list(paramFunction),
                              functionListName = functionListName,
                              noBaseThetas     = noBaseThetas,
                              dfParameters     = dfSamplesCOV,
                              dfRefRow         = refRow
)


plotData<- setupForestPlotData(dfres,plotRelative=FALSE,noVar=TRUE)
plotDataRefRow<- setupForestPlotData(dfresRefRow,plotRelative=FALSE,noVar=TRUE,sigdig=3)




old <- theme_set(theme_bw())

forestPlot <- function(plotData,
                       parameters=unique(plotData$PARAMETER),
                       ref_area=c(0.8,1.2),
                       ref_fill_col="gray",
                       ref_fill_alpha=0.5,
                       ref_line_size=1,
                       ref_line_type="dotted",
                       ref_line_col="black",
                       ci_line_type="solid",
                       ci_line_col="blue",
                       point_shape=16,
                       point_color="blue",
                       ref_subj_label = "Reference subject",
                       ref_area_label = "Reference area",
                       point_label    = "Point estimate",
                       ci_label       = "Confidence interval",
                       statlabs = paste(parameters,"statistics"),
                       xlb = "Parameter value",
                       returnPlotList = FALSE,
                       table=TRUE,
                       rightStrip=TRUE) {


  #######
  ## Create two functions, one for the error bar charts and one for the tables
  #######
  parPlot <- function(data,...) {

    ## This function is designed to only deal with one parameter at a time.
    if(length(unique(data$PARAMETER)) !=1 ) stop("Can only deal with one oparameter at a time.")


    ref_value <- data %>% distinct(GROUPNAME,PARAMETER,REF)
    rect_data <- data.frame(ref_value=ref_value$REF,
                            xmin=ref_value$REF * min(ref_area),
                            xmax = ref_value$REF * max(ref_area),
                            ymin = -Inf, ymax = Inf,
                            GROUPNAME=ref_value$GROUPNAME,
                            PARAMETER=ref_value$PARAMETER)


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

      geom_errorbarh(aes(color="CI",linetype="CI"),key_glyph = "path",height=0) +
      geom_point(aes(shape="Point estimate"),color=point_color,size=2.5) +

      scale_fill_manual(name     = NULL, values = c("Reference area"    = ref_fill_col),
                        labels=c(ref_area_label)) +
      scale_color_manual(name    = NULL, values = c("Reference subject" = ref_line_col,"CI"=ci_line_col),
                         labels=c(ci_label,ref_subj_label)) +
      scale_linetype_manual(name = NULL, values = c("Reference subject" = ref_line_type,"CI"=ci_line_type),
                            labels=c(ci_label,ref_subj_label)) +
      scale_shape_manual(name    = NULL, values = c("Point estimate"    = point_shape),
                         labels=c(point_label)) +

      guides(linetype=guide_legend(override.aes=list(size=1))) +
      facet_grid(GROUPNAME~PARAMETER,scales = "free",labeller = label_wrap_gen(width=5),space = "free") +
      ylab(NULL) +
      xlab(xlb) +
      theme(plot.margin = unit(c(5.5,0,5.5,5.5), "pt"))
  }

  names(statlabs) <- parameters

  tablePlot <- function(data,...) {

    p2a <- ggplot(data,aes(x=1,y=COVNAME)) +
      geom_text(aes(label = LABEL)) +
      facet_grid(GROUPNAME~PARAMETER,scales = "free",labeller = labeller(PARAMETER=statlabs),space = "free") +
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

    plotList[[i]] <- parPlot(subset(plotData,PARAMETER==parameters[i])) +
      theme(plot.margin = unit(c(5.5,0,5.5,5.5), "pt"))

    ## Deal with y-axis elements
    if(i==1 && i != length(parameters)) {
      plotList[[i]] <- plotList[[i]]
    }

    if(i <= length(parameters) && i!=1) {
        plotList[[i]] <- plotList[[i]] +
          theme(axis.text.y = element_blank()) +
          theme(axis.ticks.y = element_blank()) +
          theme(axis.title.y = element_blank())
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
        plotList[[i]] <- plotList[[i]]
      }
    }

  }

  ## The table plots
  for(i in 1:length(parameters)) {

    tabList[[i]] <- tablePlot(subset(plotData,PARAMETER==parameters[i])) +
      theme(plot.margin = unit(c(5.5,5.5,5.5,0), "pt"))

    if(i < length(parameters)) {
      tabList[[i]] <- tabList[[i]] +
        theme(strip.text.y=element_blank())
    }

    if(i == length(parameters) & rightStrip) {  # Last panel with right strip
      tabList[[i]] <- tabList[[i]]
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
    myPlot <- ggarrange(plotlist=totList,nrow=1,widths =rep(c(1,0.5),length(parameters)),align="h",common.legend = T)
  } else{
    myPlot <- ggarrange(plotlist=totList,nrow=1,align="h",common.legend = T)
  }

  ## return either the arranged plot or the list of plot objects
  if(!returnPlotList) {
    return(myPlot)
  } else {
    return(totList)
  }
}


forestPlot(plotData)
forestPlot(plotDataRefRow)
forestPlot(plotData,parameters="CL")
forestPlot(plotData2,parameters="CL")




