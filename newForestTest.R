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

  return(list(CL,FREL))
}

functionListName <- c("CL","FREL")


runno   <- 1
extFile <- system.file("extdata",paste0("run",runno,".ext"),package="PMXForest")
covFile <- system.file("extdata",paste0("run",runno,".cov"),package="PMXForest")
dfSamplesCOV <- getSamples(covFile,extFile=extFile,n=200)



setupPlotData <- function(dfCovs,covnames,functionList,functionListName,noBaseThetas,dfSamples,plotRelative=TRUE,noVar=FALSE,reference="func") {

  vars <- getPlotVars(plotRelative,noVar,reference)

  dfres <- getForestDFSCM(dfCovs           = dfCovs,
                          cdfCovsNames     = covnames,
                          functionList     = functionList,
                          functionListName = functionListName,
                          noBaseThetas     = noBaseThetas,
                          dfParameters     = dfSamples
  )

  sigdigits <- 2
  plotData <- dfres %>% select(GROUPNAME,COVNUM,COVNAME,PARAMETER,REF=vars["ref"],point=vars["point"],q1=vars["q1"],q2=vars["q2"]) %>%
    mutate(
      reference=reference,
      REF=ifelse(plotRelative,1,REF)) %>%
    mutate(COVNAME = factor(COVNUM,labels=unique(COVNAME))) %>%
    mutate(
      meanlabel  = signif_pad(REF, sigdigits),
      lowcilabel = signif_pad(q1, sigdigits),
      upcilabel  = signif_pad(q2,sigdigits),
      LABEL      = paste0(meanlabel, " [", lowcilabel, "-", upcilabel,"]"),
      LABEL      = str_pad(LABEL,max(str_length(LABEL)),side="right",' ')
    )

  return(plotData)
}

plotData<- setupPlotData(dfCovs,covnames,list(paramFunction),functionListName,noBaseThetas = 16,dfSamplesCOV)


# ## Adjust dfres
# dfres <- dfres %>% mutate(COVNAME = factor(COVNUM,labels=unique(COVNAME)))
# #dfres <- dfres %>% mutate(GROUP = factor(GROUP,labels=covariate))
#
# sigdigits <- 2
# dfres <- dfres %>% mutate(
#   meanlabel  = signif_pad(POINT_REL_REFFUNC, sigdigits),
#   lowcilabel = signif_pad(Q1_REL_REFFUNC, sigdigits),
#   upcilabel  = signif_pad(Q2_REL_REFFUNC,sigdigits),
#   LABEL      = paste0(meanlabel, " [", lowcilabel, "-", upcilabel,"]"),
#   LABEL      = str_pad(LABEL,max(str_length(LABEL)),side="right",' ')
#   )

## Main plot

ref_value <- 1
ref_area  <- c(0.8,1.2)
ref_fill_col <- "gray"


p1a <- ggplot(subset(plotData,PARAMETER=="CL"),aes(x=point,y=COVNAME,xmin=q1,xmax=q2)) +
  theme_bw() +
  geom_blank() +
  annotate("rect", xmin = ref_value * min(ref_area), xmax = ref_value * max(ref_area),
           ymin = -Inf, ymax = Inf, fill = ref_fill_col,alpha=0.5) +

  geom_ribbon(data = data.frame(x = ref_value, ymax = max(ref_area), ymin = min(ref_area),fill = "Reference area"),
              aes(x = x, ymax = ymax, ymin = ymin, fill=fill),
              size = 1,inherit.aes = FALSE) +

  geom_vline(aes(xintercept = REF, linetype = "Reference subject",
               color = "Reference subject"), size = 1,key_glyph = "path") +

  geom_errorbarh(aes(color="CI",linetype="CI"),key_glyph = "path",height=0) +
  geom_point(aes(color="CI",shape="Point estimate"),size=2.5) +

  scale_linetype_manual(name=NULL,values = c("CI"="solid","Reference subject"="dotted")) +
  scale_fill_manual(name=NULL,values = c("Reference area"="gray")) +
  scale_color_manual(name=NULL,values = c("CI"="blue","Reference subject"="black")) +
  scale_shape_manual(name=NULL,values = c("Point estimate"=16)) +

  guides(shape=guide_legend(override.aes = list(color=c("blue"))),
         # linetype=guide_legend(override.aes=list(linetype="Reference subject")),
         color=guide_legend(override.aes = list(shape=NA))
          ) +
  facet_grid(GROUPNAME~PARAMETER,scales = "free",labeller = label_wrap_gen(width=5),space = "free") +
  ylab(NULL)+

  theme(strip.text.y=element_blank()) +
  theme(plot.margin = unit(c(5.5,0,5.5,5.5), "pt"))


p1b <- ggplot(subset(plotData,PARAMETER=="FREL"),aes(x=point,y=COVNAME,xmin=q1,xmax=q2)) +
  theme_bw() +
  geom_blank() +
  annotate("rect", xmin = ref_value * min(ref_area), xmax = ref_value * max(ref_area),
           ymin = -Inf, ymax = Inf, fill = ref_fill_col,alpha=0.5) +

  geom_ribbon(data = data.frame(x = ref_value, ymax = max(ref_area), ymin = min(ref_area),fill = "Reference area"),
              aes(x = x, ymax = ymax, ymin = ymin, fill=fill),
              size = 1,inherit.aes = FALSE) +

  geom_vline(aes(xintercept = REF, linetype = "Reference subject",
                 color = "Reference subject"), size = 1,key_glyph = "path") +

  geom_errorbarh(aes(color="CI",linetype="CI"),key_glyph = "path",height=0) +
  geom_point(aes(color="CI",shape="Point estimate"),size=2.5) +

  scale_linetype_manual(name=NULL,values = c("CI"="solid","Reference subject"="dotted")) +
  scale_fill_manual(name=NULL,values = c("Reference area"="gray")) +
  scale_color_manual(name=NULL,values = c("CI"="blue","Reference subject"="black")) +
  scale_shape_manual(name=NULL,values = c("Point estimate"=16)) +

  guides(shape=guide_legend(override.aes = list(color=c("blue"))),
         # linetype=guide_legend(override.aes=list(linetype="Reference subject")),
         color=guide_legend(override.aes = list(shape=NA))
  ) +
  facet_grid(GROUPNAME~PARAMETER,scales = "free",labeller = label_wrap_gen(width=5),space = "free") +
  theme(strip.text.y=element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.title.y = element_blank())


## Stats plot

#p2a <- ggplot(subset(dfres,PARAMETER=="CL"),aes(x=1,y=COVNAME)) +

parlabs <- c("CL statistics","FREL statistics")
names(parlabs) <- c("CL","FREL")

p2a <- ggplot(plotData %>% filter(PARAMETER=="CL"),aes(x=1,y=COVNAME)) +
  geom_text(aes(label = LABEL)) +
  facet_grid(GROUPNAME~PARAMETER,scales = "free",labeller = labeller(PARAMETER=parlabs),space = "free") +
  theme_bw() +
  #scale_x_continuous(expand=expand_scale(mult=c(0.01,0.02)))+
  #theme(panel.spacing.x = unit(0, "lines"))+
  theme(plot.margin = unit(c(5.5,5.5,5.5,0), "pt")) +
  theme(strip.text.y=element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank())

p2b <- ggplot(subset(plotData,PARAMETER=="FREL"),aes(x=1,y=COVNAME)) +
  geom_text(aes(label = LABEL)) +
  facet_grid(GROUPNAME~PARAMETER,scales = "free",labeller = labeller(PARAMETER=parlabs),space = "free") +
  theme_bw() +
  # theme(strip.text.x=element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank())


ggpubr::ggarrange(p1a,p2a,p1b,p2b,nrow=1,widths =rep(c(1,0.5),2),align="h",common.legend = T)



