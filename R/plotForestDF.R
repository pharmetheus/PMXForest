#' plotForestDF()
#'
#' @param df the forest plot data frame, should be based on a call to getForestDF or manually created accfording to the getForestDF structure
#' @param plotRelative if the plot should be based on the relative scale to the reference (or parameter scale), default = TRUE
#' @param addTable if a CI table should be added to the plot, default = TRUE
#' @param strxlab x-label, default = "Covariate-Parameter effect"
#' @param strylab y-label, default = ""
#' @param decdig number of digits after the comma in the CI table, only used with addTable = TRUE, default=2
#' @param percentscalexmax increase of xmax, default 1.5 (fraction of x-axis), only used with addTable = TRUE
#' @param xtextoffset textoffset, default 0.25 on x-axis scale, only used with addTable = TRUE
#' @param xtextoffsetpercent default = 1.1, scale the xaxis by the fraction, only used if addTable = TRUE
#' @param textsize default = 6, text size of the CI text, only used if addTable = TRUE
#' @param vlinecol color of vertical reference line, default "lightgrey"
#' @param vlinesize size of vertical reference line, default 1.5
#' @param pointsize size of geom point (e.g. median value), default 3
#' @param errorbarsize size of errorbars, default 1
#' @param freescalex_wrap if wrap as free_x or not, default based on plotRelative and addTable
#' @param log10x if the x-axis should be set to a log10 axis, default=FALSE
#' @param fixedSpacing A boolean (TRUE/FALSE) if fixed spacing between covariate groups should be used in the Forest plot. Default is TRUE.
#' If FALSE, the y coordinates are calculated relative to the number of groups and numbers of covariates within a group.
#' If fixed spacing is used, groupdist and withingroupdist will be used as well.
#' @param groupeddist A number defining the y distance between groups of covariates.
#' @param withingroupeddist A number defining the y distance within groups of covariates.
#' @param useTrueRef A flag representing which value to use as the reference. Per default set to FALSE,
#' indicating that the Ref function used is the same as function used to calculate the "point" for each covariate line. If set to TRUE,
#' the True reference, i.e. the first sample in the dfSamples, which is often set to the final model estimates, will be used.
#' @param useRefUncertainty A flag representing if the rleative reference (no covariate effect) should take the parameter uncertainty into account. Per default set to TRUE,
#' indicating that the relative forest plot uncertainty is also accounting for the structural parameter uncertainty. If set to FALSE,
#' the relative reference value as well as the uncertainty for a covariate without an effect will be 1.


#'
#' @return a ggplot2 object with a Forest plot
#' @export
#'
#' @examples
#' #' \dontrun{
#' plotForestDF(dfresCOV,textsize = 5) +
#' xlim(0,3.3) +
#' geom_vline(xintercept=c(0.8,1.25),linetype="dashed")
#' }
plotForestDF <-function(df,
                        plotRelative=TRUE,
                        addTable=TRUE,
                        strxlab="Covariate-Parameter effect",
                        strylab="",
                        decdig=2,
                        percentscalexmax=1.5,
                        xtextoffset=0.25,
                        xtextoffsetpercent=1.1,
                        textsize=5,
                        vlinecol="lightgrey",
                        vlinesize=1.5,
                        pointsize=3,
                        errorbarsize=1,
                        freescalex_wrap=(!plotRelative || addTable),
                        log10x=FALSE,
                        fixedSpacing = TRUE,
                        groupdist = 0.3,
                        withingroupdist = 0.2,
                        useTrueRef=FALSE,
                        useRefUncertainty=TRUE)
{

  if (plotRelative && useTrueRef && useRefUncertainty==FALSE) {
    warning("True ref cannot be mixed with RefUncertainty=FALSE, \nRefUncertainty is automatically set to TRUE")
    useRefUncertainty<-TRUE
  }

  #Handle the spacing in Y based on the GROUP of covariates
  if (!fixedSpacing) {
    df$Y <- NA
    for (i in 1:length(sort(unique(df$GROUP)))) {
      dft <- df[df$GROUP == sort(unique(df$GROUP))[i], ]
      num_in_group <- nrow(dft) / length(unique(df$PARAMETER))
      for (n in 1:length(unique(df$COVNUM))) {
        df$Y[df$GROUP == sort(unique(df$GROUP))[i] & df$COVNUM == unique(dft$COVNUM)[n]] <- (i -1) + n / (num_in_group + 1)
      }
    }
  }  else {
    df$Y <- NA
    a <- 0
    for (i in 1:length(sort(unique(df$GROUP)))) {
      dft <- df[df$GROUP == sort(unique(df$GROUP))[i], ]
      num_in_group <- nrow(dft) / length(unique(df$PARAMETER))
      for (n in 1:length(unique(dft$COVNUM))) {
        df$Y[df$GROUP == sort(unique(df$GROUP))[i] & df$COVNUM == unique(dft$COVNUM)[n]] <- a
        a <- a + withingroupdist
      }
      a <- a + groupdist
    }
  }

  #Inform on what to use as a reference, TRUE (=final model estimates) or the function used in getForestDF....
  strRefVar<-"REFFUNC"
  if (useTrueRef) strRefVar<-"REFTRUE"


  if (addTable) { #Insert CI table
    df$XMAX<-NA
    df<-plyr::ddply(df,  .variables=c("PARAMETER"), .fun = function(x,pr=plotRelative,
                                                                    ru=useRefUncertainty,
                                                                    strRV=strRefVar){
      if (pr) {
        if (ru) {
          x$XMAX<-max(x$q2/x[[strRV]])
        } else {
          x$XMAX<-max(x$q2_NOVAR)
        }

      } else {
        x$XMAX<-max(x$q2)
      }
      return(x)
    })
    if (plotRelative) {
      if (useRefUncertainty) {
        df$Q1TEXT<-ifelse(df$q1/df[[strRefVar]]<df$q2/df[[strRefVar]],df$q1/df[[strRefVar]],df$q2/df[[strRefVar]])
        df$Q2TEXT<-ifelse(df$q2/df[[strRefVar]]>df$q1/df[[strRefVar]],df$q2/df[[strRefVar]],df$q1/df[[strRefVar]])
      } else {
        df$Q1TEXT<-ifelse(df$q1_NOVAR<df$q2_NOVAR,df$q1_NOVAR,df$q2_NOVAR)
        df$Q2TEXT<-ifelse(df$q2_NOVAR>df$q1_NOVAR,df$q2_NOVAR,df$q1_NOVAR)
      }
    } else {
      df$Q1TEXT<-ifelse(df$q1<df$q2,df$q1,df$q2)
      df$Q2TEXT<-ifelse(df$q2>df$q1,df$q2,df$q1)
    }
  }

  p<-ggplot(df)

  if (plotRelative) {
    #Based on relative  values
    refcol<-sym(strRefVar)
    p<-p+geom_vline(aes(xintercept=1),size=vlinesize,color=vlinecol)
    if (useRefUncertainty) { #If uncertainty for the reference
      p<-p+geom_errorbarh(aes(y=Y,xmin=q1/!!refcol,xmax=q2/!!refcol,color=as.factor(GROUP)),size=errorbarsize)
      p<-p+geom_point(aes(y=Y,x=FUNC/!!refcol,color=as.factor(GROUP)),size=pointsize)
    } else { #Without uncertainty for the reference
      p<-p+geom_errorbarh(aes(y=Y,xmin=q1_NOVAR,xmax=q2_NOVAR,color=as.factor(GROUP)),size=errorbarsize)
      p<-p+geom_point(aes(y=Y,x=FUNC_NOVAR,color=as.factor(GROUP)),size=pointsize)
    }
  } else {
    #Based on actual values
    p<-p+geom_errorbarh(aes(y=Y,xmin=q1,xmax=q2,color=as.factor(GROUP)),size=errorbarsize)
    p<-p+geom_point(aes(y=Y,x=FUNC,color=as.factor(GROUP)),size=pointsize)
  }

  p<-p + scale_y_continuous(breaks=unique(df$Y),labels = unique(df$COVNAME))
  p<-p+guides(color=FALSE)
  p<-p+theme(panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank(), axis.ticks.y = element_blank())
  p<-p+ylab(strylab)+xlab(strxlab)
  if (log10x) p<-p+scale_x_log10()

  if (addTable) { #Insert CI table
    if (plotRelative) {
        refcol<-sym(strRefVar)
        if (useRefUncertainty) {
          p<-p+geom_text(aes(y=Y,x=XMAX+xtextoffset,label=paste0(formatC(FUNC/!!refcol,format="f",decdig)," [",formatC(Q1TEXT,format="f",decdig)," - ",formatC(Q2TEXT,format="f",decdig),"]")),hjust = 0,size=textsize)
          p<-p+geom_segment(aes(y=Y,x=q1/!!refcol,yend=Y,xend=q2/!!refcol*percentscalexmax),color="NA")
        } else {
          p<-p+geom_text(aes(y=Y,x=XMAX+xtextoffset,label=paste0(formatC(FUNC_NOVAR,format="f",decdig)," [",formatC(Q1TEXT,format="f",decdig)," - ",formatC(Q2TEXT,format="f",decdig),"]")),hjust = 0,size=textsize)
          p<-p+geom_segment(aes(y=Y,x=q1_NOVAR,yend=Y,xend=q2_NOVAR*percentscalexmax),color="NA")
        }
    } else {
      p<-p+geom_text(aes(y=Y,x=XMAX*xtextoffsetpercent,label=paste0(formatC(FUNC,format="f",decdig)," [",formatC(Q1TEXT,format="f",decdig)," - ",formatC(Q2TEXT,format="f",decdig),"]")),hjust = 0,size=textsize)
      p<-p+geom_segment(aes(y=Y,x=q1,yend=Y,xend=q2*percentscalexmax),color="NA")
    }
  }

  if (freescalex_wrap) {
    p<-p+facet_wrap(~PARAMETER,scales="free_x")
  } else {
    p<-p+facet_wrap(~PARAMETER)
  }

  return(p)
}
