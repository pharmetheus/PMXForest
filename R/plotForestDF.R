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
                        freescalex_wrap=(!plotRelative || addTable))
{

  if (addTable) { #Insert CI table
    df$XMAX<-NA
    df<-plyr::ddply(df,  .variables=c("PARAMETER"), .fun = function(x,pr=plotRelative){
      if (pr) {
        x$XMAX<-max(x$q3/x$REFMEDIAN)
      } else {
        x$XMAX<-max(x$q3)
      }
      return(x)
    })
    if (plotRelative) {
      df$Q1TEXT<-ifelse(df$q1/df$REFMEDIAN<df$q3/df$REFMEDIAN,df$q1/df$REFMEDIAN,df$q3/df$REFMEDIAN)
      df$Q3TEXT<-ifelse(df$q3/df$REFMEDIAN>df$q1/df$REFMEDIAN,df$q3/df$REFMEDIAN,df$q1/df$REFMEDIAN)
    } else {
      df$Q1TEXT<-ifelse(df$q1<df$q3,df$q1,df$q3)
      df$Q3TEXT<-ifelse(df$q3>df$q1,df$q3,df$q1)
    }
  }

  p<-ggplot(df)

  if (plotRelative) {
    #Based on relative  values
    p<-p+geom_vline(aes(xintercept=1),size=vlinesize,color=vlinecol)
    p<-p+geom_errorbarh(aes(y=Y,xmin=q1/REFMEDIAN,xmax=q3/REFMEDIAN,color=as.factor(GROUP)),size=errorbarsize)
    #p<-p+geom_errorbarh(aes(y=Y,xmin=q1/REFMEDIAN,xmax=q3/REFMEDIAN,x=q2/REFMEDIAN,color=as.factor(GROUP)),size=errorbarsize)
    p<-p+geom_point(aes(y=Y,x=q2/REFMEDIAN,color=as.factor(GROUP)),size=pointsize)
  } else {
    #Based on actual values
    #p<-p+geom_vline(aes(xintercept=REFMEDIAN),size=vlinesize,color=vlinecol)
    p<-p+geom_errorbarh(aes(y=Y,xmin=q1,xmax=q3,color=as.factor(GROUP)),size=errorbarsize)
    p<-p+geom_point(aes(y=Y,x=q2,color=as.factor(GROUP)),size=pointsize)
  }

  p<-p + scale_y_continuous(breaks=unique(df$Y),labels = unique(df$COVNAME))
  p<-p+guides(color=FALSE)
  p<-p+theme(panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank(), axis.ticks.y = element_blank())
  p<-p+ylab(strylab)+xlab(strxlab)

  if (addTable) { #Insert CI table

    if (plotRelative) {
      p<-p+geom_text(aes(y=Y,x=XMAX+xtextoffset,label=paste0(formatC(q2/REFMEDIAN,format="f",decdig)," [",formatC(Q1TEXT,format="f",decdig)," - ",formatC(Q3TEXT,format="f",decdig),"]")),hjust = 0,size=textsize)
      p<-p+geom_segment(aes(y=Y,x=q1/REFMEDIAN,yend=Y,xend=q3/REFMEDIAN*percentscalexmax),color="NA")
    } else {
      p<-p+geom_text(aes(y=Y,x=XMAX*xtextoffsetpercent,label=paste0(formatC(q2,format="f",decdig)," [",formatC(Q1TEXT,format="f",decdig)," - ",formatC(Q3TEXT,format="f",decdig),"]")),hjust = 0,size=textsize)
      p<-p+geom_segment(aes(y=Y,x=q1,yend=Y,xend=q3*percentscalexmax),color="NA")
    }
  }

  if (freescalex_wrap) {
    p<-p+facet_wrap(~PARAMETER,scales="free_x")
  } else {
    p<-p+facet_wrap(~PARAMETER)
  }

  return(p)
}
