## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,cache=F------------------------------------------------------------
packageName <- "PMXForest"
suppressPackageStartupMessages(library(packageName,character.only = TRUE))
library(tools)
library(ggplot2)
library(stats)
library(boot)

theme_set(theme_bw(base_size=18))
theme_update(plot.title = element_text(hjust = 0.5))

set.seed(865765)

## -----------------------------------------------------------------------------
library(datasets)
# Create useable (long format) data set Titanic from frequency
dfT<-as.data.frame(Titanic)

#Convert frequencies to individual observations
dfData<-data.frame()
for (i in 1:nrow(dfT)) {
  if (dfT$Freq[i]!=0) dfData<-rbind(dfData,data.frame(Class=as.character(dfT$Class[i]),
                                                      Sex=as.character(dfT$Sex[i]),
                                                    Age=as.character(dfT$Age[i]),
                                                    Survived=as.character(dfT$Survived[i]),
                                                      stringsAsFactors = T)[rep(1,dfT$Freq[i]),])
}


## -----------------------------------------------------------------------------

#A GLM fit function for log regresssion
glmfitfunc <- function(fitdata,indicies,strDV,strCovs,strInteraction = "") {
  dt<-fitdata[indicies,]
  glmfit<-glm(formula = paste0(strDV," ~ ",paste(strCovs,collapse = "+"),strInteraction),
              family = binomial(link = "logit"), 
              data = dt)
  return(glmfit)
}

#A GLM pred function
glmpredfunc <- function(glmfit,preddata) {
  #Predict probability of survived
  return(predict(glmfit,newdata=preddata,type="response")) 
}

#A GLM fit & pred function
glmfitpredfunc <- function(fitdata,indicies,strDV,strCovs,strInteraction = "",preddata=NULL) {
  if (is.null(preddata)) preddata<-fitdata[indicies,]
  return(glmpredfunc(glmfitfunc(fitdata,indicies,strDV,
                                strCovs,strInteraction),
                     preddata))
}

#Fit this dataset
#glmfitfunc(dfData,1:nrow(dfData),strDV = "Survived",strCovs = c("Class","Age","Sex"),strInteraction = "")


#Predict probabilities for this dataset
#glmfitpredfunc(dfData,1:nrow(dfData),strDV = "Survived",strCovs = c("Class","Age","Sex"),strInteraction = "")


## -----------------------------------------------------------------------------
dfCovs <- createInputForestData(
  list( "Class" = c("1st","2nd","3rd","Crew"),
        "Sex" = c("Male","Female"),
        "Age" = c("Child","Adult")),
  iMiss=-99)
dfCovs$Age[dfCovs$Age==-99]<-"Child"
dfCovs$Sex[dfCovs$Sex==-99]<-"Male"
dfCovs$Class[dfCovs$Class==-99]<-"3rd"


## ---- echo=F------------------------------------------------------------------
dfCovs

## -----------------------------------------------------------------------------
covnames <- c("First class","2nd class","3rd class","Crew","Male","Female","Child","Adult")

## -----------------------------------------------------------------------------
dfRefRow <- data.frame("Age"="3rd","Sex"="Male","Class"="Child")

## -----------------------------------------------------------------------------
#Do a bootstrap based prediction of the covariate effect in dfCovs
BoostrapPreds <- boot(dfData,R=500,statistic=glmfitpredfunc,
                      strDV="Survived",strCovs=c("Age","Sex","Class"),strInteraction="",preddata=dfCovs)

#Do predictions using the original fit
orgpred<-glmfitpredfunc(dfData,1:nrow(dfData),strDV = "Survived",strCovs = c("Class","Age","Sex"),strInteraction = "")

dfres<-dfCovs
dfres$COVNAME<-covnames
iGroups<-c(1,1,1,1,5,5,7,7) #Make the grouping of the covariates
dfres$GROUP<-iGroups
dfres$COVNUM<-1:nrow(dfres)
#Set the medianreference (i.e. 3rd combination of dfCovs: Child, Male, 3rd Class)
dfres$REFFUNC=100*median(BoostrapPreds$t[,3])
#Set the orginal data median reference
dfres$REFTRUE<-100*median(orgpred)

for (i in 1:nrow(dfres)) {
  #Get 95% CI and median
  q<-quantile(BoostrapPreds$t[,i],probs=c(0.025,0.5,0.975),names=FALSE)
  dfres$Q1[i]<-100*q[1] #Lower CI
  dfres$POINT[i]<-100*q[2] #Median 
  dfres$Q2[i]<-100*q[3] #Upper CI
}

dfres$GROUPNAME<-dfres$COVARIATEGROUPS
dfres$PARAMETER="Chance of surviving Titanic %"
dfres$COVEFF=F

## ----ForestSurvival,fig.width=10,fig.asp=.5,out.width="100%",cache=F----------
forestPlot(dfres,plotRelative = F,referenceInfo=NULL,xlb = "Reference subject is a boy in 3rd class",errbartabwidth = c(1,1))


