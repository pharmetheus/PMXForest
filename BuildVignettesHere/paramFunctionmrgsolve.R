#This function calculates the concentration at any given times as well as the AUC
# assuming multiple doses with a tau and starting at time 0. Using a model defined using mrgsolve syntax
paramFunctionmrgsolve <- function(thetas, df, times, tau,model, ...) {

  #Load the model again since a different R process is used using parallelization
  loadso(model)

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

  D1FR <- TVD1
  MAT  <- TVMAT

  D1    = MAT*(1-D1FR)
  F1    = FREL
  KA    = 1 / (MAT-D1)

  dose = 80*F1 #Calculate the actual input
  rate = dose/D1 #Calculate the actual rate

  sampletimes=times

  dosetimes<-seq(0,max(sampletimes),by=tau)
  amt<-rep(dose,length(dosetimes))
  rate=rep(rate,length(dosetimes))

  #Create the design dataset
  data <-
    tibble::tibble(ID = 1,
                   cmt=1,
                   time=c(sampletimes,dosetimes),
                   evid=c(rep(0,length(sampletimes)),rep(1,length(dosetimes))),
                   amt=c(rep(0,length(sampletimes)),amt),
                   rate=c(rep(0,length(sampletimes)),rate),
                   CL=CL,
                   KA=KA,
                   V=V)

  data<-data[order(data$time),]

  #Simulate data using the model and dataset
  out<-mrgsim_q(x=model,data=data)
  tmp <-  data.frame(CONC=out$CP,AUC=out$AUC)
  #Only return the smapletimes, not the dosetimes
  return(tmp[match(sampletimes,out$time),])


}

