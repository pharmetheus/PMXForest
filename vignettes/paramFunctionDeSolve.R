
#This function calculates concentrations and AUC following a single dose with the deSolve package.
paramFunctionDeSolve <- function(thetas, df, times, ...) {

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
  S2    = V

  ## parameters for the ODE
  pars <- c(
    KA = KA,
    DOSE = 80*F1,
    D1 = D1,
    CL  = CL,
    V  = V)

  #Inline function of the ODE, note that an extra
  #compartment for AUC is added
  Fun_ODE <- function (t, y, pars) {
    with (as.list(c(y, pars)), {
      r<-0 #Input rate into depot compartment
      if (t<=D1) r<-DOSE/D1
      dA1<- -KA*A1 + r
      dA2<-  KA*A1 - CL/V*A2
      dA3<-  A2*1000/V
      return(list(dy = c(dA1, dA2,dA3),
                  CONC = A2*1000/V,
                  AUC  = A3))
    })}

  #Initial conditions
  y <- c(A1 = 0, A2 = 0, A3 = 0)
  #Sort times if they come in non-increasing order
  s1<-sort(times,index.return=TRUE)
  ## ODE model solved with ode solver
  ODE <- ode(y = y, times = s1$x, func = Fun_ODE,
             parms = pars, atol = 1e-10, rtol = 1e-10)

  #Return all compartments
  return(ODE[s1$ix,])
}
