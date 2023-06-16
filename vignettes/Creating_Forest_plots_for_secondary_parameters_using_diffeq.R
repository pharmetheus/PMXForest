## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------

packageName <- "PMXForest"
suppressPackageStartupMessages(library(packageName,character.only = TRUE))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(ggplot2))

# Used to enter diff eq models with "easy" designs
suppressPackageStartupMessages(library(deSolve))
# Used to enter diff eq with e.g. multiple dosing
suppressPackageStartupMessages(library(mrgsolve))



theme_set(theme_bw(base_size = 18))
theme_update(plot.title = element_text(hjust = 0.5))

set.seed(865765)

## -----------------------------------------------------------------------------
dfCovs <- createInputForestData(
  list(
    "FORM" = c(0, 1),
    "FOOD" = c(0, 1),
    "GENO" = c(1, 2, 3, 4),
    "RACEL" = c(1, 2, 3),
    "WT" = c(65, 115),
    "AGE" = c(27, 62),
    "CRCL" = c(83, 150),
    "SEX" = c(1, 2)
  ),
  iMiss = -99
)

## -----------------------------------------------------------------------------
covnames <- c(
  "Oral tablets", "FDC", "Fasted", "Fed", "2D6 UM", "2D6 EM", "2D6 IM", "2D6 PM", "Caucasian",
  "African American", "Asian and other", "WT 65 kg", "WT 115 kg",
  "Age 27 y", "Age 62 y", "CRCL 83 mL/min", "CRCL 150 mL/min", "Male", "Female"
)

## -----------------------------------------------------------------------------
dfRefRow <- data.frame("FORM" = 1, "FOOD" = 1, RACEL = 0, "GENO" = 2, "WT" = 75, "AGE" = 70, "CRCL" = 118, "SEX" = 1)

## -----------------------------------------------------------------------------

# This function calculates the concentration at any given time as well as the AUC
# assuming a single dose administered at time 0.
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
  if (any(names(df) == "GENO")) {
    if (df$GENO == 1) TVFREL <- TVFREL * (1 + thetas[11])
    if (df$GENO == 3) TVFREL <- TVFREL * (1 + thetas[12])
    if (df$GENO == 4) TVFREL <- TVFREL * (1 + thetas[13])
  }

  TVFREL <- FRELCOV * TVFREL

  if (!any(names(df) == "WT") || df$WT == -99) {
    TVCL <- thetas[4]
  } else {
    TVCL <- thetas[4] * (df$WT / 75)**thetas[2]
  }

  if (any(names(df) == "GENO")) {
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
  TVD1 <- thetas[7]

  FREL <- TVFREL
  CL <- TVCL
  V <- TVV

  D1FR <- TVD1
  MAT <- TVMAT

  D1 <- MAT * (1 - D1FR)
  F1 <- FREL
  KA <- 1 / (MAT - D1)

  ## parameters for the ODE
  pars <- c(
    KA = KA,
    DOSE = 80 * F1,
    D1 = D1,
    CL = CL,
    V = V
  )

  # Inline function of the ODE, note that an extra
  # compartment for AUC is added
  Fun_ODE <- function(t, y, pars) {
    with(as.list(c(y, pars)), {
      r <- 0 # Input rate into depot compartment
      if (t <= D1) r <- DOSE / D1
      dA1 <- -KA * A1 + r
      dA2 <- KA * A1 - CL / V * A2
      dA3 <- A2 * 1000 / V
      return(list(
        dy = c(dA1, dA2, dA3),
        CONC = A2 * 1000 / V,
        AUC = A3
      ))
    })
  }

  # Initial conditions
  y <- c(A1 = 0, A2 = 0, A3 = 0)
  # Sort times if they come in non-increasing order
  s1 <- sort(times, index.return = TRUE)
  ## ODE model solved with ode solver
  ODE <- ode(
    y = y, times = s1$x, func = Fun_ODE,
    parms = pars, atol = 1e-10, rtol = 1e-10
  )

  # Return all compartments
  return(ODE[s1$ix, ])
}

## -----------------------------------------------------------------------------
# Function to get the cmax within an interval [t1,t2]
Cmax <- function(thetas, df, interval) {
  source("./paramFunctionDeSolve.R")
  wrapper <- function(time, thetas, df) {
    return(paramFunctionDeSolve(thetas, df, times = c(0, time))[2, "CONC"])
  }
  return(optimize(
    f = wrapper, interval = interval,
    maximum = TRUE, thetas = thetas, df = df
  )$objective)
}

# Function to get the auc within an interval [t1,t2]
AUC <- function(thetas, df, interval) {
  source("./paramFunctionDeSolve.R")
  # Add a zero to ensure dosing starts at time 0
  sol <- paramFunctionDeSolve(thetas, df, times = c(0, interval))
  return(sol[nrow(sol), "AUC"] - sol[nrow(sol) - 1, "AUC"])
}

## ----TypicalPlotdeSolve,fig.width=12,fig.asp=0.6,out.width="100%"-------------
runno <- 1
extFile <- system.file("extdata/oldExample", paste0("run", runno, ".ext"), package = packageName)

dfExt <- subset(getExt(extFile = extFile), ITERATION == "-1000000000")

# Call the concentration time-profile function using a Ref covariate setting
sol <- as.data.frame(paramFunctionDeSolve(thetas = as.numeric(dfExt[2:17]), df = dfRefRow, times = seq(0, 72, 0.1)))


p <- ggplot() + geom_line(aes(x = sol$time, y = sol$CONC), size = 1.8) + xlab("Time (h)") + ylab("Concentration ng/mL")
p

## -----------------------------------------------------------------------------
cmax <- Cmax(thetas = as.numeric(dfExt[2:17]), df = dfRefRow, interval = seq(0, 72))

## -----------------------------------------------------------------------------
auc <- AUC(thetas = as.numeric(dfExt[2:17]), df = dfRefRow, interval = c(0, 72))

## -----------------------------------------------------------------------------
bootFile <- system.file("extdata/oldExample", "bs1_n100.dir", paste0("raw_results_run", runno, ".csv"), package = packageName)
dfSamplesBS <- getSamples(bootFile, extFile = extFile)

## ----deSolveCmaxAuc, cache=TRUE-----------------------------------------------
dfres <- getForestDFSCM(
  dfCovs = dfCovs,
  cdfCovsNames = covnames,
  functionList = list(Cmax, AUC),
  functionListName = c("CMAX", "AUC"),
  noBaseThetas = 16,
  dfParameters = dfSamplesBS,
  dfRefRow = dfRefRow,
  cstrPackages = c("deSolve","dplyr"),
  ncores = 1,
  interval = c(0, 24)
)

## ----ForestCOV,fig.width=12,fig.asp=0.6,out.width="100%"----------------------
plotList<-forestPlot(dfres,plotRelative = TRUE,keepYlabs = T,rightStrip = F,return="plotList")
#Add the log10 x-axis
plotList[[1]]<-plotList[[1]]+scale_x_continuous(trans='log10')

ggpubr::ggarrange(plotlist      = plotList,
                  ncol          = 2,
                  nrow          = 2,
                  widths        = c(1,0.5,
                                    1,0.5),
                  align         ="h",
                  common.legend = T)


## -----------------------------------------------------------------------------

# First we define the model with mrgsolve syntax
code <- "
$PARAM CL=0.0378645, V=4.59094, KA=0.45
$CMT ABS CENT AUC
$ODE
dxdt_ABS = -KA*ABS;
dxdt_CENT = KA*ABS-(CL/V)*CENT;
dxdt_AUC  = CENT*1000/V;
$TABLE double CP  = CENT*1000/V;
$CAPTURE CP
"


# This function calculates the concentration at any given times as well as the AUC
# assuming multiple doses with a tau and starting at time 0.
paramFunctionmrgsolve <- function(thetas, df, times, tau, model, ...) {

  # Load the model again since a different R process is used using parallelization
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
  if (any(names(df) == "GENO")) {
    if (df$GENO == 1) TVFREL <- TVFREL * (1 + thetas[11])
    if (df$GENO == 3) TVFREL <- TVFREL * (1 + thetas[12])
    if (df$GENO == 4) TVFREL <- TVFREL * (1 + thetas[13])
  }

  TVFREL <- FRELCOV * TVFREL

  if (!any(names(df) == "WT") || df$WT == -99) {
    TVCL <- thetas[4]
  } else {
    TVCL <- thetas[4] * (df$WT / 75)**thetas[2]
  }

  if (any(names(df) == "GENO")) {
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
  TVD1 <- thetas[7]

  FREL <- TVFREL
  CL <- TVCL
  V <- TVV

  D1FR <- TVD1
  MAT <- TVMAT

  D1 <- MAT * (1 - D1FR)
  F1 <- FREL
  KA <- 1 / (MAT - D1)

  dose <- 80 * F1 # Calculate the actual input
  rate <- dose / D1 # Calculate the actual rate

  sampletimes <- times

  dosetimes <- seq(0, max(sampletimes), by = tau)
  amt <- rep(dose, length(dosetimes))
  rate <- rep(rate, length(dosetimes))

  # Create the design dataset
  data <-
    tibble::tibble(
      ID = 1,
      cmt = 1,
      time = c(sampletimes, dosetimes),
      evid = c(rep(0, length(sampletimes)), rep(1, length(dosetimes))),
      amt = c(rep(0, length(sampletimes)), amt),
      rate = c(rep(0, length(sampletimes)), rate),
      CL = CL,
      KA = KA,
      V = V
    )

  data <- data[order(data$time), ]

  # Simulate data using the model and dataset
  out <- mrgsim_q(x = model, data = data)
  tmp <- data.frame(CONC = out$CP, AUC = out$AUC)
  # Only return the sampletimes, not the dosetimes
  return(tmp[match(sampletimes, out$time), ])
}

## -----------------------------------------------------------------------------
# Function to get the cmax within an interval [t1,t2]

Cmax <- function(thetas, df, interval, tau, model) {
  source("./paramFunctionmrgsolve.R")
  wrapper <- function(time, thetas, df, tau, model) {
    return(paramFunctionmrgsolve(thetas, df, times = c(0, time), tau = tau, model = model)[2, "CONC"])
  }
  return(optimize(
    f = wrapper, interval = interval,
    maximum = TRUE, tol = 1E-03, thetas = thetas, df = df, tau, model
  )$objective)
}

# Function to get the auc within an interval [t1,t2]
AUC <- function(thetas, df, interval, tau, model) {
  source("./paramFunctionmrgsolve.R")
  # Add a zero to ensure dosing starts at time 0
  sol <- paramFunctionmrgsolve(thetas, df, times = c(0, interval), tau = tau, model = model)
  return(sol[nrow(sol), "AUC"] - sol[nrow(sol) - 1, "AUC"])
}

## -----------------------------------------------------------------------------
mymodel <- mcode("mymodel", code, atol = 1e-8, rtol = 1e-8, maxsteps = 5000)

## ----TypicalPlotmrgsolve,fig.width=12,fig.asp=0.6,out.width="100%"------------

# Call the concentration time-profile function using a Ref covariate setting
out <- as.data.frame(paramFunctionmrgsolve(thetas = as.numeric(dfExt[2:17]), df = dfRefRow, times = seq(0, 72, 0.1), tau = 24, model = mymodel))


p <- ggplot() + geom_line(aes(x = seq(0, 72, 0.1), y = out$CONC), size = 1.8) + xlab("Time (h)") + ylab("Concentration ng/mL")
p

## -----------------------------------------------------------------------------
cmax1tau <- Cmax(thetas = as.numeric(dfExt[2:17]), df = dfRefRow, interval = seq(0, 24), tau = 24, model = mymodel)
cmax3tau <- Cmax(thetas = as.numeric(dfExt[2:17]), df = dfRefRow, interval = seq(0, 72), tau = 24, model = mymodel)

## -----------------------------------------------------------------------------

auc1dose <- AUC(thetas = as.numeric(dfExt[2:17]), df = dfRefRow, interval = c(0, 72), tau = 72, model = mymodel)
auc3dose <- AUC(thetas = as.numeric(dfExt[2:17]), df = dfRefRow, interval = c(0, 72), tau = 24, model = mymodel)

## ----Cmaxmrgsolve, cache=TRUE-------------------------------------------------
# Cmax
dfresmrgsolveCmax <- getForestDFSCM(
  dfCovs = dfCovs,
  cdfCovsNames = covnames,
  functionList = list(Cmax),
  functionListName = c("CMAX"),
  noBaseThetas = 16,
  dfParameters = dfSamplesBS,
  dfRefRow = dfRefRow,
  cstrPackages = c("mrgsolve", "dplyr"),
  ncores = 1,
  interval = c(0, 72),
  tau = 24,
  model = mymodel
)

## ----AUCmrgsolve, cache=TRUE--------------------------------------------------
# AUC
dfresmrgsolveAuc <- getForestDFSCM(
  dfCovs = dfCovs,
  cdfCovsNames = covnames,
  functionList = list(AUC),
  functionListName = c("AUC"),
  noBaseThetas = 16,
  dfParameters = dfSamplesBS,
  dfRefRow = dfRefRow,
  cstrPackages = c("mrgsolve", "dplyr"),
  ncores = 1,
  interval = c(0, 72),
  tau = 24,
  model = mymodel
)

## ----ForestCOVMRGSOLVE,fig.width=12,fig.height=30-----------------------------
plotList<-forestPlot(rbind(dfresmrgsolveCmax, dfresmrgsolveAuc),plotRelative = TRUE,keepYlabs = T,rightStrip = F,return="plotList",
                     parameterLabels = rep(c("CMAX","AUC"),each=nrow(dfresmrgsolveCmax)))
#Add the log10 x-axis
plotList[[1]]<-plotList[[1]]+scale_x_continuous(trans='log10')
plotList[[3]]<-plotList[[3]]+scale_x_continuous(trans='log10')

ggpubr::ggarrange(plotlist      = plotList,
                  ncol          = 2,
                  nrow          = 2,
                  widths        = c(1,0.5,
                                    1,0.5),
                  align         ="h",
                  common.legend = T)

