## ---------------------------------------------------------------------------
## Compare the hand-written run7.mod parameter function with the one generated
## by PMXForest::createParamFunction().
##
## Both Forest plots are built from the same covariate structure, the same
## reference row and the SAME parameter samples, so any visible difference is
## attributable to the parameter function alone. The script reports the largest
## absolute difference in the Forest plot data and writes the two plots side by
## side; compare-paramFunction.png in this directory is the reference output.
##
## Run from any working directory:
##
##   Rscript $(Rscript -e 'cat(system.file("validation/compare-paramFunction.R",
##                                         package = "PMXForest"))')
##
## or, from a source checkout:  Rscript inst/validation/compare-paramFunction.R
##
## The plot is written to the current working directory.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(PMXForest)
  library(ggplot2)
})

nSamples <- 100
seed     <- 865765
outFile  <- "compare-paramFunction.png"

modFile  <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")
extFile  <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")
covFile  <- system.file("extdata", "SimVal/run7.cov", package = "PMXForest")
dataFile <- system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")

## ---------------------------------------------------------------------------
## 1. The hand-written parameter function
##
## This is now the only place the hand-written run7 function survives - the
## Walkthrough it was transcribed from generates its function with
## createParamFunction(). That is the right home for it: the script exists to
## compare a hand-written function against a generated one, so keeping the
## hand-written half here rather than in a tutorial is the point.
##
## Note how the -99 guard is repeated at every covariate reference, and how the
## same covariate is tested more than once.
## ---------------------------------------------------------------------------

manualFunction <- function(thetas, df, ...) {

  CLFOOD <- 1
  if (any(names(df) == "FOOD") && df$FOOD != -99 && df$FOOD == 0) CLFOOD <- 1 + thetas[11]

  FRELFORM <- 1
  if (any(names(df) == "FORM") && df$FORM != -99 && df$FORM == 0) FRELFORM <- 1 + thetas[12]

  FRELSEX <- 1
  if (any(names(df) == "SEX") && df$SEX != -99 && df$SEX == 2) FRELSEX <- 1 + thetas[14]

  FRELGENO4 <- 1
  if (any(names(df) == "GENO4") && df$GENO4 != -99 && df$GENO4 == 1) FRELGENO4 <- 1 + thetas[13]

  FREL <- thetas[1] * FRELSEX * FRELFORM * FRELGENO4

  if (any(names(df) == "WT") && df$WT != -99) {
    TVCL <- thetas[4] * (df$WT / 75)^thetas[2]
  } else {
    TVCL <- thetas[4]
  }
  if (any(names(df) == "GENO1") && df$GENO1 != -99 && df$GENO1 == 1) TVCL <- TVCL * (1 + thetas[8])
  if (any(names(df) == "GENO3") && df$GENO3 != -99 && df$GENO3 == 1) TVCL <- TVCL * (1 + thetas[9])
  if (any(names(df) == "GENO4") && df$GENO4 != -99 && df$GENO4 == 1) TVCL <- TVCL * (1 + thetas[10])
  CL <- CLFOOD * TVCL

  if (any(names(df) == "WT") && df$WT != -99) {
    V <- thetas[5] * (df$WT / 75)^thetas[3]
  } else {
    V <- thetas[5]
  }

  list(CL = CL, FREL = FREL, V = V)
}

## ---------------------------------------------------------------------------
## 2. The generated parameter function
##
## GENO1 and GENO3 are 0/1 dummies whose reference level the control stream does
## not state outright, so createParamFunction() proposes 0 and warns. That is
## correct here - with every dummy at 0 the genotype is the reference category -
## so confirm it explicitly through covRef and the warning goes away.
## ---------------------------------------------------------------------------

generated <- createParamFunction(
  modFile,
  parameters = c("CL", "FREL", "V"),
  covRef     = list(GENO1 = 0, GENO3 = 0),
  quiet      = FALSE
)

cat("\n---------------- generated source ----------------\n")
print(generated$code)
cat("--------------------------------------------------\n\n")

generatedFunction <- eval(parse(text = generated$code))

## ---------------------------------------------------------------------------
## 3. Shared inputs
## ---------------------------------------------------------------------------

dfData <- read.csv(dataFile)
covs   <- c("WT", "SEX", "FOOD", "FORM", "GENO1", "GENO3", "GENO4")

dfCovs <- setupDfCovs(dfData, covariates = covs, idVar = "ID")

## Take the reference from the control stream rather than from the data.
##
## Both parameter functions fall back to the model's own reference when a
## covariate is inactive - WT at the 75 kg normalisation weight, FORM at the
## level marked "; Most common". A reference row built from data statistics
## would instead sit at the median weight (85.4 kg) and the modal formulation,
## and every row where those covariates are inactive would be displaced from 1:
## (75/85.4)^theta2 = 0.907 for CL and (75/85.4)^theta3 = 0.878 for V. Passing
## the generated object as `model` guarantees the reference row and the
## parameter function come from one derivation.
dfRefRow <- setupDfRefRow(dfCovs, dfData, covariates = covs, singleRef = TRUE,
                          idVar = "ID", contRef = "model", catRef = "model",
                          model = generated)

## The same samples for both plots, so the comparison isolates the function.
set.seed(seed)
dfSamples <- getSamples(covFile, extFile, n = nSamples)

forestData <- function(fun) {
  getForestDFSCM(
    dfCovs           = dfCovs,
    functionList     = list(fun),
    functionListName = c("CL", "Frel", "V"),
    noBaseThetas     = generated$noBaseThetas,
    dfParameters     = dfSamples,
    dfRefRow         = dfRefRow,
    ncores           = 1
  )
}

dfresManual    <- forestData(manualFunction)
dfresGenerated <- forestData(generatedFunction)

## ---------------------------------------------------------------------------
## 4. Numeric comparison - the plots are the illustration, this is the evidence
## ---------------------------------------------------------------------------

numCols <- c("POINT", "Q1", "Q2", "POINT_REL_REFFUNC",
             "Q1_REL_REFFUNC", "Q2_REL_REFFUNC", "REFFUNC")

stopifnot(nrow(dfresManual) == nrow(dfresGenerated),
          identical(as.character(dfresManual$COVNAME),
                    as.character(dfresGenerated$COVNAME)),
          identical(as.character(dfresManual$PARAMETER),
                    as.character(dfresGenerated$PARAMETER)))

diffs <- vapply(numCols, function(cl)
  max(abs(dfresManual[[cl]] - dfresGenerated[[cl]])), numeric(1))

cat("Rows compared: ", nrow(dfresManual), "\n", sep = "")
cat("Largest absolute difference, by column:\n")
print(data.frame(column = names(diffs), maxAbsDiff = as.numeric(diffs),
                 row.names = NULL))

if (max(diffs) == 0) {
  cat("\nThe two parameter functions give bit-for-bit identical Forest plot data.\n\n")
} else {
  cat(sprintf("\nLargest difference anywhere: %.3e\n\n", max(diffs)))
}

## ---------------------------------------------------------------------------
## 5. The two Forest plots, side by side
## ---------------------------------------------------------------------------

## forestPlot() already returns a ggpubr composite, so the title goes on with
## annotate_figure() rather than ggtitle().
plotOne <- function(dfres, title) {
  ggpubr::annotate_figure(
    forestPlot(dfres, parameters = c("CL", "Frel", "V"), sigdigits = 3),
    top = ggpubr::text_grob(title, face = "bold", size = 15)
  )
}

pManual    <- plotOne(dfresManual,    "Hand-written parameter function")
pGenerated <- plotOne(dfresGenerated, "Generated by createParamFunction()")

combined <- ggpubr::ggarrange(pManual, pGenerated, ncol = 2, nrow = 1)

## Each Forest plot is itself six facets wide, so the canvas has to be generous
## or the statistics column gets clipped.
ggsave(outFile, combined, width = 26, height = 9, dpi = 140, bg = "white")
cat("Side-by-side plot written to ", normalizePath(outFile), "\n", sep = "")
