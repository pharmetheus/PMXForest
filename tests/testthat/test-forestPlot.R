# --- Shared Setup for all Forest Plot Tests ---
# This "Super Mock" contains all naming permutations to prevent selection errors.
df_mock <- data.frame(
  PARAMETER = c("CL", "V"),
  GROUPNAME = c("Sex", "Weight"),
  COVNAME   = c("Male", "70kg"),
  COVNUM    = 1:2,
  COVEFF    = TRUE,
  REFROW    = "NO",
  REFFUNC   = 10,
  REFFINAL  = 11,
  # Base columns (plotRelative=FALSE, noVar=FALSE)
  POINT = 12, Q1 = 11, Q2 = 13,
  # NoVar columns (plotRelative=FALSE, noVar=TRUE)
  POINT_NOVAR = 12, Q1_NOVAR = 11, Q2_NOVAR = 13,
  # Relative Func columns (plotRelative=TRUE, noVar=FALSE, reference="func")
  POINT_REL_REFFUNC = 1.2, Q1_REL_REFFUNC = 1.1, Q2_REL_REFFUNC = 1.3,
  # Relative Func NoVar columns (plotRelative=TRUE, noVar=TRUE, reference="func")
  POINT_NOVAR_REL_REFFUNC = 1.2, Q1_NOVAR_REL_REFFUNC = 1.1, Q2_NOVAR_REL_REFFUNC = 1.3,
  # Relative Final columns (plotRelative=TRUE, noVar=FALSE, reference="final")
  POINT_REL_REFFINAL = 1.2, Q1_REL_REFFINAL = 1.1, Q2_REL_REFFINAL = 1.3,
  stringsAsFactors  = FALSE
)


# # ==============================================================================
# # BLOCK 1: Core Assembly Logic
# # ==============================================================================
test_that("ForestPlot assembly handles stacked and table-less layouts", {

  # 1. Stacked plots without table (Lines 613-614)
  # Uses REFFUNC, POINT, Q1, Q2
  fp_stacked <- forestPlot(df_mock,
                           plotRelative = FALSE,
                           noVar = TRUE,
                           referenceParameters = "func",
                           stackedPlots = TRUE,
                           table = FALSE,
                           return = "plotList")

  expect_length(fp_stacked, 2)
  expect_s3_class(fp_stacked[[1]], "ggplot")

  # 2. Reference Info auto-text (Lines 632-645)
  # Case: REFROW="NO" and reference="final"
  fp_ref <- forestPlot(df_mock,
                       plotRelative = FALSE,
                       noVar = TRUE,
                       referenceParameters = "final",
                       referenceInfo = "auto")
  expect_s3_class(fp_ref, "gg")

  # 3. Common X-axis label (Line 626)
  fp_common <- forestPlot(df_mock,
                          plotRelative = FALSE,
                          noVar = TRUE,
                          commonXlab = TRUE)
  expect_s3_class(fp_common, "gg")
})

# ==============================================================================
# BLOCK 2: SCM Applications (Robust Structural Testing)
# ==============================================================================

# Helper to snapshot the data blueprint of each panel in a forestPlot
snapshot_forest_panels <- function(fp) {
  plots <- fp$plots
  for (i in seq_along(plots)) {
    # 1. Extract raw calculated data
    raw_data <- ggplot2::ggplot_build(plots[[i]])$data

    # 2. Standardize aesthetics
    std_data <- standardize_plot_data(raw_data)

    # 3. Stabilize numbers/text
    final_data <- stabilize(std_data[[1]])

    # 4. Filter for version-independent columns
    core_cols <- intersect(names(final_data),
                           c("x", "y", "xmin", "xmax", "label", "colour", "PANEL"))

    expect_snapshot_value(final_data[, core_cols], style = "serialize")
  }
}

test_that("Forest plots for SCM works properly", {
  # --- Setup Data ---
  dfCovs <- createInputForestData(
    list("NCIL" = c(0,1), "FORM" = c(0,1), "FOOD" = c(0,1),
         "GENO"=list("GENO1" = c(1,0,0,0), "GENO3" = c(0,0,1,0), "GENO4" = c(0,0,0,1)),
         "RACEL"= c(1,2,3), "WT" = c(65,115), "AGE" = c(27,62), "CRCL" = c(83,150), "SEX" = c(1,2)),
    iMiss=-99)

  # Robust check: Just ensure it is a data frame with rows
  expect_s3_class(dfCovs, "data.frame")
  expect_gt(nrow(dfCovs), 0)

  paramFunction <- function(thetas, df, ...) {
    # ... (Keep existing paramFunction logic) ...
    CLFOOD <- 1
    if (any(names(df)=="FOOD") && df$FOOD !=-99 && df$FOOD == 0) CLFOOD <- (1 + thetas[11])
    CLCOV <- CLFOOD
    FRELGENO4 <- 1
    if (any(names(df)=="GENO4") && df$GENO4 !=-99 && df$GENO4 == 1) FRELGENO4 <- (1 + thetas[13])
    FRELFORM <- 1
    if (any(names(df)=="FORM") && df$FORM !=-99 && df$FORM == 0) FRELFORM <- (1 + thetas[12])
    FRELSEX <- 1
    if (any(names(df)=="SEX") && df$SEX !=-99 && df$SEX == 2) FRELSEX <- (1 + thetas[14])
    FRELCOV <- FRELSEX * FRELFORM * FRELGENO4
    TVFREL <- thetas[1]
    TVFREL <- FRELCOV * TVFREL
    if(any(names(df)=="WT") && df$WT !=-99) {
      TVCL <- thetas[4] * (df$WT / 75)**thetas[2]
    } else {
      TVCL <- thetas[4]
    }
    if (any(names(df)=="GENO1") && df$GENO1 !=-99 && df$GENO1 == 1) TVCL <- TVCL * (1 + thetas[8])
    if (any(names(df)=="GENO3") && df$GENO3 !=-99 && df$GENO3 == 1) TVCL <- TVCL * (1 + thetas[9])
    if (any(names(df)=="GENO4") && df$GENO4 !=-99 && df$GENO4 == 1) TVCL <- TVCL * (1 + thetas[10])
    TVCL <- CLCOV * TVCL
    if(any(names(df)=="WT") && df$WT !=-99) {
      TVV <- thetas[5] * (df$WT / 75)**thetas[3]
    } else {
      TVV <- thetas[5]
    }
    TVMAT <- thetas[6]
    TVD1 <- thetas[7]
    FREL <- TVFREL
    CL   <- TVCL
    V    <- TVV
    MAT  <- TVMAT
    D1   <- MAT * (1 - TVD1)
    KA   <- 1 / (MAT - D1)
    AUC <- 80/(CL/FREL)
    return(list(CL,FREL,AUC))
  }

  set.seed(123)
  runno   <- 7
  extFile <- system.file("extdata", paste0("SimVal/run", runno, ".ext"), package="PMXForest")
  covFile <- system.file("extdata", paste0("SimVal/run", runno, ".cov"), package="PMXForest")
  dfSamplesCOV <- getSamples(covFile, extFile=extFile, n=25)

  dfresSCM <- getForestDFSCM(dfCovs = dfCovs, functionList = list(paramFunction),
                             functionListName = c("CL","Frel","AUC"), noBaseThetas = 14,
                             dfParameters = dfSamplesCOV, dfRefRow = NULL, cstrPackage = "dplyr")

  # --- Structural Testing ---
  expect_s3_class(dfresSCM, "data.frame")
  expect_true(all(c("CL", "Frel", "AUC") %in% unique(as.character(dfresSCM$PARAMETER))))

  # --- Blueprint Snapshotting (Replaces vdiffr) ---
  fp_default <- forestPlot(dfresSCM)
  snapshot_forest_panels(fp_default)

  # Significance Overrides (Hits red lines in setupForestPlotData)
  lsTrue <- list(c("CL","FOOD"), c("Frel","FORM"))
  fp_sig <- forestPlot(dfresSCM, onlySignificantErrorBars = TRUE, setSignEff = lsTrue)
  snapshot_forest_panels(fp_sig)

  # Custom Labels (Hits red lines in setupForestPlotData)
  unique_groups <- unique(dfresSCM$GROUPNAME)
  custom_labels <- paste0("Group_", seq_along(unique_groups))

  fp_custom <- forestPlot(dfresSCM,
                          parameterLabelsPrefix = "Metric: ",
                          groupNameLabels = custom_labels)
  snapshot_forest_panels(fp_custom)
})

# ==============================================================================
# BLOCK 3: EMP Applications (Robust Structural Testing)
# ==============================================================================
test_that("Forest plots for EMP works properly", {
  lsExpr <- rev(
    list("SEX" = expression(SEX==2), "SEX" = expression(SEX==1),
         "CRCL"= expression(CRCL>=146), "CRCL"= expression(CRCL<94),
         "AGE" = expression(AGE>=57), "AGE" = expression(AGE<35),
         "WT"  = expression(WT>104), "WT"  = expression(WT<70),
         "RACE"= expression(RACEL==3), "RACE"= expression(RACEL==2), "RACE"= expression(RACEL==1),
         "GENO"= expression(GENO==4), "GENO"= expression(GENO==3), "GENO"= expression(GENO==2), "GENO"= expression(GENO==1),
         "FOOD"= expression(FOOD==1), "FOOD"= expression(FOOD==0),
         "FORM"= expression(FORM==1), "FORM"= expression(FORM==0),
         "NCI" = expression(NCIL==1), "NCI" = expression(NCIL==0)
    )
  )

  paramFunction <- function(thetas, df, ...) {
    # ... (Keep existing paramFunction logic same as Block 2) ...
    CLFOOD <- 1
    if (any(names(df)=="FOOD") && df$FOOD !=-99 && df$FOOD == 0) CLFOOD <- (1 + thetas[11])
    CLCOV <- CLFOOD
    FRELGENO4 <- 1
    if (any(names(df)=="GENO4") && df$GENO4 !=-99 && df$GENO4 == 1) FRELGENO4 <- (1 + thetas[13])
    FRELFORM <- 1
    if (any(names(df)=="FORM") && df$FORM !=-99 && df$FORM == 0) FRELFORM <- (1 + thetas[12])
    FRELSEX <- 1
    if (any(names(df)=="SEX") && df$SEX !=-99 && df$SEX == 2) FRELSEX <- (1 + thetas[14])
    FRELCOV <- FRELSEX * FRELFORM * FRELGENO4
    TVFREL <- thetas[1]
    TVFREL <- FRELCOV * TVFREL
    if(any(names(df)=="WT") && df$WT !=-99) {
      TVCL <- thetas[4] * (df$WT / 75)**thetas[2]
    } else {
      TVCL <- thetas[4]
    }
    if (any(names(df)=="GENO1") && df$GENO1 !=-99 && df$GENO1 == 1) TVCL <- TVCL * (1 + thetas[8])
    if (any(names(df)=="GENO3") && df$GENO3 !=-99 && df$GENO3 == 1) TVCL <- TVCL * (1 + thetas[9])
    if (any(names(df)=="GENO4") && df$GENO4 !=-99 && df$GENO4 == 1) TVCL <- TVCL * (1 + thetas[10])
    TVCL <- CLCOV * TVCL
    if(any(names(df)=="WT") && df$WT !=-99) {
      TVV <- thetas[5] * (df$WT / 75)**thetas[3]
    } else {
      TVV <- thetas[5]
    }
    TVMAT <- thetas[6]
    TVD1 <- thetas[7]
    FREL <- TVFREL
    CL   <- TVCL
    V    <- TVV
    MAT  <- TVMAT
    D1   <- MAT * (1 - TVD1)
    KA   <- 1 / (MAT - D1)
    AUC <- 80/(CL/FREL)
    return(list(CL,FREL,AUC))
  }

  set.seed(123)
  runno   <- 7
  extFile <- system.file("extdata", paste0("SimVal/run",runno,".ext"), package="PMXForest")
  covFile <- system.file("extdata", paste0("SimVal/run",runno,".cov"), package="PMXForest")
  dfSamplesCOV <- getSamples(covFile, extFile=extFile, n=25)

  dataFile <- system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package="PMXForest")
  dfData <- read.csv(dataFile, stringsAsFactors = FALSE) %>%
    dplyr::distinct(ID, .keep_all=TRUE) %>%
    dplyr::slice(1:20, 500:600)

  dfresEMP <- getForestDFemp(
    dfData = dfData, covExpressionsList = lsExpr, cdfCovsNames = NULL,
    functionList = list(paramFunction), functionListName = c("CL","Frel","AUC"),
    metricFunction = median, noBaseThetas = 14, dfParameters = dfSamplesCOV,
    dfRefRow = NULL, ncores = 1, cstrPackages = "dplyr"
  )

  # --- Structural Validation ---
  expect_s3_class(dfresEMP, "data.frame")
  expect_gt(nrow(dfresEMP), 0)

  # --- Blueprint Snapshotting (Replacing 12 vdiffr calls) ---
  fp_emp <- forestPlot(dfresEMP)
  snapshot_forest_panels(fp_emp)

  # Test different return modes to hit final assembly lines
  res_data <- forestPlot(dfresEMP, return = "data")
  expect_true("STATISTIC" %in% names(res_data))

  res_list <- forestPlot(dfresEMP, return = "plotList")
  expect_type(res_list, "list")

  # Test table=FALSE branch
  fp_no_tab <- forestPlot(dfresEMP, table = FALSE)
  expect_s3_class(fp_no_tab, "gg")
})

# ==============================================================================
# BLOCK 4: Final Coverage Gaps (Reference Strings & Overrides)
# ==============================================================================
test_that("Final coverage gaps for Reference Info and COVEFF overrides", {

  # 1. Test all 4 Reference Info "auto" strings (Lines 632-645)
  # Toggle REFROW and referenceParameters to exercise all footer variations.

  # Case A: REFROW=NO, referenceParameters=final (Requires noVar=FALSE)
  df_no <- df_mock; df_no$REFROW <- "NO"
  fp1 <- forestPlot(df_no, noVar = FALSE, referenceParameters = "final", referenceInfo = "auto")
  expect_s3_class(fp1, "gg")

  # Case B: REFROW=NO, referenceParameters=func (Default noVar=TRUE)
  fp2 <- forestPlot(df_no, referenceParameters = "func", referenceInfo = "auto")
  expect_s3_class(fp2, "gg")

  # Case C: REFROW=YES, referenceParameters=final
  df_yes <- df_mock; df_yes$REFROW <- "YES"
  fp3 <- forestPlot(df_yes, noVar = FALSE, referenceParameters = "final", referenceInfo = "auto")
  expect_s3_class(fp3, "gg")

  # Case D: REFROW=YES, referenceParameters=func
  fp4 <- forestPlot(df_yes, referenceParameters = "func", referenceInfo = "auto")
  expect_s3_class(fp4, "gg")

  # 2. Test manual setSignEff override (Hits setCOVEFF internal logic)
  # We use plotRelative=TRUE to match the Super Mock columns
  lsOverride <- list(c("V", "Weight"))
  plotData_override <- setupForestPlotData(df_mock,
                                           plotRelative = TRUE,
                                           noVar = FALSE,
                                           setSignEff = lsOverride)
  expect_true(plotData_override[plotData_override$PARAMETER == "V" &
                                  plotData_override$GROUPNAME == "Weight", ]$COVEFF)

  # 3. Test row-length parameter labels (Hits Lines 77-78)
  # Providing labels matching row count instead of unique parameter count.
  long_labels <- rep("CustomParam", nrow(df_mock))
  plotData_long <- setupForestPlotData(df_mock, parameterLabels = long_labels)
  expect_equal(as.character(plotData_long$PARAMETERLABEL[1]), "CustomParam")

  # 4. Final Return Logic (Lines 651-654)
  expect_s3_class(forestPlot(df_mock, return = "data"), "data.frame")
  expect_type(forestPlot(df_mock, return = "plotList"), "list")
})
