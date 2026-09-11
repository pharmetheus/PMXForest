test_that("setupDfRefRow generates correct singleRef = TRUE geometry", {
  mock_data <- data.frame(
    ID = 1:5,
    WT = c(60, 70, 70, 80, 90), # Median 70
    SEX = c(1, 1, 1, 0, 0) # Mode 1
  )

  df_covs <- setupDfCovs(mock_data, covariates = c("WT", "SEX"), minLevels = 3)

  ref_single <- setupDfRefRow(
    dfCovs = df_covs,
    data = mock_data,
    covariates = c("WT", "SEX"),
    minLevels = 3
  )

  expect_equal(nrow(ref_single), 1, info = "singleRef = TRUE did not return a single row.")
  expect_equal(ref_single$WT, 70, info = "Continuous reference incorrect.")
  expect_equal(ref_single$SEX, 1, info = "Categorical reference incorrect.")
  expect_equal(ref_single$COVARIATEGROUPS, "Reference")
})

test_that("setupDfRefRow generates correct singleRef = FALSE geometry", {
  mock_data <- data.frame(
    ID = 1:5,
    WT = c(60, 70, 70, 80, 90), # Median 70
    SEX = c(1, 1, 1, 0, 0) # Mode 1
  )

  # WT and SEX have -99s in inactive cells
  df_covs <- setupDfCovs(mock_data, covariates = c("WT", "SEX"), minLevels = 3)

  ref_matrix <- setupDfRefRow(
    dfCovs = df_covs,
    data = mock_data,
    covariates = c("WT", "SEX"),
    singleRef = FALSE,
    minLevels = 3
  )

  expect_equal(nrow(ref_matrix), nrow(df_covs), info = "singleRef = FALSE did not match dfCovs rows.")

  # Extract the WT varying block.
  # WT should be 70 everywhere in this block. SEX should remain -99.
  wt_block <- ref_matrix[ref_matrix$COVARIATEGROUPS == "WT", ]
  expect_true(all(wt_block$WT == 70), info = "Active continuous cells were not overwritten with reference.")
  expect_true(all(wt_block$SEX == -99), info = "Inactive missVal cells were not preserved.")
})

test_that("setupDfRefRow honours refLevels so its columns match setupDfCovs", {
  mock_data <- data.frame(
    ID   = 1:8,
    WT   = c(60, 70, 70, 80, 90, 65, 75, 72),
    GENO = c(1, 2, 2, 3, 4, 2, 3, 2) # mode is 2
  )

  df_covs <- setupDfCovs(
    mock_data,
    covariates = c("WT", "GENO"), catRef = list(GENO = 2)
  )

  # refLevels is deprecated; it must still work and forward to catRef.
  expect_warning(
    ref <- setupDfRefRow(
      dfCovs     = df_covs,
      data       = mock_data,
      covariates = c("WT", "GENO"),
      refLevels  = list(GENO = 2)
    ),
    "`refLevels` is deprecated"
  )

  # Column names must line up with df_covs (GENO_1, GENO_3, GENO_4)
  expect_true(all(c("GENO_1", "GENO_3", "GENO_4") %in% names(ref)))
  # Mode genotype is 2 (the reference) -> all GENO dummies at 0
  expect_equal(ref$GENO_1, 0)
  expect_equal(ref$GENO_3, 0)
  expect_equal(ref$GENO_4, 0)
})

test_that("setupDfRefRow errors when a covariate is entirely missing", {
  mock_data <- data.frame(
    ID = 1:4,
    WT = c(60, 70, 80, 90),
    GONE = c(-99, -99, -99, -99)
  )
  df_covs <- setupDfCovs(mock_data, covariates = "WT", idVar = "ID")

  expect_error(
    setupDfRefRow(df_covs,
      data = mock_data, covariates = c("WT", "GONE"),
      idVar = "ID"
    ),
    "contains only missing values"
  )
})

test_that("catRef and sep are inherited from the dfCovs they were used to build", {
  # Requiring the same two arguments in two calls, kept in step by hand, is a
  # standing invitation to a silent mismatch: a different sep gives GENO_1 where
  # dfCovs has GENO1, and the reference for that covariate is quietly lost.
  d <- read.csv(system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv",
    package = "PMXForest"
  ))
  dfCovs <- setupDfCovs(d,
    covariates = c("WT", "SEX", "GENO"),
    catRef = list(GENO = 2), sep = "", idVar = "ID"
  )

  inherited <- setupDfRefRow(dfCovs, d,
    covariates = c("WT", "SEX", "GENO"), idVar = "ID"
  )
  explicit <- setupDfRefRow(dfCovs, d,
    covariates = c("WT", "SEX", "GENO"),
    catRef = list(GENO = 2), sep = "", idVar = "ID"
  )
  expect_identical(inherited, explicit)
  expect_true(all(c("GENO1", "GENO3", "GENO4") %in% names(inherited)))

  # an argument given explicitly still wins over the recorded one
  override <- suppressMessages(setupDfRefRow(dfCovs, d,
    covariates = c("WT", "SEX", "GENO"), sep = "_", idVar = "ID"
  ))
  expect_true(all(override$GENO1 == -99)) # looked for GENO_1, found nothing

  # and a dfCovs without the attribute behaves exactly as before
  plain <- dfCovs
  attr(plain, "pmxCovSetup") <- NULL
  expect_identical(
    suppressMessages(setupDfRefRow(plain, d,
      covariates = c("WT", "SEX", "GENO"), idVar = "ID"
    )),
    suppressMessages(setupDfRefRow(dfCovs, d,
      covariates = c("WT", "SEX", "GENO"), sep = "_", idVar = "ID"
    ))
  )
})

test_that("a covariate in dfCovs but not in covariates falls back to missVal", {
  # covariates need not cover every column of dfCovs. An unresolved column gets
  # missVal, which every parameter function already understands as "not active
  # on this row" and answers with its own reference. It must not be NA: the
  # generated preamble tests `df[[cov]] != missVal`, and `if (NA)` is an error.
  d <- read.csv(system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv",
    package = "PMXForest"
  ))
  dfCovs <- setupDfCovs(d,
    covariates = c("WT", "SEX", "GENO"),
    catRef = list(GENO = 2), sep = "", idVar = "ID"
  )

  expect_message(
    partial <- setupDfRefRow(dfCovs, d, covariates = c("WT", "SEX"), idVar = "ID"),
    "set to -99"
  )
  expect_equal(unname(unlist(partial[, c("GENO1", "GENO3", "GENO4")])), rep(-99, 3))
  expect_false(anyNA(partial))

  # the row is usable: a generated parameter function runs on it
  gen <- createParamFunction(
    system.file("extdata", "SimVal/run7.mod", package = "PMXForest"),
    parameters = c("CL", "V"),
    covRef = list(GENO1 = 0, GENO3 = 0), quiet = TRUE
  )
  fun <- eval(parse(text = gen$code))
  ext <- getExt(system.file("extdata", "SimVal/run7.ext", package = "PMXForest"))
  thetas <- as.numeric(ext[ext$ITERATION == -1000000000, 2:15])
  vals <- fun(thetas, partial[, setdiff(names(partial), "COVARIATEGROUPS"), drop = FALSE])
  expect_false(anyNA(unlist(vals)))

  # singleRef = FALSE: a cell is replaced by the reference only where the
  # covariate is active on that row. With no reference resolved for GENO, its
  # own rows keep the dfCovs values and every other row stays at missVal - so
  # the GENO rows end up equal to the reference and plot at 1, which is the
  # honest answer for a covariate no reference was asked for.
  wide <- suppressMessages(setupDfRefRow(dfCovs, d,
    covariates = c("WT", "SEX"),
    idVar = "ID", singleRef = FALSE
  ))
  isGeno <- wide$COVARIATEGROUPS == "GENO"
  expect_true(all(wide$GENO1[!isGeno] == -99))
  expect_equal(wide$GENO1[isGeno], dfCovs$GENO1[dfCovs$COVARIATEGROUPS == "GENO"])
  expect_false(anyNA(wide))
})
