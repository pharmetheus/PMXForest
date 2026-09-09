library(testthat)

make_data <- function() {
  data.frame(
    ID   = 1:6,
    GENO = c(1, 2, 3, 4, 2, -99),
    RACE = c(1, 1, 2, 3, 2, 1),
    WT   = c(60, 70, 80, 90, 65, 75)
  )
}

test_that("character spec encodes with the lowest level as reference", {
  out <- oneHotEncode(make_data(), spec = "GENO")

  expect_true(all(c("GENO_2", "GENO_3", "GENO_4") %in% names(out)))
  expect_false("GENO_1" %in% names(out))

  # Row values: level 1 -> all zero; level k -> GENO_k == 1
  expect_equal(out$GENO_2, c(0, 1, 0, 0, 1, -99))
  expect_equal(out$GENO_3, c(0, 0, 1, 0, 0, -99))
  expect_equal(out$GENO_4, c(0, 0, 0, 1, 0, -99))

  # Raw column retained by default
  expect_true("GENO" %in% names(out))
})

test_that("missing rows are kept as missVal by default and imputed when asked", {
  d <- make_data()

  keep <- oneHotEncode(d, spec = "GENO")
  expect_equal(keep$GENO_2[6], -99)

  imp <- oneHotEncode(d, spec = "GENO", imputeMissing = TRUE)
  expect_equal(imp$GENO_2[6], 0)
  expect_equal(imp$GENO_3[6], 0)
  expect_equal(imp$GENO_4[6], 0)
})

test_that("custom missVal is honoured", {
  d <- make_data()
  d$GENO[6] <- -999
  out <- oneHotEncode(d, spec = "GENO", missVal = -999)
  expect_equal(out$GENO_2[6], -999)
})

test_that("explicit reference level changes which dummy columns are created", {
  out <- oneHotEncode(make_data(), spec = list(GENO = list(ref = 2)))

  expect_true(all(c("GENO_1", "GENO_3", "GENO_4") %in% names(out)))
  expect_false("GENO_2" %in% names(out))

  expect_equal(out$GENO_1, c(1, 0, 0, 0, 0, -99))
  expect_equal(out$GENO_3, c(0, 0, 1, 0, 0, -99))
})

test_that("a bare scalar spec element is treated as the reference level", {
  a <- oneHotEncode(make_data(), spec = list(GENO = 2))
  b <- oneHotEncode(make_data(), spec = list(GENO = list(ref = 2)))
  expect_equal(a, b)
})

test_that("levels sub-spec limits the dummy columns", {
  out <- oneHotEncode(make_data(), spec = list(GENO = list(ref = 1, levels = c(3, 4))))
  expect_true(all(c("GENO_3", "GENO_4") %in% names(out)))
  expect_false("GENO_2" %in% names(out))
})

test_that("includeReference adds the reference dummy", {
  out <- oneHotEncode(make_data(), spec = "GENO", includeReference = TRUE)
  expect_true(all(c("GENO_1", "GENO_2", "GENO_3", "GENO_4") %in% names(out)))
  expect_equal(out$GENO_1, c(1, 0, 0, 0, 0, -99))
})

test_that("sep controls the column name separator", {
  out <- oneHotEncode(make_data(), spec = "GENO", sep = "")
  expect_true(all(c("GENO2", "GENO3", "GENO4") %in% names(out)))
})

test_that("dropOriginal removes the raw column", {
  out <- oneHotEncode(make_data(), spec = "GENO", dropOriginal = TRUE)
  expect_false("GENO" %in% names(out))
  expect_true("GENO_2" %in% names(out))
})

test_that("encoding is idempotent when a consistent dummy column already exists", {
  d <- oneHotEncode(make_data(), spec = "GENO")
  again <- oneHotEncode(d, spec = "GENO")
  expect_equal(again, d)
})

test_that("a pre-existing FREM-style imputed column is accepted", {
  d <- make_data()
  # Imputed encoding: missing row -> 0 rather than missVal
  d$GENO_2 <- c(0, 1, 0, 0, 1, 0)
  out <- oneHotEncode(d, spec = "GENO")
  # existing GENO_2 left untouched, the other dummies added
  expect_equal(out$GENO_2, c(0, 1, 0, 0, 1, 0))
  expect_equal(out$GENO_3, c(0, 0, 1, 0, 0, -99))
})

test_that("an inconsistent pre-existing column is an error", {
  d <- make_data()
  d$GENO_2 <- c(1, 1, 1, 1, 1, 1) # wrong
  expect_error(oneHotEncode(d, spec = "GENO"), "inconsistent")
})

test_that("unknown covariate and single-level covariate are skipped with a warning", {
  d <- make_data()
  d$CONST <- 1
  expect_warning(oneHotEncode(d, spec = "NOPE"), "not found")
  expect_warning(oneHotEncode(d, spec = "CONST"), "fewer than two")
})

test_that("COVARIATEGROUPS is kept in the last position", {
  d <- make_data()
  d$COVARIATEGROUPS <- "GENO"
  out <- oneHotEncode(d, spec = "GENO")
  expect_equal(names(out)[length(names(out))], "COVARIATEGROUPS")
})

test_that("character-valued categorical levels are supported", {
  d <- data.frame(
    ID    = 1:4,
    Class = c("1st", "2nd", "3rd", "1st")
  )
  out <- oneHotEncode(d, spec = list(Class = list(ref = "1st")))
  expect_true(all(c("Class_2nd", "Class_3rd") %in% names(out)))
  expect_equal(out$Class_2nd, c(0, 1, 0, 0))
})

test_that("non-data.frame input errors", {
  expect_error(oneHotEncode(list(a = 1), spec = "a"), "data.frame")
})

test_that("a spec that is neither a character vector nor a named list is an error", {
  d <- make_data()
  expect_error(oneHotEncode(d, spec = 5), "character vector .* or a named list")
  expect_error(oneHotEncode(d, spec = list(1, 2)), "character vector .* or a named list")
})

test_that("a spec element that is a multi-value vector is an error", {
  d <- make_data()
  expect_error(
    oneHotEncode(d, spec = list(GENO = c(1, 2))),
    "must be NULL, a single reference value, or a list"
  )
})

test_that("includeReference adds the reference to an explicit levels sub-spec", {
  d <- make_data()
  out <- oneHotEncode(
    d,
    spec = list(GENO = list(ref = 1, levels = c(3, 4))),
    includeReference = TRUE
  )
  # ref (1) prepended to c(3, 4)
  expect_true(all(c("GENO_1", "GENO_3", "GENO_4") %in% names(out)))
  expect_false("GENO_2" %in% names(out))
  expect_equal(out$GENO_1, c(1, 0, 0, 0, 0, -99))
})
