test_that("getCovStats works", {

  data <- read.csv(system.file("extdata/SimVal/DAT-1-MI-PMX-2.csv",package = "PMXForest"))
  covariates <- c("WT","HT","SEX")

  expect_snapshot(getCovStats(data,covariates,missVal=-99))
  expect_snapshot(getCovStats(data,covariates,missVal=-99,probs=c(0.1,0.9)))

})
