test_that("mvrnorm_vector samples correctly", {

  expect_equal(mvrnorm_vector(1,sigma=1,iSampleIndex=0),1)
  expect_equal(mvrnorm_vector(c(1,2),sigma=matrix(c(0,0,0,1),ncol=2),fixed_mu = c(0,1)),c(1,2))

  ## Read the ext-file
  extFile <- getExt(system.file("extdata","SimVal/run7.ext",package="PMXForest"))

  ## Read the cov-file and extract the covmatrix
  dfcov     <- read.table(system.file("extdata","SimVal/run7.cov",package="PMXForest"),
                          fill = TRUE, header = TRUE, sep = "",
                          skip = 1, stringsAsFactors = FALSE)
  sigma     <- data.matrix(dfcov[, 2:ncol(dfcov)])

  ## Get the final parameter estimates
  finPar  <- subset(extFile, ITERATION == "-1000000000")
  mu      <- as.numeric(finPar[, -(c(1, ncol(finPar)))])

  #Get the parameters which are fixed based on the covariance
  fixedmu <-  rep(FALSE,1,ncol(sigma))
  for (j in 1:ncol(sigma)) {
    fixedmu[j]<-all(sigma[,j]==0)
  }

  ## Test iSampleIndex
  dfParameters <- as.data.frame(mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = 10,fixed_mu=fixedmu))
  expect_equal(nrow(dfParameters),10)

  dfParameters <- mvrnorm_vector(mu = mu, sigma = sigma, iSampleIndex = 1,fixed_mu = fixedmu)
  expect_type(dfParameters,"double")

})
