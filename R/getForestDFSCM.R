#' getForestDFSCM
#'
#' @description Get a data frame with Forest border for each univariate or multivariate covariate
#' (and value(s)) in the input data frame. If a list a data frame will be created from the list,
#' see function dfCreateInputForestData
#'
#' @details
#'   **Every row is conditioned on the reference subject.** A row of `dfCovs`
#'   changes the covariates it sets and leaves the others at `-99`, which a
#'   generated parameter function reads as the reference of the model. Each row
#'   is divided by the value for `dfRefRow`. A row therefore shows how the
#'   parameter of the reference subject changes when only that covariate
#'   changes, not the parameter of a typical subject with that covariate value.
#'   To condition every row on another value, set that covariate on every row of
#'   `dfCovs` and give `dfRefRow` the same value. [getForestDFemp()] instead
#'   evaluates the subjects in the data who match each row, with all their own
#'   covariate values.
#'
#'   **Covariate names are checked against a generated parameter function.** A
#'   function made by [createParamFunction()] carries the names of the data
#'   columns it reads, as the attribute `pmxCovariates`. `getForestDFSCM()`
#'   compares them with the columns of `dfCovs`, after the encoding set by
#'   `oneHot` and without `COVARIATEGROUPS`. It warns when two differences
#'   occur together: a column that no function in `functionList` reads, so
#'   that the rows on which only that column is set have no effect, and a
#'   covariate that a function reads and `dfCovs` lacks, so that it is held at
#'   the reference of the function on every row. A name spelled differently
#'   from `$PK` always gives that pair, for example `GENO_1` in `dfCovs` for
#'   `GENO1` in the model (the separator, `sep` in [setupDfCovs()]), or a
#'   `GENO2` column for `GENO1` (the reference level, `catRef`). Either
#'   difference alone is often a deliberate choice of rows and gives no
#'   warning. A hand-written function carries no names, so no check is made
#'   when `functionList` holds one. The warning has the class
#'   `pmxCovariateMismatch`, so an intended pair can be silenced alone:
#'   `withCallingHandlers(getForestDFSCM(...), pmxCovariateMismatch = function(w) invokeRestart("muffleWarning"))`.
#'
#'   **The columns of `dfRefRow` are checked against `dfCovs`.** A column of
#'   `dfRefRow` that `dfCovs` lacks, after the encoding set by `oneHot`, gives
#'   a warning of the same class, with the names in the field `refOnly`. With
#'   a name spelled differently, such as `wt` in `dfRefRow` for `WT` in
#'   `dfCovs`, the reference subject is held at the reference of the parameter
#'   function for that covariate, and every row is shifted. This check is made
#'   for hand-written functions too. A reference row from [setupDfRefRow()]
#'   has the columns of `dfCovs` and is never reported.
#'
#' @import doParallel
#' @import foreach
#' @import dplyr
#'
#' @param dfCovs A data frame with covariates to include, if a covariate value is set -99 or NA,
#' they are assumed missing and will not be included in any FFEM transformations. If a dfCovs is
#' a list an attempt will be made to create the appropriate data frame with the createInputForestData function.
#' @param cdfCovsNames A string vector with names of the rows in dfCovs, if not used, names will be
#' automatically assigned based on the covariate values and column names in dfCovs.
#' @param functionList A list of functions with input (basethetas, covthetas,dfrow and ...) from which the change
#' from the reference value will be calculated. If the function returns a vector of values, each value will be used
#' but functionListName must contain the names with a length of all return for all functions in the functionList
#' @param functionListName A vector of strings (names) of the parameters for each function in the functionList
#' @param noBaseThetas the number of parameters from `dfParameters` to pass to the functions in `functionList`.
#' @param dfParameters A data frame with parameter samples from the uncertainty distribution.
#' The vector of final parameter estimates is assumed to be in the first row.
#' The column order is assumed the same as in the NONMEM ext file except the ITERATION and OBJ
#' columns whichshould not be included.
#' @param quiet If output should be allowed during the function call, default= TRUE. (This option
#' is mainly for debugging purposes.)
#' @param probs A vector of probabilities that should be computed for each of the parameters from
#' functionList. These will be used as the
#' as the uncertainties in the Forest plots. The probs vector position one and two will be used
#' for plotting the uncertanties (i.e. columns q1 and q2). Default is c(0.05, 0.95).
#' @param pointFunction The function used to calculate the point for each covariate in the forest
#' plot. default=median
#' This function is also used for the reference covariate combination
#' @param dfRefRow A data frame  (one row or equal number of rows as dfCovs) with the covariate
#' values that will be used as the reference, if NULL the typical subject is used as reference.
#' @param cGrouping A vector of numbers defining how to group the y-axis of the Forest plot, the
#' length of the vector should match the number of rows in dfCovs.
#' If NULL (default) an educated guess of the grouping will be set
#' @param ncores the number of cores to use for the calculations, default = 1 which means no parallellization.
#'   On Windows the workers are new R sessions. They load the same package
#'   versions as this session, including a package loaded from a folder that
#'   is not on `.libPaths()`, and the run stops if they cannot.
#' @param cstrPackages a character vector with package names needed to run the calculations in parallel, default = NULL
#' @param cstrExports a character vector with variables needed to run the calculations in parallel, default = NULL
#' @param iMiss The missing value number. -99 by default.
#' @param oneHot An optional one-hot encoding specification (see
#' \code{\link{oneHotEncode}}). When supplied, raw multi-level categorical columns
#' are one-hot encoded before the parameter functions are evaluated - in `dfCovs`
#' and `dfRefRow` for `getForestDFSCM()`, in `dfData` for `getForestDFemp()` - so
#' a single parameter function written against the dummy columns can serve both
#' the parametric and empirical workflows. Default `NULL` (no encoding). For
#' `getForestDFSCM()` the raw column is dropped after encoding.
#' @param oneHotSep The separator between the covariate name and the level in the
#' dummy column names created when `oneHot` is used. Defaults to `"_"`. Set to
#' `""` to match dummy columns named like `GENO1`, `GENO2`.
#' @param ... additional variables to be forwarded to the the functionList functions
#'
#'
#' @return A data frame with summary statistics for each parameters and covariate combinations:
#'
#' @seealso
#'   `vignette("Part2-walkthrough", package = "PMXForest")` for how this fits the whole workflow.
#'
#' @export
#'
#' @examples
#' dfData <- read.csv(
#'   system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
#' )
#' extFile <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")
#' covFile <- system.file("extdata", "SimVal/run7.cov", package = "PMXForest")
#'
#' # Covariate rows and uncertainty samples (run7 has 14 THETAs)
#' dfCovs <- setupDfCovs(dfData,
#'   covariates = c("WT", "AGE", "CRCL", "FOOD"),
#'   idVar = "ID"
#' )
#' dfSamples <- getSamples(covFile, extFile, n = 50)
#'
#' # Rebuild the part of the model needed for CL
#' paramFunction <- function(thetas, df, ...) {
#'   TVCL <- thetas[4]
#'   if (any(names(df) == "WT") && df$WT != -99) {
#'     TVCL <- thetas[4] * (df$WT / 75)^thetas[2]
#'   }
#'   if (any(names(df) == "FOOD") && df$FOOD != -99 && df$FOOD == 0) {
#'     TVCL <- TVCL * (1 + thetas[11])
#'   }
#'   list(CL = TVCL)
#' }
#'
#' dfres <- getForestDFSCM(
#'   dfCovs           = dfCovs,
#'   functionList     = list(paramFunction),
#'   functionListName = "CL",
#'   noBaseThetas     = 14,
#'   dfParameters     = dfSamples
#' )
#' head(dfres[, c("COVNAME", "GROUPNAME", "PARAMETER", "POINT", "Q1", "Q2")])
getForestDFSCM <- function(dfCovs,
                           cdfCovsNames = NULL,
                           functionList = list(
                             function(basethetas, dfrow, ...) {
                               return(basethetas[1])
                             }
                           ),
                           functionListName = "PAR1",
                           noBaseThetas,
                           dfParameters,
                           quiet = TRUE,
                           probs = c(0.05, 0.95),
                           pointFunction = median,
                           dfRefRow = NULL,
                           cGrouping = NULL,
                           ncores = 1,
                           cstrPackages = NULL,
                           cstrExports = NULL,
                           iMiss = -99,
                           oneHot = NULL,
                           oneHotSep = "_",
                           ...) {
  if (!is.null(dfRefRow) && nrow(dfRefRow) != 1 && nrow(dfRefRow) != nrow(dfCovs)) {
    stop(
      "The number of reference rows (dfRefRow) should be either NULL (missing used as reference), one ",
      "(this row used as reference) or equal to dfCovs (change reference for each covariate ",
      "combination)"
    )
  }

  ## Remove samples with problems. Will use THETA1 == NA as an indicator for a problematic sample
  dfParameters <- dfParameters %>% filter(!is.na(THETA1))

  ## Try to make dfCovs into a data.frame if it isn't that already
  if (!is.data.frame(dfCovs)) {
    dfCovs <- createInputForestData(dfCovs)
  }
  ## Coerce to a plain data.frame: the internal `[` code assumes base-R
  ## drop-to-vector behaviour, which a tibble does not provide.
  dfCovs <- as.data.frame(dfCovs)
  if (!is.null(dfRefRow) && is.data.frame(dfRefRow)) dfRefRow <- as.data.frame(dfRefRow)
  dfCovs[is.na(dfCovs)] <- iMiss

  ## Optionally one-hot encode raw multi-level categorical columns so a single
  ## parameter function can be written against the dummy columns.
  if (!is.null(oneHot)) {
    ## Resolve the encoding (reference level and dummy columns) once, against
    ## dfCovs, and reuse it for dfRefRow. Re-deriving it from dfRefRow would use
    ## only the levels that row happens to carry - a single-row reference has
    ## one, so it would get no dummy columns at all while keeping the raw column
    ## dfCovs has just dropped, and the parameter function would silently
    ## compute the baseline as if every dummy were absent.
    oneHotSpec <- normalizeOneHotSpec(oneHot, dfCovs, iMiss, FALSE)
    dfCovs <- oneHotEncode(dfCovs,
      spec = oneHotSpec, sep = oneHotSep,
      missVal = iMiss, dropOriginal = TRUE
    )
    if (!is.null(dfRefRow) && is.data.frame(dfRefRow)) {
      ## A reference row built by setupDfRefRow() may already carry the dummy
      ## columns and not the raw one; encoding is then a no-op bar the warning.
      dfRefRow <- suppressWarnings(
        oneHotEncode(dfRefRow,
          spec = oneHotSpec, sep = oneHotSep,
          missVal = iMiss, dropOriginal = TRUE
        )
      )
    }
  }

  groupnames <- NULL # Store the temp groupnames
  if (any(names(dfCovs) == "COVARIATEGROUPS")) {
    groupnames <- dfCovs[, "COVARIATEGROUPS"]
    dfCovs[, "COVARIATEGROUPS"] <- NULL
  }

  ## After the one-hot encoding, so the dummy columns are compared and not the
  ## raw covariate, and with COVARIATEGROUPS gone.
  checkCovariateNames(functionList, names(dfCovs))
  if (!is.null(dfRefRow) && is.data.frame(dfRefRow)) checkRefRowNames(names(dfRefRow), names(dfCovs))

  ## Define the grouping
  getGroups <- function(df) {
    cGroups <- c()
    cUnique <- c()
    iGroup <- 0
    for (i in seq_len(nrow(df))) {
      tmp <- paste0(names(dfCovs[i, , drop = FALSE])[as.numeric(dfCovs[i, , drop = FALSE]) != iMiss], collapse = ",")
      if (tmp %in% cUnique) {
        tmpl <- which(tmp == cUnique)
        cGroups <- c(cGroups, tmpl)
      } else {
        iGroup <- iGroup + 1
        cGroups <- c(cGroups, iGroup)
        cUnique <- c(cUnique, tmp)
      }
    }
    return(cGroups)
  }
  if (is.null(cGrouping)) cGrouping <- getGroups(dfCovs)


  ## Register to allow for paralell computing. Tear the cluster down on exit
  ## (including on error) rather than only at the end of a successful run.
  if (ncores > 1) {
    stopWorkers <- startParallelWorkers(ncores, pkgs = unique(c("PMXForest", cstrPackages)))
    on.exit(stopWorkers(), add = TRUE)
  }

  ## Calculate the parameters. Result rows per parameter vector k: one per
  ## (dfCovs row) x (element returned by the functionList). Fill typed column
  ## vectors and build the data.frame once - a per-cell data.frame() +
  ## bind_rows() loop is ~50-80x slower for closed-form parameter functions.
  nOut <- length(functionListName)
  internalCalc <- function(k) {
    thetas <- as.numeric(dfParameters[k, 1:noBaseThetas])
    nRow <- nrow(dfCovs) * nOut
    ITER <- integer(nRow)
    COVS <- integer(nRow)
    NAME <- character(nRow)
    VALUE <- numeric(nRow)
    VALUEBASE <- numeric(nRow)
    p <- 0L

    for (i in seq_len(nrow(dfCovs))) {
      n <- 1L
      for (j in seq_along(functionList)) {
        val <- functionList[[j]](thetas = thetas, df = dfCovs[i, , drop = FALSE], ...)
        if (!is.null(dfRefRow)) {
          indi <- min(i, nrow(dfRefRow))
          valbase <- functionList[[j]](thetas = thetas, df = dfRefRow[indi, , drop = FALSE], ...)
        } else {
          dfMissing <- as.data.frame(dfCovs[1, , drop = FALSE])
          dfMissing[, ] <- iMiss
          valbase <- functionList[[j]](thetas = thetas, df = dfMissing, ...)
        }
        for (l in seq_along(val)) {
          p <- p + 1L
          ITER[p] <- k
          COVS[p] <- i
          NAME[p] <- functionListName[n]
          VALUE[p] <- val[[l]]
          VALUEBASE[p] <- valbase[[l]]
          n <- n + 1L
        }
      }
    }
    if (p != nRow) { # a functionList element returned length != 1 somewhere
      idx <- seq_len(p)
      ITER <- ITER[idx]
      COVS <- COVS[idx]
      NAME <- NAME[idx]
      VALUE <- VALUE[idx]
      VALUEBASE <- VALUEBASE[idx]
    }
    data.frame(
      ITER = ITER, COVS = COVS, NAME = NAME, VALUE = VALUE,
      VALUEBASE = VALUEBASE, stringsAsFactors = FALSE
    )
  }

  if (ncores > 1) {
    parts <- foreach(
      k = seq_len(nrow(dfParameters)), .packages = cstrPackages,
      ## Bundle the whole local environment for PSOCK workers (Windows). foreach's
      ## static global detection does not reliably follow `internalCalc`'s free
      ## variables; `cstrExports` remains available for anything outside this frame.
      .export = c(ls(environment()), cstrExports),
      .verbose = !quiet
    ) %dopar% {
      internalCalc(k)
    }
  } else {
    parts <- vector("list", nrow(dfParameters))
    for (k in seq_len(nrow(dfParameters))) parts[[k]] <- internalCalc(k)
  }
  dfres <- bind_rows(parts)

  getCovNameString <- function(dfrow) {
    strName <- ""
    colnames <- names(dfrow)
    for (i in seq_len(ncol(dfrow))) {
      if (dfrow[1, i] != iMiss) {
        if (strName == "") {
          strName <- paste0(colnames[i], "=", dfrow[1, i])
        } else {
          strName <- paste0(strName, ", ", colnames[i], "=", dfrow[1, i])
        }
      }
    }
    if (strName == "") {
      strName <- "Ref cov"
    }
    return(strName)
  }


  dfret <- data.frame()
  for (i in seq_len(nrow(dfCovs))) {
    if (is.null(cdfCovsNames)) {
      covname <- getCovNameString(dfCovs[i, , drop = FALSE])
    } else {
      covname <- cdfCovsNames[i]
    }
    group <- cGrouping[i]
    for (j in seq_along(functionListName)) {
      dft <- dfres[dfres$COVS == i & dfres$NAME == functionListName[j], ]

      quant <- quantile(dft$VALUE, probs = probs, names = FALSE, na.rm = TRUE)
      # Calculate the point value of the forest plot
      FUNCVAL <- pointFunction(dft$VALUE)

      # Define the relative without parameter uncertainty
      dft$RELINTERNAL <- dft$VALUE / dft$VALUEBASE
      # Get quantile and pointvalue when uncertainty is not taken into account
      quantrel <- quantile(dft$RELINTERNAL, probs = probs, names = FALSE, na.rm = TRUE)
      FUNCNOVAR <- pointFunction(dft$RELINTERNAL)

      # Calculate reference value based one the pointFunction
      func_base <- pointFunction(dft$VALUEBASE)
      true_base <- dft$VALUEBASE[dft$ITER == 1]

      # Quantiles of the value relative to the reference. Computed on the ratio
      # directly (rather than dividing `quant` by the scalar base) so the endpoints
      # stay ordered as lower/upper even when the reference value is negative.
      # Fall back to quant/base when the base is NA so an all-NA reference yields
      # NA columns rather than erroring in quantile().
      quant_reffunc <- if (is.na(func_base)) {
        quant / func_base
      } else {
        quantile(dft$VALUE / func_base, probs = probs, names = FALSE, na.rm = TRUE)
      }
      quant_reffinal <- if (is.na(true_base)) {
        quant / true_base
      } else {
        quantile(dft$VALUE / true_base, probs = probs, names = FALSE, na.rm = TRUE)
      }
      groupname <- group
      if (!is.null(groupnames)) groupname <- groupnames[i]
      dfrow <- cbind(dfCovs[i, , drop = FALSE], data.frame(
        GROUP = group,
        GROUPNAME = groupname,
        COVNUM = i, COVNAME = covname, PARAMETER = functionListName[j],
        REFFUNC = func_base, REFFINAL = true_base, POINT = FUNCVAL, POINT_NOVAR_REL_REFFUNC = FUNCNOVAR,
        POINT_REL_REFFUNC = FUNCVAL / func_base, POINT_REL_REFFINAL = FUNCVAL / true_base,
        COVEFF = !all(dft$RELINTERNAL == 1)
      ))
      for (k in seq_along(probs)) {
        dfp <- data.frame(X1 = 1)
        dfp[[paste0("Q", k)]] <- quant[k]
        dfp[[paste0("Q", k, "_REL_REFFUNC")]] <- quant_reffunc[k]
        dfp[[paste0("Q", k, "_REL_REFFINAL")]] <- quant_reffinal[k]
        dfrow <- cbind(dfrow, dfp[, 2:4])
      }
      for (k in seq_along(probs)) {
        dfp <- data.frame(X1 = 1)
        dfp[[paste0("Q", k, "_NOVAR_REL_REFFUNC")]] <- quantrel[k]
        dfrow <- cbind(dfrow, dfp[, 2])
        names(dfrow)[ncol(dfrow)] <- paste0("Q", k, "_NOVAR_REL_REFFUNC")
      }
      dfret <- rbind(dfret, dfrow)
    }
  }

  ## (cluster teardown is handled by on.exit() registered above)

  ## Add a column with YES/NO depending on if refRow was provided or if it was set to the default NULL
  dfret <- dfret %>% mutate(REFROW = ifelse(is.null(dfRefRow), "NO", "YES"))

  ## Make sure GROUPNAME, COVNAME, PARAMETER are factors
  dfret$GROUPNAME <- factor(dfret$GROUPNAME, levels = unique(dfret$GROUPNAME))
  dfret$COVNAME <- factor(dfret$COVNAME, levels = unique(dfret$COVNAME))
  dfret$PARAMETER <- factor(dfret$PARAMETER, levels = unique(dfret$PARAMETER))

  return(dfret)
}

#' Check the covariate columns of dfCovs against generated parameter functions
#'
#' A function made by [createParamFunction()] carries the columns it reads as
#' the attribute `pmxCovariates`. A name in `dfCovs` spelled differently from
#' `$PK`, such as `GENO_1` for `GENO1`, gives two differences at once: a column
#' that no function reads, so its rows have no effect, and a covariate that a
#' function reads and `dfCovs` lacks, so it stays at the reference of the
#' function on every row. Only that pair is reported, because either difference
#' alone is often a deliberate choice of rows. The check is made only when
#' every function carries the list.
#'
#' @param functionList The parameter functions.
#' @param covNames The column names of `dfCovs` after one-hot encoding.
#' @return `NULL`, invisibly. Warns with class `pmxCovariateMismatch`, whose
#'   fields `unread` and `absent` hold the names.
#' @noRd
checkCovariateNames <- function(functionList, covNames) {
  reads <- lapply(functionList, attr, which = "pmxCovariates", exact = TRUE)
  ## A function without the list may read any column, so no column can be
  ## called unread, and the pair cannot be found.
  if (any(vapply(reads, is.null, logical(1)))) {
    return(invisible(NULL))
  }
  readAll <- unique(unlist(reads))
  unread <- setdiff(covNames, readAll)
  absent <- setdiff(readAll, covNames)
  ## Either alone is an extra row or a deliberate subset; a misspelled name
  ## always gives both.
  if (length(unread) == 0 || length(absent) == 0) {
    return(invisible(NULL))
  }

  msg <- c(
    paste0(
      "The covariate columns of dfCovs and the covariates the parameter function ",
      "reads differ in both directions, as they do when a name in dfCovs is ",
      "spelled differently from $PK."
    ),
    paste0(
      "No parameter function reads these columns of dfCovs: ",
      paste(unread, collapse = ", "),
      ". Rows on which only these columns are set have no effect."
    ),
    paste0(
      "dfCovs has no column for these covariates the parameter function reads: ",
      paste(absent, collapse = ", "),
      ". These covariates stay at the reference value of the function on every row."
    ),
    paste0(
      "Check the separator of indicator columns (sep in setupDfCovs(), oneHotSep ",
      "in getForestDFSCM()), which gives GENO_1 against GENO1, and the reference ",
      "level (catRef in setupDfCovs()), which gives GENO2 against GENO1. If both ",
      "lists are intended, the warning can be ignored."
    )
  )
  warning(structure(
    class = c("pmxCovariateMismatch", "warning", "condition"),
    list(
      message = paste(msg, collapse = "\n"), call = NULL,
      unread = unread, absent = absent
    )
  ))
  invisible(NULL)
}

#' Warn when dfRefRow has a column that dfCovs lacks
#'
#' A column of `dfRefRow` that `dfCovs` lacks is either a name spelled
#' differently from `dfCovs`, which leaves the reference subject at the
#' reference of the parameter function for that covariate, or a covariate set
#' in the reference subject and in no row. Either way every row of the Forest
#' plot is shifted, unless the value is the reference of the function. The
#' check needs no list of the columns a function reads, so it is made for
#' hand-written functions too. A column of `dfCovs` that `dfRefRow` lacks is
#' left at the reference of the function in both, and is not reported.
#'
#' @param refNames The column names of `dfRefRow` after one-hot encoding.
#' @param covNames The column names of `dfCovs` after one-hot encoding and
#'   without `COVARIATEGROUPS`.
#' @return `NULL`, invisibly. Warns with class `pmxCovariateMismatch`, whose
#'   field `refOnly` holds the names.
#' @noRd
checkRefRowNames <- function(refNames, covNames) {
  ## setupDfRefRow() copies COVARIATEGROUPS from dfCovs, where it has already
  ## been taken out.
  refOnly <- setdiff(refNames, c(covNames, "COVARIATEGROUPS"))
  if (length(refOnly) == 0) {
    return(invisible(NULL))
  }

  msg <- c(
    paste0("dfRefRow has columns that dfCovs lacks: ", paste(refOnly, collapse = ", "), "."),
    paste0(
      "With a name spelled differently from dfCovs, such as wt for WT, the reference ",
      "subject is held at the reference of the parameter function for that covariate, ",
      "and every row of the Forest plot is shifted."
    ),
    paste0(
      "A covariate set in dfRefRow and in no row of dfCovs has a different value in ",
      "the reference subject than in every row, unless the value is the reference of ",
      "the parameter function."
    ),
    paste0(
      "setupDfRefRow() builds dfRefRow with the columns of dfCovs. If the columns are ",
      "intended, the warning can be ignored."
    )
  )
  warning(structure(
    class = c("pmxCovariateMismatch", "warning", "condition"),
    list(message = paste(msg, collapse = "\n"), call = NULL, refOnly = refOnly)
  ))
  invisible(NULL)
}
