#' Setup the Expression List for Empirical Forest Plots
#'
#' @description The empirical-workflow counterpart of [setupDfCovs()]. It turns a
#'   data frame plus a vector of covariate names into the named
#'   `covExpressionsList` (and a matching `cdfCovsNames` label vector) consumed by
#'   [getForestDFemp()], reusing the deduplicated level and quantile logic of
#'   [getCovStats()]. Each covariate becomes one or more forest-plot rows, each row
#'   an `expression()` that selects the subjects it summarises.
#'
#' @details
#'   Statistics are computed on **deduplicated baseline data** (one record per
#'   `idVar`), matching [setupDfCovs()] and [setupDfRefRow()]. Rows equal to
#'   `missVal`, and `NA`, are dropped per covariate before any calculation.
#'
#'   **Continuous covariates** (more than `minLevels` unique values) produce two
#'   rows. With `contSplit = "quantile"` (default) they are `cov < q_low` and
#'   `cov >= q_high`, where `q_low` / `q_high` are the `probs` quantiles rounded to
#'   `nsig` significant digits; subjects between the two quantiles fall in neither
#'   row (a deliberate tails-versus-reference contrast). With
#'   `contSplit = "median"` the split is `cov < m` / `cov >= m` at the median, so
#'   every subject lands in exactly one row.
#'
#'   **Binary covariates** produce one row per level, `cov == lo` and `cov == hi`.
#'
#'   **Multi-level categorical covariates** (`minLevels` or fewer unique values,
#'   more than two) produce one `cov == level` row per level. With
#'   `includeReference = TRUE` (default) every level gets a row, since in an
#'   empirical plot each level is a real subject subset. With
#'   `includeReference = FALSE` the reference level is dropped (from
#'   `catRef[[cov]]` if named, otherwise the lowest level), mirroring the
#'   parametric `dfCovs`.
#'
#'   **Additional covariates.** Each entry of `additionalCovs` gets its own rows,
#'   exactly like a primary covariate, appended after the primary blocks. In
#'   addition, a fixed condition on that covariate is combined with `&` into every
#'   other covariate's expression: a level for a categorical covariate
#'   (`list(FOOD = 1)` gives `FOOD == 1`), or a `prob` / `value` split with a
#'   direction for a continuous covariate
#'   (`list(CRCL = list(prob = 0.5, dir = "gt"))` gives `CRCL > q_0.5`,
#'   `list(AGE = list(value = 65, dir = "lt"))` gives `AGE < 65`). A primary
#'   covariate's rows carry every additional-covariate condition; an additional
#'   covariate's own rows carry every *other* additional-covariate condition, not
#'   its own.
#'
#'   **Subject-count safeguard.** Every generated expression is evaluated against
#'   the deduplicated `data`. If any selects fewer than `minSubjects` subjects the
#'   function stops, so no empty or near-empty subset reaches [getForestDFemp()].
#'   The check uses the `data` given here; if a different or subset `dfData` is
#'   later passed to [getForestDFemp()], re-verify.
#'
#'   **Labels.** `cdfCovsNames` holds terse generated labels (`"WT <63.2"`,
#'   `"SEX 1"`, `"GENO 2"`). Publication plots normally override them with domain
#'   labels; the additional-covariate conditions are not repeated in the labels.
#'
#'   For a fixed reference row, build one with [setupDfCovs()] and
#'   [setupDfRefRow()] (`singleRef = TRUE`) and pass it to
#'   `getForestDFemp(dfRefRow = )`, exactly as for a parametric forest plot.
#'
#' @param data A data frame that includes the covariates to summarise. Only the
#'   first record per subject (identified by `idVar`) is used.
#' @param covariates A character vector of primary covariate names to evaluate.
#' @param additionalCovs An optional named list, one entry per additional
#'   covariate. A scalar value is the level for a categorical covariate; a list
#'   `list(dir = "lt"/"gt", prob = )` or `list(dir = "lt"/"gt", value = )` is the
#'   split for a continuous covariate (supply exactly one of `prob` or `value`;
#'   `prob` must lie in `(0, 1)`). Defaults to `NULL`. A name may not also appear
#'   in `covariates`.
#' @param contSplit How to split continuous covariates into two rows. `"quantile"`
#'   (default) uses the `probs` quantile tails; `"median"` splits at the median.
#' @param includeReference Logical. If `TRUE` (default), a multi-level categorical
#'   covariate emits a row for every level. If `FALSE`, the reference level is
#'   dropped.
#' @param minSubjects The minimum number of subjects an expression may select in
#'   the deduplicated `data`. The function stops if any expression selects fewer.
#'   Default is 10.
#' @param probs A numeric vector of two probabilities used to split continuous
#'   covariates and, when requested, to place a continuous additional-covariate
#'   condition. Defaults to `c(0.05, 0.95)`.
#' @param minLevels The maximum number of unique values a covariate can have to be
#'   treated as categorical. Default is 10.
#' @param idVar The name of the subject identifier column. Defaults to `"ID"`.
#' @param missVal The numeric value indicating missing data, excluded per
#'   covariate before any calculation. Defaults to -99.
#' @param nsig The number of significant digits for rounding continuous split
#'   thresholds. Defaults to 3.
#' @param catRef The reference level for multi-level categorical covariates,
#'   e.g. `list(GENO = 2)`; also accepts `"lowest"` (the default), `"mode"` or
#'   `"model"`, either as a single setting or in a named list with a `default`
#'   component. Only consulted when `includeReference = FALSE`. Replaces
#'   `refLevels`.
#' @param model The NONMEM control stream, required when `catRef` is `"model"`.
#'   Either a path to the `.mod` file or the list returned by
#'   [createParamFunction()].
#' @param refLevels Deprecated. Use `catRef`.
#'
#' @return A list of two elements, named for the [getForestDFemp()] arguments they
#'   feed:
#'   \itemize{
#'     \item `covExpressionsList` -- a named list of `expression()` objects, one
#'       per forest-plot row; the names are the covariate names and drive the
#'       `GROUPNAME` grouping.
#'     \item `cdfCovsNames` -- a character vector of row labels, the same length as
#'       `covExpressionsList`.
#'   }
#'
#' @seealso [getForestDFemp()], [getCovStats()], [setupDfCovs()],
#'   [setupDfRefRow()]
#'
#' @export
#'
#' @examples
#' dfData <- read.csv(
#'   system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
#' )
#'
#' # Continuous -> quantile-tail rows; SEX -> both levels; GENO -> one row per level
#' out <- setupCovExpressionsList(dfData, covariates = c("WT", "SEX", "GENO"),
#'                                idVar = "ID")
#' out$covExpressionsList
#' out$cdfCovsNames
#'
#' # Median split instead of quantile tails
#' setupCovExpressionsList(dfData, covariates = "WT", contSplit = "median",
#'                         idVar = "ID")$covExpressionsList
#'
#' # Drop the reference genotype (level 2) from the GENO rows
#' setupCovExpressionsList(dfData, covariates = "GENO", includeReference = FALSE,
#'                         catRef = list(GENO = 2),
#'                         idVar = "ID")$covExpressionsList
#'
#' # FOOD as a categorical additional covariate and CRCL as a continuous one:
#' # every WT/SEX row is also conditioned on FOOD == 1 and CRCL above its lower quartile
#' setupCovExpressionsList(
#'   dfData, covariates = c("WT", "SEX"), contSplit = "median",
#'   additionalCovs = list(FOOD = 1, CRCL = list(prob = 0.25, dir = "gt")),
#'   idVar = "ID"
#' )$covExpressionsList
#'
#' # Feed the result into getForestDFemp(): pass `covExpressionsList` and
#' # `cdfCovsNames` from the returned list. Each row's POINT is the empirical
#' # median CL over the subjects that row's expression selects.
#' dfDataEmp <- dfData[!duplicated(dfData$ID), ][1:100, ]
#' dfSamples <- getSamples(
#'   system.file("extdata", "SimVal/run7.cov", package = "PMXForest"),
#'   system.file("extdata", "SimVal/run7.ext", package = "PMXForest"),
#'   n = 20
#' )
#' paramFunction <- function(thetas, df, ...) {
#'   TVCL <- thetas[4]
#'   if (df$WT != -99) TVCL <- TVCL * (df$WT / 75)^thetas[2]
#'   if (df$FOOD != -99 && df$FOOD == 0) TVCL <- TVCL * (1 + thetas[11])
#'   list(CL = TVCL)
#' }
#' out <- setupCovExpressionsList(dfData, covariates = c("WT", "FOOD"), idVar = "ID")
#' dfresEmp <- getForestDFemp(
#'   dfData             = dfDataEmp,
#'   covExpressionsList = out$covExpressionsList,
#'   cdfCovsNames       = out$cdfCovsNames,
#'   functionList       = list(paramFunction),
#'   functionListName   = "CL",
#'   noBaseThetas       = 14,
#'   dfParameters       = dfSamples,
#'   ncores             = 1
#' )
#' dfresEmp[, c("COVNAME", "GROUPNAME", "PARAMETER", "POINT", "POINT_REL_REFFUNC")]
setupCovExpressionsList <- function(data, covariates, additionalCovs = NULL,
                                    contSplit = c("quantile", "median"),
                                    includeReference = TRUE, minSubjects = 10,
                                    probs = c(0.05, 0.95), minLevels = 10,
                                    idVar = "ID", missVal = -99, nsig = 3,
                                    catRef = NULL, model = NULL,
                                    refLevels = NULL) {

  catRef <- refLevelsToCatRef(refLevels, catRef, "setupCovExpressionsList")

  contSplit <- match.arg(contSplit)
  data <- as.data.frame(data)

  # ---- validate additionalCovs shape -------------------------------------
  addNames <- names(additionalCovs)
  if (!is.null(additionalCovs)) {
    if (!is.list(additionalCovs) || is.null(addNames) ||
        any(addNames == "") || anyDuplicated(addNames)) {
      stop("`additionalCovs` must be a named list with one entry per covariate, ",
           "e.g. list(FOOD = 1).")
    }
    clash <- intersect(addNames, covariates)
    if (length(clash)) {
      stop("A covariate cannot be both primary and additional: ",
           paste(clash, collapse = ", "), ".")
    }
  }

  allNames <- c(covariates, addNames)
  if (!all(allNames %in% names(data))) {
    stop("Not all covariates are present in the data.")
  }
  if (length(probs) != 2) stop("`probs` must have length 2.")
  if (!is.numeric(minSubjects) || length(minSubjects) != 1 || minSubjects < 1) {
    stop("`minSubjects` must be a single number >= 1.")
  }

  dedup <- data %>% dplyr::distinct(!!rlang::sym(idVar), .keep_all = TRUE)

  # non-missing, NA-safe values for one covariate (deduplicated)
  getVals <- function(cov) {
    v <- dedup[[cov]]
    v <- v[v != missVal & !is.na(v)]
    if (length(v) == 0) stop("Covariate ", cov, " contains only missing values.")
    v
  }

  # covariate type from the unique-value count: continuous | binary | multi | single
  covType <- function(cov) {
    nLev <- length(unique(getVals(cov)))
    if (nLev > minLevels) "continuous"
    else if (nLev == 2)   "binary"
    else if (nLev > 2)    "multi"
    else                  "single"
  }

  # ---- one block of rows for a covariate --------------------------------
  # Returns list(exprs = <list of language>, labels = <character>).
  makeBlock <- function(cov) {
    v    <- getVals(cov)
    type <- covType(cov)

    if (type == "continuous") {
      if (contSplit == "quantile") {
        qs <- signif(stats::quantile(v, probs = probs, names = FALSE), nsig)
      } else {
        m  <- signif(stats::median(v), nsig)
        qs <- c(m, m)
      }
      exprs  <- list(bquote(.(as.name(cov)) <  .(qs[1])),
                     bquote(.(as.name(cov)) >= .(qs[2])))
      labels <- c(paste0(cov, " <", qs[1]), paste0(cov, " >=", qs[2]))

    } else if (type == "single") {
      lev <- sort(unique(v))[1]
      warning("Covariate ", cov, " has only one non-missing level; emitting a ",
              "single row with no contrast.")
      exprs  <- list(bquote(.(as.name(cov)) == .(lev)))
      labels <- paste0(cov, " ", lev)

    } else { # binary or multi
      levs <- sort(unique(v))
      if (is.numeric(levs)) levs <- as.numeric(levs) # drop integer literal suffix
      if (type == "multi" && !includeReference) {
        refLev <- refEncodingLevel(catRef, cov, levs, model, missVal)
        if (is.null(refLev)) refLev <- refMode(v)
        if (!refLev %in% levs) {
          stop("Reference level ", refLev, " for covariate '", cov,
               "' is not present in the data.")
        }
        levs <- setdiff(levs, refLev)
      }
      exprs  <- lapply(levs, function(L) bquote(.(as.name(cov)) == .(L)))
      labels <- paste0(cov, " ", levs)
    }
    list(exprs = exprs, labels = labels)
  }

  # ---- fixed condition fragment for each additional covariate ----------
  fragList <- list()
  for (a in addNames) {
    spec <- additionalCovs[[a]]
    if (covType(a) == "continuous") {
      if (!is.list(spec)) {
        stop("`additionalCovs$", a, "` is a continuous covariate; supply ",
             "list(prob = <p>, dir = \"lt\"/\"gt\") or ",
             "list(value = <v>, dir = \"lt\"/\"gt\").")
      }
      if (is.null(spec$dir) || !spec$dir %in% c("lt", "gt")) {
        stop("`additionalCovs$", a, "$dir` must be \"lt\" or \"gt\".")
      }
      hasProb  <- !is.null(spec$prob)
      hasValue <- !is.null(spec$value)
      if (hasProb == hasValue) {
        stop("`additionalCovs$", a, "` needs exactly one of `prob` or `value`.")
      }
      if (hasProb) {
        if (spec$prob <= 0 || spec$prob >= 1) {
          stop("`additionalCovs$", a, "$prob` must be in (0, 1).")
        }
        p <- signif(stats::quantile(getVals(a), probs = spec$prob, names = FALSE), nsig)
      } else {
        p <- spec$value
      }
      fragList[[a]] <- if (spec$dir == "lt") bquote(.(as.name(a)) <  .(p)) else
                                             bquote(.(as.name(a)) >  .(p))
    } else {
      if (is.list(spec) || length(spec) != 1) {
        stop("`additionalCovs$", a, "` is a categorical covariate; supply a ",
             "single level value.")
      }
      if (!spec %in% getVals(a)) {
        stop("Level ", spec, " for `additionalCovs$", a,
             "` is not present in the data.")
      }
      fragList[[a]] <- bquote(.(as.name(a)) == .(spec))
    }
  }

  # combine a base call with a set of fragments using `&`
  andFrags <- function(base, frags) {
    Reduce(function(x, y) bquote(.(x) & .(y)), frags, base)
  }

  # ---- assemble every block -------------------------------------------
  blocks <- c(
    lapply(covariates, function(cov)
      list(cov = cov, blk = makeBlock(cov), frags = fragList)),
    lapply(addNames, function(a)
      list(cov = a, blk = makeBlock(a), frags = fragList[setdiff(addNames, a)]))
  )

  parts <- lapply(blocks, function(b) {
    exprs <- lapply(b$blk$exprs, function(e) as.expression(andFrags(e, b$frags)))
    list(exprs = exprs, labels = b$blk$labels,
         groups = rep(b$cov, length(exprs)))
  })

  covExpressionsList <- unlist(lapply(parts, `[[`, "exprs"), recursive = FALSE)
  cdfCovsNames       <- unlist(lapply(parts, `[[`, "labels"), use.names = FALSE)
  names(covExpressionsList) <- unlist(lapply(parts, `[[`, "groups"),
                                      use.names = FALSE)

  # ---- subject-count safeguard ---------------------------------------
  for (i in seq_along(covExpressionsList)) {
    n <- nrow(subset(dedup, eval(covExpressionsList[[i]][[1]])))
    if (n < minSubjects) {
      stop("Expression `", as.character(covExpressionsList[[i]]), "` selects ", n,
           " subject(s) in `data`, fewer than `minSubjects` (", minSubjects,
           "). Relax `minSubjects`, widen `probs`, or adjust the split / ",
           "`additionalCovs` conditions.")
    }
  }

  list(covExpressionsList = covExpressionsList, cdfCovsNames = cdfCovsNames)
}
