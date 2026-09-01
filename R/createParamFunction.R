#' Generate a Parameter Function from a NONMEM Control Stream
#'
#' @description Translates the `$PK` block of a NONMEM control stream into R
#'   source code for a `paramFunction` of the form
#'   `function(thetas, df, ...)`, the shape required by `getForestDFSCM()` and
#'   `getForestDFemp()`. The result is **source text for you to read, check and
#'   edit** - nothing is evaluated. Writing this function by hand duplicates code
#'   that already exists in the control stream and is an easy place to introduce
#'   a silent error.
#'
#' @details
#'   **Typical values.** Every `ETA(n)` in `$PK` is set to 0, so the generated
#'   function returns typical values, which is what a Forest plot needs.
#'
#'   **Covariate references.** `getForestDFSCM()` evaluates the parameter
#'   function on a row where every covariate equals `missVal` whenever it is
#'   called without a `dfRefRow`, so each covariate needs a reference value. All
#'   of that handling is hoisted into a single preamble block at the top of the
#'   generated function - one line per covariate, each annotated with where its
#'   reference came from - after which the `$PK` algebra is a straight
#'   transliteration. The references are taken from the control stream, using
#'   four rules in decreasing order of confidence:
#'   \enumerate{
#'     \item explicit handling in the code, `IF(WT.EQ.-99) WT = 75`;
#'     \item the branch that PsN's scm marks with `; Most common`;
#'     \item the branch whose right-hand side is the identity value (`1` for a
#'       multiplicative term, `0` for an additive one);
#'     \item the normalisation constant in `(WT/75)` or `(AGE-50)`.
#'   }
#'   A fifth, weaker rule proposes the level that no `IF()` on the covariate
#'   tests, and **warns**; confirm it or override it through `covRef`. A
#'   covariate that matches no rule is an error - supply its reference through
#'   `covRef`. The generator never guesses in silence.
#'
#'   **Secondary parameters are not generated.** AUC, Cmax, event probabilities
#'   and anything else reached through `$ERROR` or `$DES` are not derivable from
#'   `$PK` - the dose used for an AUC, for instance, appears nowhere in the
#'   control stream. The emitted source carries a marked extension point for
#'   you to add them.
#'
#'   **Accepted syntax.** Assignments, one-line `IF(...) VAR = ...`,
#'   `IF/ELSE IF/ELSE/END IF` blocks, arithmetic (`**` becomes `^`), the
#'   `.EQ.`/`.AND.` operator family, and the usual intrinsic functions.
#'   Anything else - `$DES`, compartment amounts `A(n)`, verbatim FORTRAN, `DO`
#'   loops, `CALL` - raises an error naming the file and line.
#'
#' @param modFile Path to the NONMEM control stream (`.mod` or `.ctl`).
#' @param parameters A character vector of `$PK` variables the function should
#'   return. Defaults to `NULL`, meaning every variable assigned in `$PK`; trim
#'   it to the parameters you actually want to plot.
#' @param covRef An optional named list of covariate reference values, e.g.
#'   `list(WT = 75)`. Overrides the values derived from the control stream, and
#'   supplies them for covariates where no rule fires.
#' @param functionName The name given to the generated function. Defaults to
#'   `"paramFunction"`.
#' @param extFile An optional path to the model's `.ext` file. When supplied, the
#'   THETA count is read from its header, which is authoritative, instead of
#'   being counted from the `$THETA` records.
#' @param file An optional path to write the generated source to. The source is
#'   returned either way.
#' @param missVal The value marking an inactive covariate in `dfCovs`. Defaults
#'   to -99, matching `getForestDFSCM()`.
#' @param quiet Logical. If `FALSE` (default), reports the covariates found and
#'   the reference value chosen for each.
#'
#' @return A list of four elements:
#'   \itemize{
#'     \item `code` - the generated R source, a character vector of lines with
#'       class `"pmxParamFunction"` so that printing it renders the source.
#'     \item `functionListName` - a character vector matching the order of the
#'       returned parameters, for `getForestDFSCM()`.
#'     \item `noBaseThetas` - the number of THETAs in the model.
#'     \item `covRef` - the reference value used for each covariate, with the
#'       rule and control-stream line it came from.
#'     \item `etaMap` - for each returned parameter written as
#'       `P = <expr> * EXP(ETA(n))`, the ETA index `n`. Used by
#'       [verifyParamFunction()] to recover typical values from a NONMEM table.
#'     \item `modFile`, `missVal` - as supplied.
#'   }
#'
#' @seealso [verifyParamFunction()] to check the generated function against a
#'   NONMEM `$TABLE`; [getForestDFSCM()], [getForestDFemp()].
#'
#' @export
#'
#' @examples
#' modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")
#'
#' ## Generate and inspect the source
#' out <- createParamFunction(modFile, parameters = c("CL", "FREL", "V"))
#' cat(out$code, sep = "\n")
#'
#' ## Where each covariate reference came from
#' out$covRef$WT
#' out$noBaseThetas
#'
#' ## Read it in and use it like any other parameter function
#' paramFunction <- eval(parse(text = out$code))
#' paramFunction(thetas = rep(1, out$noBaseThetas),
#'               df     = data.frame(WT = 90, FOOD = 0))
createParamFunction <- function(modFile, parameters = NULL, covRef = NULL,
                                functionName = "paramFunction", extFile = NULL,
                                file = NULL, missVal = -99, quiet = FALSE) {

  mod <- nmReadModel(modFile)
  pk  <- nmRecord(mod, "\\$PK\\b")
  if (nrow(pk) == 0) {
    stop("No $PK record found in ", basename(modFile),
         ". createParamFunction() translates $PK blocks; a $PRED model must be ",
         "written by hand.", call. = FALSE)
  }

  stmts <- nmParseStatements(pk, modFile)
  if (length(stmts) == 0) {
    stop("The $PK record in ", basename(modFile), " contains no statements.",
         call. = FALSE)
  }
  # Recorded before folding, while the ETA() references are still in the tree.
  etaMap <- nmEtaMap(stmts)
  # Setting ETA() to 0 leaves artefacts such as `TVCL * exp(0)`; fold them away
  # so the emitted source stays diffable against the control stream.
  stmts <- nmSimplifyStmts(stmts)

  ## Covariates: named in $INPUT and never assigned in $PK.
  syms       <- nmSymbols(stmts)
  inputNames <- nmInputNames(mod)
  covariates <- intersect(syms$used, setdiff(inputNames, syms$assigned))
  covariates <- covariates[order(match(covariates, syms$used))]

  ## References, control stream first, then user overrides.
  derived <- nmCovRef(stmts, covariates, missVal)
  for (cov in names(covRef)) {
    derived[[cov]] <- list(value = covRef[[cov]], line = NA_integer_,
                           confident = TRUE, source = "supplied through covRef")
  }

  missingRef <- setdiff(covariates, names(derived))
  if (length(missingRef) > 0) {
    stop("No reference value could be derived from ", basename(modFile),
         " for: ", paste(missingRef, collapse = ", "),
         ".\nSupply them through covRef, e.g. covRef = list(",
         paste(paste0(missingRef, " = <value>"), collapse = ", "), ").",
         call. = FALSE)
  }

  unsure <- covariates[!vapply(derived[covariates], `[[`, logical(1), "confident")]
  if (length(unsure) > 0) {
    warning("The reference value for ", paste(unsure, collapse = ", "),
            " was inferred rather than read from the control stream. Check the ",
            "preamble of the generated function and override through covRef if ",
            "it is wrong.", call. = FALSE)
  }

  ## Parameters to return.
  if (is.null(parameters)) {
    parameters <- syms$assigned
  } else {
    unknown <- setdiff(parameters, syms$assigned)
    if (length(unknown) > 0) {
      stop("Not assigned in the $PK block of ", basename(modFile), ": ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
  }

  ## THETA count: the .ext header when available, else the $THETA records.
  noBaseThetas <- if (!is.null(extFile)) {
    nms <- names(getExt(extFile))
    length(grep("^THETA", nms))
  } else {
    nmCountThetas(mod)
  }
  maxTheta <- nmMaxTheta(stmts)
  if (noBaseThetas < maxTheta) {
    stop("The model declares ", noBaseThetas, " THETA(s) but $PK references ",
         "THETA(", maxTheta, ") in ", basename(modFile),
         ". Refusing to continue: the theta indices would be wrong.",
         call. = FALSE)
  }

  code <- nmEmit(stmts, derived, covariates, parameters, functionName,
                 modFile, missVal, noBaseThetas)
  class(code) <- c("pmxParamFunction", "character")

  if (!is.null(file)) writeLines(code, file)

  if (!quiet) {
    message("Translated $PK of ", basename(modFile), ": ", length(stmts),
            " statement(s), ", length(covariates), " covariate(s), ",
            noBaseThetas, " theta(s).")
    for (cov in covariates) {
      message("  ", cov, " reference ", nmFormatNum(derived[[cov]]$value),
              " - ", derived[[cov]]$source)
    }
    if (!is.null(file)) message("Written to ", file)
  }

  list(code = code, functionListName = parameters,
       noBaseThetas = noBaseThetas, covRef = derived,
       etaMap = etaMap[intersect(names(etaMap), parameters)],
       modFile = modFile, missVal = missVal)
}

#' Print generated parameter-function source
#'
#' @param x A `pmxParamFunction` object, the `code` element of a
#'   [createParamFunction()] result.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @export
print.pmxParamFunction <- function(x, ...) {
  cat(unclass(x), sep = "\n")
  invisible(x)
}
