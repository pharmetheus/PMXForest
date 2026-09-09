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
#'   **Secondary parameters.** AUC, Cmax, event probabilities and anything else
#'   reached through `$ERROR` or `$DES` are not derivable from `$PK` - the dose
#'   used for an AUC, for instance, appears nowhere in the control stream.
#'   Supply them through `secondary`, a named list. Each entry's value is
#'   either a line of R code (`secondary = list(AUC = "df$DOSE / CL")`), the
#'   path to an `.R` file of arbitrary code - including a `deSolve` or
#'   `mrgsolve` simulation (`secondary = list(CMAX = "cmax.R")`) - or a list
#'   carrying that `source` plus constants the code needs
#'   (`secondary = list(CMAX = list(source = "cmax.R", dose = 100, tau = 12))`).
#'   Each entry is spliced into the generated function inside a `local()` block -
#'   so it sees `thetas`, `df`, `...`, any constants passed alongside `source`,
#'   and every structural parameter by name, with covariate columns reached as
#'   `df$NAME` - and its value is added to the return list under the entry's
#'   name. The names also appear in `functionListName`, so `getForestDFSCM()`
#'   picks the secondary parameters up automatically. Entries are evaluated in
#'   order, so a later one may use an earlier one. Without `secondary`, the
#'   emitted source carries a marked extension point for you to add them by
#'   hand instead.
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
#' @param secondary An optional named list of secondary parameters to append to
#'   the generated function's return list. Each entry's value is a single string
#'   (R code whose last value is the result, `list(AUC = "df$DOSE / CL")`, or
#'   the path to an `.R` file of arbitrary code such as an `mrgsolve`
#'   simulation, `list(CMAX = "cmax.R")`), or a list `list(source = <string>,
#'   ...)` carrying that `source` plus named atomic constants the code needs
#'   (`list(CMAX = list(source = "cmax.R", dose = 100, tau = 12))`). A file's
#'   text is inlined into the generated source, so the result stays
#'   self-contained. See Details.
#'
#' @return A list of:
#'   \itemize{
#'     \item `code` - the generated R source, a character vector of lines with
#'       class `"pmxParamFunction"` so that printing it renders the source.
#'     \item `functionListName` - a character vector matching the order of the
#'       returned parameters, for `getForestDFSCM()`. Includes the `secondary`
#'       names, appended after the `$PK` parameters.
#'     \item `primaryNames` - the `$PK` parameters only.
#'     \item `secondaryNames` - the `secondary` parameter names (`character(0)`
#'       when none). [verifyParamFunction()] skips these by default.
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
#'
#' ## Append secondary parameters. An entry is a snippet, an .R file path, or
#' ## a list with `source` plus constants the code needs:
#' out2 <- createParamFunction(
#'   modFile, parameters = c("CL", "V"), quiet = TRUE,
#'   secondary = list(
#'     AUC24 = list(source = "dose / CL", dose = 100),   # constant `dose`
#'     KEL   = "CL / V"))
#' out2$functionListName        # c("CL", "V", "AUC24", "KEL")
#' cat(out2$code, sep = "\n")
createParamFunction <- function(modFile, parameters = NULL, covRef = NULL,
                                functionName = "paramFunction", extFile = NULL,
                                file = NULL, missVal = -99, quiet = FALSE,
                                secondary = NULL) {

  p <- nmParsePK(modFile, parameters = parameters, covRef = covRef,
                 extFile = extFile, missVal = missVal)

  sec      <- nmResolveSecondary(secondary, quiet = quiet)
  secNames <- unname(vapply(sec, `[[`, "", "name"))

  # Setting ETA() to 0 leaves artefacts such as `TVCL * exp(0)`; fold them away
  # so the emitted source stays diffable against the control stream.
  stmts <- nmSimplifyStmts(p$statements)

  code <- nmEmit(stmts, p$covRef, p$covariates, p$parameters, functionName,
                 modFile, missVal, p$noBaseThetas, secondary = sec)
  class(code) <- c("pmxParamFunction", "character")

  if (!is.null(file)) writeLines(code, file)

  if (!quiet) {
    message("Translated $PK of ", basename(modFile), ": ", length(stmts),
            " statement(s), ", length(p$covariates), " covariate(s), ",
            p$noBaseThetas, " theta(s)",
            if (length(secNames))
              paste0(", ", length(secNames), " secondary parameter(s)") else "",
            ".")
    for (cov in p$covariates) {
      message("  ", cov, " reference ", nmFormatNum(p$covRef[[cov]]$value),
              " - ", p$covRef[[cov]]$source)
    }
    if (!is.null(file)) message("Written to ", file)
  }

  list(code = code,
       functionListName = c(p$parameters, secNames),
       primaryNames     = p$parameters,
       secondaryNames   = secNames,
       noBaseThetas = p$noBaseThetas, covRef = p$covRef,
       etaMap = p$etaMap[intersect(names(p$etaMap), p$parameters)],
       modFile = modFile, missVal = missVal)
}

#' Parse a NONMEM `$PK` block into a reusable structure
#'
#' @description The shared front end of [createParamFunction()]: it reads the
#'   control stream, parses `$PK` into a statement tree, and works out the
#'   covariates, their reference values, the THETA count and the ETA that
#'   carries each parameter's between-subject variability. [createParamFunction()]
#'   emits SCM-style parameter-function source from this structure; other
#'   packages emit their own (PMXFrem's FREM parameter functions wrap each
#'   covariate parameter as `TV * exp(covthetas + eta)` instead).
#'
#'   Nothing is evaluated and no source is generated - this returns the parse
#'   only. The statement tree keeps `ETA()` references intact, so a downstream
#'   emitter can decide what to do with them.
#'
#' @inheritParams createParamFunction
#'
#' @return A list:
#'   \itemize{
#'     \item `statements` - the parsed `$PK` statement tree. Each element is an
#'       `assign` (`lhs`, `rhs`, `lineno`, `comment`) or an `if`
#'       (`cond`, `then`, `elifs`, `else_`, `lineno`); `rhs`/`cond` are
#'       expression nodes. Use [nmDeparse()] to render a node as R source.
#'     \item `covariates` - covariate names: `$INPUT` columns that `$PK` reads
#'       before assigning them, in first-use order. Being assigned later does
#'       not disqualify a column, since NONMEM populates the data items before
#'       `$PK` runs - that is what makes `IF(WT.EQ.-99) WT = 75` a covariate
#'       rather than a local.
#'     \item `covRef` - a named list, one entry per covariate, each with
#'       `value`, `line`, `confident` (logical) and `source` (how the reference
#'       was found). `covRef` overrides and user-supplied entries are merged in.
#'     \item `parameters` - the `$PK` variables to return: `parameters` as given,
#'       or every assigned variable when `NULL`.
#'     \item `noBaseThetas` - the THETA count (from the `.ext` header when
#'       `extFile` is supplied, otherwise the `$THETA` records).
#'     \item `etaMap` - a named integer vector: for each parameter written
#'       `P = <expr> * EXP(ETA(n))`, the ETA index `n`.
#'     \item `inputNames` - the `$INPUT` column names.
#'     \item `modFile`, `missVal` - as supplied.
#'   }
#'
#' @seealso [createParamFunction()], [nmDeparse()].
#'
#' @export
#'
#' @examples
#' modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")
#' p <- nmParsePK(modFile, parameters = c("CL", "V"))
#' p$covariates
#' p$noBaseThetas
#' p$etaMap
nmParsePK <- function(modFile, parameters = NULL, covRef = NULL,
                      extFile = NULL, missVal = -99) {

  mod <- nmReadModel(modFile)
  pk  <- nmRecord(mod, "\\$PK\\b")
  if (nrow(pk) == 0) {
    stop("No $PK record found in ", basename(modFile),
         ". nmParsePK() parses $PK blocks; a $PRED model must be ",
         "handled by hand.", call. = FALSE)
  }

  stmts <- nmParseStatements(pk, modFile)
  if (length(stmts) == 0) {
    stop("The $PK record in ", basename(modFile), " contains no statements.",
         call. = FALSE)
  }
  # Recorded while the ETA() references are still in the tree.
  etaMap <- nmEtaMap(stmts)
  # A folded copy (ETA() -> 0, constants collapsed) for the analyses below; the
  # raw tree is what we return so a downstream emitter keeps the ETA()s.
  folded <- nmSimplifyStmts(stmts)

  ## Covariates are discovered, not filtered: a $PK that reads a name before
  ## binding it is reading a data item, which is exactly what NONMEM does. Each
  ## such $INPUT name gets a preamble binding and the walk is repeated, so the
  ## covariates come out in first-use order. Anything read before it is bound
  ## that is *not* in $INPUT cannot be supplied, and is an error here rather
  ## than an "object not found" from inside getForestDFSCM() later.
  syms       <- nmSymbols(folded)
  inputNames <- nmInputNames(mod)
  covariates <- character(0)
  repeat {
    unbound <- nmFirstUnboundUse(folded, covariates, syms$assigned)
    if (is.null(unbound)) break
    if (unbound$name %in% inputNames) {
      covariates <- c(covariates, unbound$name)
      next
    }
    at <- if (is.na(unbound$lineno)) "" else paste0(" on line ", unbound$lineno)
    stop(
      "The $PK block of ", basename(modFile), " reads ", unbound$name, at,
      if (unbound$everAssigned) {
        paste0(
          " before assigning it.\nNONMEM does not initialise $PK variables and",
          " does not clear them between data records, so the model reads",
          " whatever the previous record left there. A parameter function is",
          " evaluated one row at a time and cannot reproduce that.",
          "\nMove the assignment of ", unbound$name, " above the line that",
          " reads it."
        )
      } else {
        paste0(
          ", but never assigns it and it is not in $INPUT.\nIf it is a NONMEM",
          " reserved variable such as NEWIND, it has no value a parameter",
          " function could supply."
        )
      },
      call. = FALSE
    )
  }

  ## References, control stream first, then user overrides.
  derived <- nmCovRef(folded, covariates, missVal)
  if (length(covRef) > 0) {
    if (is.null(names(covRef)) || any(!nzchar(names(covRef)))) {
      stop("Every element of covRef must be named.", call. = FALSE)
    }
    ## A typo here used to be silently ignored, leaving the derived reference
    ## in place, and a non-numeric value emitted source that parsed but failed
    ## when called.
    stray <- setdiff(names(covRef), covariates)
    if (length(stray) > 0) {
      stop("covRef names a covariate that ", basename(modFile),
           " does not use: ", paste(stray, collapse = ", "),
           ".\nCovariates in this $PK: ",
           if (length(covariates) == 0) "none" else paste(covariates, collapse = ", "),
           ".", call. = FALSE)
    }
    bad <- names(covRef)[!vapply(covRef, function(v) {
      is.numeric(v) && length(v) == 1L && !is.na(v) && is.finite(v)
    }, logical(1))]
    if (length(bad) > 0) {
      stop("Each covRef value must be a single finite number; check: ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    }
  }
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
  maxTheta <- nmMaxTheta(folded)
  if (noBaseThetas < maxTheta) {
    stop("The model declares ", noBaseThetas, " THETA(s) but $PK references ",
         "THETA(", maxTheta, ") in ", basename(modFile),
         ". Refusing to continue: the theta indices would be wrong.",
         call. = FALSE)
  }

  list(statements = stmts, covariates = covariates, covRef = derived,
       parameters = parameters, noBaseThetas = noBaseThetas, etaMap = etaMap,
       inputNames = inputNames, modFile = modFile, missVal = missVal)
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
