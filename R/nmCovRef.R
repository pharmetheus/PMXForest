## Derivation of covariate reference values from the control stream, and
## emission of the R source. Internal; see `createParamFunction()`.
##
## A generated parameter function is evaluated on an all-`missVal` row whenever
## `getForestDFSCM()` is called without `dfRefRow`, so every covariate needs a
## reference value. Those values are taken from the NONMEM code rather than from
## the data, using four rules of decreasing confidence (see `nmCovRef()`).

## ---------------------------------------------------------------------------
## Reference derivation
## ---------------------------------------------------------------------------

#' Derive a reference value for each covariate
#'
#' Returns a named list; each element is
#' `list(value = <numeric>, source = <character>, line = <integer>,
#'       confident = <logical>)`. Covariates for which no rule fires are absent
#' from the result and are reported by the caller.
#'
#' @noRd
nmCovRef <- function(stmts, covariates, missVal) {
  flat <- nmFlatten(stmts)
  out  <- list()

  for (cov in covariates) {
    r <- nmRefExplicit(flat, cov, missVal)
    if (is.null(r)) r <- nmRefMarkedBranch(flat, cov)
    if (is.null(r)) r <- nmRefIdentityBranch(flat, cov)
    if (is.null(r)) r <- nmRefNormalisation(flat, cov)
    if (is.null(r)) r <- nmRefUntestedLevel(flat, cov)
    if (!is.null(r)) out[[cov]] <- r
  }
  out
}

#' Numeric value of a literal node, or `NA`
#'
#' A negative literal parses as a unary minus applied to a number, so the rules
#' below cannot simply test `type == "num"`. Handling that here keeps them
#' correct whether or not constant folding has already run.
#'
#' @noRd
nmAsNumber <- function(node) {
  if (is.null(node)) return(NA_real_)
  if (node$type == "num") return(node$value)
  if (node$type == "unop" && node$arg$type == "num") {
    return(if (node$op == "-") -node$arg$value else node$arg$value)
  }
  NA_real_
}

#' Flatten nested statements into a list of records for pattern matching
#'
#' Each record is `list(kind, cond, stmt, lineno, comment)` where `kind` is
#' `"assign"` for an unconditional assignment and `"cond-assign"` for one guarded
#' by a condition.
#'
#' @noRd
nmFlatten <- function(stmts, cond = NULL) {
  out <- list()
  for (s in stmts) {
    if (s$type == "assign") {
      out[[length(out) + 1L]] <- list(
        kind = if (is.null(cond)) "assign" else "cond-assign",
        cond = cond, stmt = s, lineno = s$lineno, comment = s$comment
      )
    } else {
      out <- c(out, nmFlatten(s$then, s$cond))
      for (e in s$elifs) out <- c(out, nmFlatten(e$stmts, e$cond))
      if (!is.null(s$else_)) out <- c(out, nmFlatten(s$else_, NULL))
    }
  }
  out
}

#' Rule 1: explicit handling, `IF(X.EQ.-99) X = <literal>`
#' @noRd
nmRefExplicit <- function(flat, cov, missVal) {
  for (f in flat) {
    if (f$kind != "cond-assign" || f$stmt$lhs != cov) next
    eq <- nmEqualityTest(f$cond, cov)
    if (is.null(eq) || !isTRUE(all.equal(eq, missVal))) next
    val <- nmAsNumber(f$stmt$rhs)
    if (is.na(val)) next
    return(list(value = val, line = f$lineno, confident = TRUE,
                source = paste0("explicit missing-value handling in the control stream")))
  }
  NULL
}

#' Rule 2a: the branch PsN's scm marks with "Most common"
#' @noRd
nmRefMarkedBranch <- function(flat, cov) {
  for (f in flat) {
    if (f$kind != "cond-assign") next
    if (!grepl("most\\s*common", f$comment, ignore.case = TRUE)) next
    eq <- nmEqualityTest(f$cond, cov)
    if (is.null(eq)) next
    return(list(value = eq, line = f$lineno, confident = TRUE,
                source = "reference level of the \";  Most common\" branch"))
  }
  NULL
}

#' Rule 2b: the branch whose right-hand side is the literal identity
#'
#' In an SCM block one branch assigns a bare `1` (multiplicative) or `0`
#' (additive); that branch is the reference category.
#'
#' @noRd
nmRefIdentityBranch <- function(flat, cov) {
  for (f in flat) {
    if (f$kind != "cond-assign") next
    val <- nmAsNumber(f$stmt$rhs)
    if (is.na(val) || !(val %in% c(0, 1))) next
    eq <- nmEqualityTest(f$cond, cov)
    if (is.null(eq)) next
    return(list(value = eq, line = f$lineno, confident = TRUE,
                source = paste0("branch assigning the identity value ",
                                nmFormatNum(val))))
  }
  NULL
}

#' Rule 3: the normalisation constant in `(X/c)` or `(X-c)`
#' @noRd
nmRefNormalisation <- function(flat, cov) {
  hit <- NULL
  scan <- function(node, lineno) {
    if (!is.null(hit)) return(invisible(NULL))
    if (node$type == "binop") {
      const <- nmAsNumber(node$rhs)
      if (node$op %in% c("/", "-") &&
          node$lhs$type == "sym" && node$lhs$name == cov && !is.na(const)) {
        hit <<- list(
          value = const, line = lineno, confident = TRUE,
          source = paste0("normalisation constant in (", cov, " ", node$op, " ",
                          nmFormatNum(const), ")")
        )
        return(invisible(NULL))
      }
      scan(node$lhs, lineno); scan(node$rhs, lineno)
    } else if (node$type == "unop") {
      scan(node$arg, lineno)
    } else if (node$type == "call") {
      for (a in node$args) scan(a, lineno)
    }
    invisible(NULL)
  }
  for (f in flat) {
    scan(f$stmt$rhs, f$lineno)
    if (!is.null(hit)) return(hit)
  }
  NULL
}

#' Rule 4: a level never tested by any equality on the covariate
#'
#' For a dummy referenced only as `IF(X.EQ.1) ...` with an unconditional default
#' assignment above it, the reference is the untested level. Proposed, not
#' trusted: the caller warns and asks for confirmation.
#'
#' @noRd
nmRefUntestedLevel <- function(flat, cov) {
  tested <- c()
  lineno <- NA_integer_
  for (f in flat) {
    if (f$kind != "cond-assign") next
    eq <- nmEqualityTest(f$cond, cov)
    if (is.null(eq)) next
    tested <- c(tested, eq)
    if (is.na(lineno)) lineno <- f$lineno
  }
  if (length(tested) == 0) return(NULL)

  for (candidate in c(0, 1)) {
    if (!any(vapply(tested, function(t) isTRUE(all.equal(t, candidate)), logical(1)))) {
      return(list(
        value = candidate, line = lineno, confident = FALSE,
        source = paste0("level not tested by any IF() on ", cov,
                        " (proposed, please confirm)")
      ))
    }
  }
  NULL
}

#' If `cond` is an equality test of `cov` against a literal, return the literal
#'
#' Recurses through `&` so that `IF(SEX.EQ.1.AND.STUDY.EQ.2)` still yields 1 for
#' SEX. Returns `NULL` when there is no such test.
#'
#' @noRd
nmEqualityTest <- function(cond, cov) {
  if (is.null(cond)) return(NULL)
  if (cond$type == "binop") {
    if (cond$op == "==") {
      if (cond$lhs$type == "sym" && cond$lhs$name == cov) {
        v <- nmAsNumber(cond$rhs)
        if (!is.na(v)) return(v)
      }
      if (cond$rhs$type == "sym" && cond$rhs$name == cov) {
        v <- nmAsNumber(cond$lhs)
        if (!is.na(v)) return(v)
      }
      return(NULL)
    }
    if (cond$op == "&") {
      l <- nmEqualityTest(cond$lhs, cov)
      if (!is.null(l)) return(l)
      return(nmEqualityTest(cond$rhs, cov))
    }
  }
  NULL
}

## ---------------------------------------------------------------------------
## Emitter
## ---------------------------------------------------------------------------

#' Emit the R source for a parameter function
#'
#' @noRd
nmEmit <- function(stmts, covRef, covariates, parameters, functionName,
                   modFile, missVal, noBaseThetas) {
  base <- basename(modFile)
  ind  <- "  "
  L    <- character(0)
  add  <- function(...) L <<- c(L, ...)

  add(paste0("## Generated by PMXForest::createParamFunction() from ", base, "."),
      "## Typical values: every ETA() in $PK has been set to 0.",
      "## Review this against the control stream before use.",
      "",
      paste0(functionName, " <- function(thetas, df, ...) {"))

  ## -- covariate reference preamble ------------------------------------------
  if (length(covariates) > 0) {
    add("",
        paste0(ind, "## ---- Covariate references, taken from ", base,
               " ", strrep("-", max(0, 46 - nchar(base)))))
    width <- max(nchar(covariates))
    for (cov in covariates) {
      r   <- covRef[[cov]]
      loc <- if (is.na(r$line)) "" else paste0("  [", base, ":", r$line, "]")
      add(paste0(ind, "## ", formatC(cov, width = width, flag = "-"), "  ",
                 r$source, loc))
      pad <- formatC(cov, width = width, flag = "-")
      add(paste0(ind, pad, " <- if (!is.null(df$", cov, ") && df$", cov,
                 " != ", nmFormatNum(missVal), ") df$", cov, " else ",
                 nmFormatNum(r$value)))
    }
  }

  ## -- transliterated $PK ----------------------------------------------------
  add("",
      paste0(ind, "## ---- $PK, transliterated ", strrep("-", 51)))
  add(nmEmitStmts(stmts, ind, base))

  ## -- extension point -------------------------------------------------------
  add("",
      paste0(ind, "## ---- Secondary parameters: add yours below ", strrep("-", 33)),
      paste0(ind, "## Quantities such as AUC, Cmax or event probabilities are not"),
      paste0(ind, "## derivable from $PK and are left to you, for example:"),
      paste0(ind, "##   AUC <- 80 / (CL / FREL)"))

  ## -- return ----------------------------------------------------------------
  add("",
      paste0(ind, "list("),
      paste0(ind, ind,
             paste(paste0(parameters, " = ", parameters), collapse = ",\n    ")),
      paste0(ind, ")"),
      "}")

  add("",
      paste0("## functionListName <- c(",
             paste(paste0('"', parameters, '"'), collapse = ", "), ")"),
      paste0("## noBaseThetas     <- ", noBaseThetas))

  L
}

#' Emit a statement list at a given indent
#' @noRd
nmEmitStmts <- function(stmts, ind, base) {
  out <- character(0)
  for (s in stmts) {
    if (s$type == "assign") {
      etaNote <- if (isTRUE(s$hadEta)) "  # ETA() -> 0 (typical value)" else ""
      out <- c(out, paste0(ind, s$lhs, " <- ", nmDeparse(s$rhs), etaNote))
    } else if (isTRUE(s$oneline)) {
      # Keep a one-line IF on one line, so the source reads like the $PK block.
      inner   <- s$then[[1]]
      etaNote <- if (isTRUE(inner$hadEta)) "  # ETA() -> 0 (typical value)" else ""
      out <- c(out, paste0(ind, "if (", nmDeparse(s$cond), ") ",
                           inner$lhs, " <- ", nmDeparse(inner$rhs), etaNote))
    } else {
      out <- c(out, paste0(ind, "if (", nmDeparse(s$cond), ") {"))
      out <- c(out, nmEmitStmts(s$then, paste0(ind, "  "), base))
      for (e in s$elifs) {
        out <- c(out, paste0(ind, "} else if (", nmDeparse(e$cond), ") {"))
        out <- c(out, nmEmitStmts(e$stmts, paste0(ind, "  "), base))
      }
      if (!is.null(s$else_)) {
        out <- c(out, paste0(ind, "} else {"))
        out <- c(out, nmEmitStmts(s$else_, paste0(ind, "  "), base))
      }
      out <- c(out, paste0(ind, "}"))
    }
  }
  out
}

#' @noRd
nmHasEta <- function(node) {
  switch(node$type,
    eta   = TRUE,
    call  = any(vapply(node$args, nmHasEta, logical(1))),
    unop  = nmHasEta(node$arg),
    binop = nmHasEta(node$lhs) || nmHasEta(node$rhs),
    FALSE
  )
}
