#' Apply a Model's $DATA IGNORE and ACCEPT Statements to a Data Frame
#'
#' @description Reads the `IGNORE` and `ACCEPT` statements from a NONMEM control
#'   stream's `$DATA` record and applies them to a data frame, returning the rows
#'   the model actually used. Summarising the raw analysis file instead of the
#'   filtered one gives covariate quantiles and reference values for a population
#'   the model never saw.
#'
#' @details
#'   **Columns are matched by position, not by name.** NONMEM skips the header
#'   line and takes its column names from `$INPUT`, so the names in the data file
#'   need not agree with the names the model uses - and where they disagree,
#'   filtering by name reads the wrong column. This function maps the first
#'   `length($INPUT)` columns of `data` onto the `$INPUT` names before evaluating
#'   any condition. Columns beyond that are ignored, as NONMEM ignores them, and
#'   `DROP`/`SKIP` columns still occupy a position. A `SYNONYM=REAL` pair can be
#'   referred to by either name.
#'
#'   **What is supported.** `IGNORE=(list)`, `IGNORE(list)`, `ACCEPT=(list)` and
#'   `ACCEPT(list)`, with the `.EQ.`/`.NE.`/`.GT.`/`.GE.`/`.LT.`/`.LE.` operator
#'   family, their `.EQN.`/`.NEN.` variants, and a bare `=` meaning equality.
#'   Conditions separated by commas are combined with OR, as are multiple
#'   statements. A record is dropped if it matches any `IGNORE` condition, and
#'   kept only if it matches some `ACCEPT` condition. NONMEM does not allow both
#'   forms in one `$DATA` record and neither does this function.
#'
#'   **What is not.** The single-character form (`IGNORE=@`, `IGNORE=C`) is a rule
#'   about the first non-blank character of the raw record rather than a condition
#'   on the data, so it cannot be applied to a data frame. `IGNORE=@` is the usual
#'   way of skipping a header line, which `read.csv()` has already done, so this
#'   is normally harmless - a warning is issued and the statement is skipped.
#'
#' @param data A data frame holding at least the first `length($INPUT)` columns of
#'   the model's data file, in that order.
#' @param modFile Path to the NONMEM control stream.
#' @param useInputNames Logical. If `FALSE` (default), the filtered rows are
#'   returned with the caller's own column names. If `TRUE`, the result is the
#'   data as NONMEM sees it: the first `length($INPUT)` columns only, renamed to
#'   the `$INPUT` names.
#' @param quiet Logical. If `FALSE` (default), reports the condition applied and
#'   how many records and subjects it removed.
#' @param idVar The subject identifier, used only for that report. Defaults to
#'   `"ID"`.
#'
#' @return A data frame containing the rows the model used.
#'
#' @seealso [createParamFunction()], [setupDfCovs()], [getCovStats()]
#'
#' @export
#'
#' @examples
#' modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")
#' dfData  <- read.csv(
#'   system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
#' )
#'
#' ## run7.mod drops BLQ records, TYPE 2 records and one subject entirely
#' dfUsed <- filterByModel(dfData, modFile)
#' c(raw = length(unique(dfData$ID)), used = length(unique(dfUsed$ID)))
#'
#' ## Covariate quantiles differ between the two populations
#' getCovStats(dfData, "CRCL", idVar = "ID")
#' getCovStats(dfUsed, "CRCL", idVar = "ID")
#'
#' ## The names $INPUT uses are not the names in this csv file
#' head(names(filterByModel(dfData, modFile, useInputNames = TRUE, quiet = TRUE)), 10)
filterByModel <- function(data, modFile, useInputNames = FALSE, quiet = FALSE,
                          idVar = "ID") {

  data <- as.data.frame(data)
  mod  <- nmReadModel(modFile)
  pos  <- nmInputPositions(mod)
  nInp <- length(pos$names)

  if (nInp == 0) {
    stop("No $INPUT record found in ", basename(modFile),
         ", so the data columns cannot be matched to the model.", call. = FALSE)
  }
  if (ncol(data) < nInp) {
    stop("$INPUT declares ", nInp, " columns but `data` has only ", ncol(data),
         ". Columns are matched by position, so the data frame must hold at ",
         "least the columns the model reads.", call. = FALSE)
  }
  if (ncol(data) > nInp && !quiet) {
    message("$INPUT declares ", nInp, " columns; the ", ncol(data) - nInp,
            " beyond that are not read by the model.")
  }

  ## The model's view of the data: first nInp columns, under the $INPUT names.
  work <- data[, seq_len(nInp), drop = FALSE]
  names(work) <- pos$names
  for (alias in names(pos$aliases)) work[[alias]] <- work[[pos$aliases[[alias]]]]

  keep <- nmDataFilter(mod, work, modFile, quiet)

  if (!quiet) {
    nrRec <- nrow(data) - sum(keep)
    nrSub <- if (idVar %in% names(data)) {
      length(unique(data[[idVar]])) - length(unique(data[[idVar]][keep]))
    } else NA_integer_
    message("Removed ", nrRec, " of ", nrow(data), " record(s)",
            if (!is.na(nrSub)) paste0(" and ", nrSub, " subject(s)") else "", ".")
  }

  if (useInputNames) work[keep, , drop = FALSE] else data[keep, , drop = FALSE]
}

#' Evaluate a $DATA record's IGNORE / ACCEPT statements
#'
#' Returns a logical vector, one per row of `work`, `TRUE` for rows the model
#' keeps. `work` must already carry the `$INPUT` names.
#'
#' @noRd
nmDataFilter <- function(mod, work, modFile, quiet) {
  rec <- nmRecord(mod, "\\$DAT(A)?\\b")
  if (nrow(rec) == 0) return(rep(TRUE, nrow(work)))

  txt <- gsub("\\s+", " ", paste(rec$code, collapse = " "))

  ## The single-character form is a rule about the raw record, not the data, so
  ## it cannot be applied here. "@" and "#" are the conventional header markers
  ## and read.csv() has already dealt with them, so only warn about the rest.
  m <- regmatches(txt, regexpr("IGNORE\\s*=\\s*['\"]?([^\\s(=])['\"]?", txt,
                               perl = TRUE, ignore.case = TRUE))
  if (length(m) == 1) {
    ch <- sub("(?i)^IGNORE\\s*=\\s*['\"]?", "", m, perl = TRUE)
    ch <- substr(ch, 1, 1)
    if (!ch %in% c("@", "#")) {
      warning("The $DATA record of ", basename(modFile), " has IGNORE=", ch,
              ", which applies to the first character of each raw record and ",
              "cannot be applied to a data frame. It has been skipped; filter ",
              "those records yourself if they matter.", call. = FALSE)
    }
  }

  grab <- function(kw) {
    m <- gregexpr(paste0(kw, "\\s*=?\\s*\\(([^)]*)\\)"), txt, ignore.case = TRUE)
    hits <- regmatches(txt, m)[[1]]
    vapply(hits, function(h) sub(paste0("(?i)^", kw, "\\s*=?\\s*\\((.*)\\)$"),
                                 "\\1", h, perl = TRUE), character(1),
           USE.NAMES = FALSE)
  }
  acceptLists <- grab("ACCEPT")
  ignoreLists <- grab("IGNORE")

  if (length(acceptLists) > 0 && length(ignoreLists) > 0) {
    stop("An ACCEPT=(list) and an IGNORE=(list) cannot both appear in the ",
         "$DATA record of ", basename(modFile), ".", call. = FALSE)
  }
  lists <- if (length(acceptLists) > 0) acceptLists else ignoreLists
  if (length(lists) == 0) return(rep(TRUE, nrow(work)))

  ## Commas separate alternatives within a list, and several statements are
  ## alternatives too, so everything is combined with OR.
  conds <- unlist(strsplit(lists, ",", fixed = TRUE))
  conds <- trimws(conds)
  conds <- conds[nzchar(conds)]

  rExprs <- vapply(conds, nmConditionToR, character(1), modFile = modFile,
                   USE.NAMES = FALSE)

  ## Every referenced column must exist and be comparable.
  used <- unique(unlist(lapply(conds, nmConditionSymbols, modFile = modFile)))
  missing <- setdiff(used, names(work))
  if (length(missing) > 0) {
    stop("The $DATA filter in ", basename(modFile), " refers to ",
         paste(missing, collapse = ", "),
         ", which $INPUT does not declare.", call. = FALSE)
  }
  chr <- used[vapply(used, function(u) is.character(work[[u]]), logical(1))]
  if (length(chr) > 0) {
    stop("Column(s) ", paste(chr, collapse = ", "), " used by the $DATA filter ",
         "read as text rather than numbers. Check the data file for non-numeric ",
         "placeholders before filtering.", call. = FALSE)
  }

  full <- paste0("(", rExprs, ")", collapse = " | ")
  if (!quiet) {
    message(if (length(acceptLists) > 0) "Applying ACCEPT: " else "Applying IGNORE: ",
            full)
  }

  hit <- eval(parse(text = full), envir = work)
  hit[is.na(hit)] <- FALSE

  if (length(acceptLists) > 0) hit else !hit
}

#' Translate one $DATA condition into R
#'
#' Uses the same lexer and expression parser as the `$PK` translation, so there
#' is one NONMEM condition grammar in the package rather than two.
#'
#' @noRd
nmConditionToR <- function(cond, modFile) {
  nmDeparse(nmParseCondExpr(cond, modFile))
}

#' @noRd
nmConditionSymbols <- function(cond, modFile) {
  syms <- character(0)
  walk <- function(node) {
    switch(node$type,
      sym   = syms <<- c(syms, node$name),
      call  = lapply(node$args, walk),
      unop  = walk(node$arg),
      binop = { walk(node$lhs); walk(node$rhs) },
      NULL
    )
    invisible(NULL)
  }
  walk(nmParseCondExpr(cond, modFile))
  unique(syms)
}

#' Parse one $DATA condition into an expression tree
#'
#' In `$DATA` a bare `=` means equality, so it is promoted to `==` before lexing;
#' there are no assignments in this position for that to be confused with.
#'
#' @noRd
nmParseCondExpr <- function(cond, modFile) {
  cond <- gsub("(?<![<>!=])=(?!=)", "==", cond, perl = TRUE)
  p <- nmParser(nmLex(cond, 1L, modFile), 1L, modFile)
  e <- nmParseExpr(p)
  if (p$pos <= length(p$toks)) {
    stop("Could not parse the $DATA condition '", cond, "' in ",
         basename(modFile), ".", call. = FALSE)
  }
  e
}
