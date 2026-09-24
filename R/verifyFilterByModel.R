#' Check `filterByModel()` against NONMEM's own account of the data
#'
#' @description Every NONMEM run prints how many records, subjects and
#'   observations it read once the `$DATA` record had been applied. Those three
#'   numbers are an independent check on [filterByModel()]: they come from
#'   NONMEM, not from a reading of the control stream, so they catch a filter
#'   that is wrong in the same way your own reading of `$DATA` is wrong.
#'
#'   It needs nothing but artefacts already on disk - the control stream, its
#'   `.lst`, and the data file - and never re-runs anything.
#'
#' @details The observation count follows NONMEM, using the rule described
#'   under **Subjects without observations** in [filterByModel()]: `MDV == 0`
#'   where there is an `MDV` column, otherwise `EVID == 0`, otherwise records
#'   with no non-zero `AMT`, `RATE` or `SS`. Columns are read by position under
#'   the `$INPUT` names.
#'
#'   The `.lst` describes the model NONMEM read, which is not always the model
#'   on disk. Where its echoed `$INPUT` or `$DATA` filter clauses disagree with
#'   `modFile`, the counts no longer describe it and the comparison is noise, so
#'   this warns. Columns are matched by position, so a changed `$INPUT` matters
#'   even when the names are all still there.
#'
#'   A `$DATA` record that removes nothing agrees with the `.lst` whatever
#'   [filterByModel()] does. That is reported as `attr(., "informative")` rather
#'   than hidden, because such a result is not evidence of anything.
#'
#' @param modFile Path to the control stream.
#' @param data The data set as a data frame. Optional: when absent the `$DATA`
#'   record is resolved relative to `modFile` and read.
#' @param lstFile Path to the run's `.lst`. Optional: `<model>.lst`, `.res`,
#'   `.out` and `NM_run1/psn.lst` are tried in turn.
#' @param idVar Name of the subject column, for the subject count.
#' @param quiet Suppress the progress messages from [filterByModel()].
#'
#' @return A length-1 logical, usable directly in an `if`, with attributes:
#'   \itemize{
#'     \item `checks` - a data frame of `CHECK`, `FILTERED`, `NONMEM`, `DIFF`
#'       and `PASS`, one row per comparison.
#'     \item `informative` - `FALSE` when the `$DATA` record removed nothing, so
#'       the comparison could not have failed.
#'     \item `removed` - how many records the `$DATA` record removed, measured
#'       against NONMEM's own count (`rawRows` minus `NO. OF DATA RECS`) and so
#'       independent of what [filterByModel()] did. On a failing check it will
#'       not reconcile with the filtered count, which is the point.
#'     \item `rawRows` - rows in the data set as read, before filtering.
#'     \item `obsBasis` - the column(s) observations were identified from:
#'       `"MDV"`, `"EVID"`, the dose items present (e.g. `"AMT/RATE"`), or
#'       `"NONE"` when every record is an observation.
#'   }
#'
#' @seealso [filterByModel()] for the filtering itself, and
#'   [verifyParamFunction()], which does the same job for a generated parameter
#'   function.
#'
#' @examples
#' modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")
#' v <- verifyFilterByModel(modFile, quiet = TRUE)
#' v
#' if (v) message("filterByModel() reproduced the run")
#' attr(v, "checks")
#' @export
verifyFilterByModel <- function(modFile, data = NULL, lstFile = NULL,
                                idVar = "ID", quiet = FALSE) {
  if (!file.exists(modFile)) {
    stop("Model file not found: ", modFile, call. = FALSE)
  }
  if (is.null(lstFile)) lstFile <- nmFindLst(modFile)
  if (is.null(lstFile) || !file.exists(lstFile)) {
    stop("No .lst found beside ", basename(modFile),
      ". NONMEM's record of what it read is the whole check here; ",
      "supply lstFile.",
      call. = FALSE
    )
  }

  want <- nmLstCounts(lstFile)
  if (is.na(want$records)) {
    stop(basename(lstFile), " does not carry the record counts, so there is ",
      "nothing to compare against.",
      call. = FALSE
    )
  }
  if (!is.na(want$problems) && want$problems > 1L) {
    warning(basename(lstFile), " covers ", want$problems, " $PROBLEMs; the ",
      "counts below are the first one's.",
      call. = FALSE
    )
  }
  agree <- nmCheckLstMatches(modFile, lstFile)

  if (is.null(data)) data <- nmReadModelData(modFile)

  ## Where the model and the run agree on the column count and the data file
  ## does not, the data is the odd one out. filterByModel() would refuse a
  ## moment later, but pointing at the data frame reads as "pass a wider one"
  ## when what happened is that this file is not the one the run used - and
  ## every position after the gap is shifted, so the columns that do line up
  ## are not the ones they appear to be.
  nInput <- attr(agree, "positions")
  if (!is.null(nInput) && isTRUE(agree) && ncol(data) < nInput) {
    stop(basename(modFile), " and the run it produced both declare ", nInput,
      " columns in $INPUT; this data file has ", ncol(data),
      ". It is not the data the run used. Columns are matched by position, so ",
      "everything after the ", nInput - ncol(data),
      "-column gap would be read from the wrong place - a =DROP item counts, ",
      "so removing one from the data set shifts the rest.",
      call. = FALSE
    )
  }

  used <- filterByModel(data, modFile, quiet = quiet, idVar = idVar)

  obs <- nmObsCount(used, modFile)
  got <- c(
    RECORDS = nrow(used),
    SUBJECTS = if (idVar %in% names(used)) length(unique(used[[idVar]])) else NA_integer_,
    OBSERVATIONS = obs$n
  )
  ref <- c(
    RECORDS = want$records, SUBJECTS = want$ids,
    OBSERVATIONS = want$obs
  )

  checks <- data.frame(
    CHECK = names(got), FILTERED = as.numeric(got), NONMEM = as.numeric(ref),
    DIFF = as.numeric(got) - as.numeric(ref),
    PASS = as.numeric(got) == as.numeric(ref),
    stringsAsFactors = FALSE, row.names = NULL
  )

  removed <- nrow(data) - want$records
  structure(
    isTRUE(all(checks$PASS[!is.na(checks$PASS)])) &&
      !all(is.na(checks$PASS)),
    class = "pmxFilterVerify", checks = checks,
    informative = isTRUE(removed > 0), removed = removed,
    rawRows = nrow(data), obsBasis = obs$basis
  )
}

#' Print a filter verification
#'
#' @param x A `pmxFilterVerify` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.pmxFilterVerify <- function(x, ...) {
  d <- attr(x, "checks")
  cat(if (isTRUE(unclass(x)[1])) "PASS" else "FAIL",
    " - verifyFilterByModel: ", sum(d$PASS, na.rm = TRUE), "/",
    sum(!is.na(d$PASS)), " check(s)",
    if (!isTRUE(attr(x, "informative"))) {
      " - but $DATA removed nothing, so this proves nothing"
    } else {
      paste0(" (", attr(x, "removed"), " record(s) removed)")
    },
    "\n",
    sep = ""
  )
  print(d, row.names = FALSE)
  invisible(x)
}

## The three numbers NONMEM prints about the data it read.
##
## Fortran carriage control puts a stray "0" in front of some labels, so these
## are deliberately not anchored at the start of the line.
##
## @noRd
nmLstCounts <- function(lstFile, maxLines = 8000L) {
  L <- readLines(lstFile, n = maxLines, warn = FALSE)
  num <- function(pat) {
    hit <- grep(pat, L, value = TRUE)
    if (!length(hit)) {
      return(NA_integer_)
    }
    suppressWarnings(as.integer(sub(".*?:\\s*([0-9]+).*", "\\1", hit[1])))
  }
  list(
    records = num("NO\\. OF DATA RECS IN DATA SET:"),
    obs = num("TOT\\. NO\\. OF OBS RECS:"),
    ids = num("TOT\\. NO\\. OF INDIVIDUALS:"),
    problems = length(grep("^\\s*PROBLEM NO\\.:", L)),
    lines = L
  )
}

## Where a run's .lst usually is.
## @noRd
nmFindLst <- function(modFile) {
  dir <- dirname(modFile)
  base <- sub("\\.[^.]*$", "", basename(modFile))
  cand <- c(
    file.path(dir, paste0(base, c(".lst", ".res", ".out", ".LST"))),
    file.path(dir, "NM_run1", "psn.lst")
  )
  cand <- cand[file.exists(cand)]
  if (length(cand)) cand[1] else NULL
}

## NONMEM's observation records in the filtered data, read by position under
## the $INPUT names - the file's own names need not agree.
## @noRd
nmObsCount <- function(used, modFile) {
  mod <- nmReadModel(modFile)
  pos <- nmInputPositions(mod)
  inp <- used[, seq_along(pos$names), drop = FALSE]
  names(inp) <- pos$names
  obs <- nmObsRecords(mod, inp, pos)
  list(n = sum(obs$obs), basis = obs$basis)
}

## Does the .lst describe the control stream on disk?
##
## NM-TRAN echoes its input at the top of the .lst. If $INPUT or the
## IGNORE/ACCEPT clauses have moved since the run, the counts describe a
## different model and comparing against them says nothing. $PK may differ
## freely - it has no bearing on which records are read.
##
## @noRd
nmCheckLstMatches <- function(modFile, lstFile) {
  L <- readLines(lstFile, n = 20000L, warn = FALSE)
  p <- grep("^\\s*\\$PRO", L)
  if (!length(p)) {
    warning("The control stream echo in ", basename(lstFile),
      " could not be read, so it is not known whether the counts describe ",
      basename(modFile), ".",
      call. = FALSE
    )
    return(invisible(structure(NA, positions = NULL)))
  }
  ## Where the echo ends. Not every .lst carries "NM-TRAN MESSAGES" - a real
  ## one does not - so take the earliest of several markers that begin NONMEM's
  ## own output, and fall back to the last $-record line rather than to a fixed
  ## window, which would truncate a long control stream mid-record.
  marks <- c(
    "NM-TRAN MESSAGES", "WARNINGS AND ERRORS", "CREATING MUMODEL",
    "^\\s*License", "^1NONLINEAR MIXED EFFECTS MODEL",
    "^\\s*PROBLEM NO\\.:", "DATA CHECKOUT RUN"
  )
  hits <- unlist(lapply(marks, function(m) grep(m, L)))
  hits <- hits[hits > p[1]]
  lastRec <- grep("^\\s*\\$", L)
  lastRec <- lastRec[lastRec >= p[1]]
  end <- if (length(hits)) {
    min(hits) - 1L
  } else if (length(lastRec)) {
    min(length(L), max(lastRec) + 20L)
  } else {
    length(L)
  }
  echo <- L[p[1]:end]

  modRec <- nmReadModel(modFile)
  ## nmReadModel() reads a path, so the echoed block goes through a temp file.
  echoFile <- tempfile(fileext = ".mod")
  on.exit(unlink(echoFile), add = TRUE)
  writeLines(echo, echoFile)
  echoRec <- nmReadModel(echoFile)

  ## Compare positions, not names. A DROP item still occupies a data column -
  ## the NONMEM guide's own example has DAT1=DROP as the second column of the
  ## file - and nmInputNames() strips them, so a $INPUT that gained or lost a
  ## =DROP entry since the run looks identical by name while every position
  ## after it has moved.
  modPos <- nmInputPositions(modRec)$names
  echoPos <- nmInputPositions(echoRec)$names
  if (!identical(modPos, echoPos)) {
    extra <- if (length(modPos) != length(echoPos)) {
      paste0(
        " It declares ", length(modPos), " column(s) where the run read ",
        length(echoPos), "."
      )
    } else {
      ""
    }
    warning("The $INPUT in ", basename(modFile), " is not the $INPUT the run ",
      "read, as echoed in ", basename(lstFile), ".", extra,
      " Columns are matched by position, so the counts in the .lst describe a ",
      "different reading of the data file. Note a =DROP item still occupies a ",
      "column, so adding or removing one shifts everything after it.",
      call. = FALSE
    )
    return(invisible(structure(FALSE, positions = length(modPos))))
  }
  filterOf <- function(rec) {
    txt <- paste(nmRecord(rec, "\\$DAT(A)?\\b")$code, collapse = " ")
    m <- gregexpr("(IGNORE|ACCEPT)\\s*=?\\s*(\\([^)]*\\)|[^ \t]+)", txt,
      ignore.case = TRUE
    )
    toupper(gsub("\\s+", "", sort(unlist(regmatches(txt, m)))))
  }
  if (!identical(filterOf(modRec), filterOf(echoRec))) {
    warning("The $DATA filter in ", basename(lstFile), " is not the one in ",
      basename(modFile), ", so the counts describe a different subset.",
      call. = FALSE
    )
    return(invisible(structure(FALSE, positions = length(modPos))))
  }
  invisible(structure(TRUE, positions = length(modPos)))
}

## Read the data file a control stream points at.
## @noRd
nmReadModelData <- function(modFile) {
  mod <- nmReadModel(modFile)
  rec <- nmRecord(mod, "\\$DAT(A)?\\b")
  if (nrow(rec) == 0) {
    stop("No $DATA record in ", basename(modFile), ".", call. = FALSE)
  }
  txt <- trimws(sub("^\\s*\\$[A-Za-z]+", "", paste(rec$code, collapse = " ")))
  raw <- gsub("^['\"]|['\"]$", "", strsplit(txt, "\\s+")[[1]][1])
  cand <- c(
    raw, file.path(dirname(modFile), raw),
    file.path(dirname(modFile), basename(raw))
  )
  cand <- cand[file.exists(cand) & !dir.exists(cand)]
  if (!length(cand)) {
    stop("The $DATA file of ", basename(modFile), " (", raw,
      ") is not on disk. Pass the data frame through `data`.",
      call. = FALSE
    )
  }
  head2 <- readLines(cand[1], n = 2L, warn = FALSE)
  sep <- if (length(strsplit(head2[1], ",")[[1]]) > 1) "," else ""
  utils::read.table(cand[1],
    sep = sep, header = TRUE, comment.char = "",
    stringsAsFactors = FALSE
  )
}
