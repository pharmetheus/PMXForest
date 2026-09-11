## verifyPMXForest.R -----------------------------------------------------------
##
## Stress-test PMXForest against your own NONMEM runs.
##
## PMXForest 1.3.0 adds two pieces of new code that read control streams:
##
##   createParamFunction()  translates a $PK block into R source
##   filterByModel()        applies the $DATA IGNORE/ACCEPT statements
##
## Both have only ever run against a handful of in-house models, and those
## models are nearly all the same shape. This script points them at yours.
##
## It NEVER re-runs a model. Everything is read from artefacts already on disk:
## the control stream, the .lst, the .ext, and any $TABLE output.
##
## It NEVER writes into a run directory. All output goes to one directory you
## name, plus the R session's tempdir().
##
## Usage
## -----
##   Rscript verifyPMXForest.R ~/projects/*/Models
##   Rscript verifyPMXForest.R --install ~/projects/xyz/Models/run7.mod
##
##   source("verifyPMXForest.R")
##   verifyPMXForest("~/projects/xyz/Models")
##
## What to send back
## -----------------
## Three files are written. Send back ONLY the one ending "-report.csv" - it is
## redacted: no paths, no covariate names, no subject counts. The file ending
## "-local.rds" holds everything verbatim and is for your eyes; the script will
## tell you what is in it. Run with --full only if you are content to share the
## detail.
##
## Contact: <maintainer>
## ------------------------------------------------------------------------ ##

PMXF_SCRIPT_VERSION <- "1.0.0"
PMXF_WANT_VERSION   <- "1.2.15.9009"
PMXF_MAX_DATA_BYTES <- 500e6
PMXF_MAX_TABLE_ROWS <- 2e5
PMXF_MODEL_TIMEOUT  <- 300


# -- small utilities ----------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a[1])) b else a

pmxfMsg <- function(...) cat(..., "\n", sep = "")

pmxfRule <- function(title = NULL) {
  if (is.null(title)) {
    cat(strrep("-", 78), "\n")
  } else {
    cat("\n", title, "\n", strrep("-", 78), "\n", sep = "")
  }
}

## FNV-1a, 64 bit, folded to 12 hex characters. No external dependency, and we
## only need collision resistance across a few thousand strings.
pmxfHash <- function(x, salt = "") {
  x <- paste0(salt, paste(x, collapse = "\n"))
  bytes <- as.integer(charToRaw(x))
  h1 <- 2166136261
  h2 <- 16777619
  a <- 0
  b <- 0
  for (byte in bytes) {
    a <- bitwXor(a %% 2^31, byte)
    a <- (a * 31 + 7) %% 2^31
    b <- bitwXor(b %% 2^31, bitwAnd(byte * 131 + a, 255))
    b <- (b * 131 + 17) %% 2^31
  }
  sprintf("%06x%06x", a %% 16^6, b %% 16^6)
}

## Every deliberate stop() in PMXForest passes call. = FALSE, so a deliberate
## refusal has a NULL condition call and a base-R crash does not. We keep that
## signal separately from the message match so the taxonomy can be re-derived
## from a returned CSV without asking anyone to re-run anything.
pmxfSafely <- function(expr, timeout = PMXF_MODEL_TIMEOUT) {
  warns <- character(0)
  msgs <- character(0)
  t0 <- proc.time()[["elapsed"]]
  setTimeLimit(cpu = timeout, elapsed = timeout, transient = TRUE)
  on.exit(setTimeLimit())
  res <- withCallingHandlers(
    tryCatch(list(ok = TRUE, value = force(expr), cond = NULL),
      error = function(e) list(ok = FALSE, value = NULL, cond = e)
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, sub("\n$", "", conditionMessage(m)))
      invokeRestart("muffleMessage")
    }
  )
  res$warnings <- warns
  res$messages <- msgs
  res$seconds <- round(proc.time()[["elapsed"]] - t0, 2)
  res
}


# -- control-stream reading ---------------------------------------------------
#
# Deliberately independent of PMXForest: we must not use the code under test to
# decide whether the code under test is right. These parsers only ROUTE - they
# decide what is feasible. Every quantity that is actually compared comes from
# PMXForest on one side and NONMEM's own .lst on the other.

pmxfReadLines <- function(path, n = -1L) {
  x <- suppressWarnings(readLines(path, n = n, warn = FALSE, skipNul = TRUE))
  if (any(is.na(nchar(x, allowNA = TRUE)))) {
    x <- iconv(x, from = "latin1", to = "UTF-8", sub = "?")
  }
  x <- sub("^﻿", "", x)
  sub("\r$", "", x)
}

pmxfEolKind <- function(path) {
  raw <- readBin(path, "raw", n = min(file.size(path), 65536))
  cr <- sum(raw == as.raw(13))
  lf <- sum(raw == as.raw(10))
  if (cr == 0 && lf > 0) "LF" else if (cr == lf && cr > 0) "CRLF" else if (lf == 0 && cr > 0) "CR" else "MIXED"
}

pmxfStripComments <- function(lines) trimws(sub(";.*$", "", lines), which = "right")

## NONMEM allows record names to be abbreviated to three characters.
pmxfRecordName <- function(line) {
  nm <- toupper(sub("^\\s*\\$([A-Za-z]+).*$", "\\1", line))
  known <- c(
    PRO = "PROBLEM", INP = "INPUT", DAT = "DATA", SUB = "SUBROUTINE",
    MOD = "MODEL", PK = "PK", DES = "DES", ERR = "ERROR", PRE = "PRED",
    THE = "THETA", OME = "OMEGA", SIG = "SIGMA", EST = "ESTIMATION",
    COV = "COVARIANCE", TAB = "TABLE", SIM = "SIMULATION", MIX = "MIX",
    ABB = "ABBREVIATED", AES = "AES", INF = "INFN", MSF = "MSFI",
    SCA = "SCATTERPLOT", PRI = "PRIOR", LEV = "LEVEL"
  )
  key <- substr(nm, 1, 3)
  if (nm == "PK") return("PK")
  if (key %in% names(known)) known[[key]] else nm
}

## Every occurrence of every record, continuations joined. No window limit:
## a $TABLE's FILE= is routinely four or five lines below the token.
pmxfRecords <- function(lines) {
  starts <- grep("^\\s*\\$", lines)
  if (!length(starts)) {
    return(data.frame(record = character(0), lineno = integer(0), text = character(0)))
  }
  ends <- c(starts[-1] - 1L, length(lines))
  data.frame(
    record = vapply(lines[starts], pmxfRecordName, ""),
    lineno = starts,
    text = vapply(seq_along(starts), function(i) {
      paste(trimws(lines[starts[i]:ends[i]]), collapse = " ")
    }, ""),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

pmxfRecordAll <- function(recs, name) recs$text[recs$record == name]

pmxfInputNames <- function(recs) {
  txt <- pmxfRecordAll(recs, "INPUT")
  if (!length(txt)) return(character(0))
  toks <- strsplit(trimws(sub("^\\s*\\$[A-Za-z]+", "", paste(txt, collapse = " "))), "\\s+")[[1]]
  toks[nzchar(toks)]
}


# -- the .lst: NONMEM's own record of what it did -----------------------------

pmxfLstEcho <- function(L) {
  p <- grep("^\\s*\\$PRO", L)
  if (!length(p)) return(character(0))
  stop_at <- grep("NM-TRAN MESSAGES|^License |^ License ", L)
  stop_at <- stop_at[stop_at > p[1]]
  L[p[1]:(if (length(stop_at)) stop_at[1] - 1L else min(length(L), p[1] + 400L))]
}

pmxfParseLst <- function(lstFile, maxLines = 8000L) {
  L <- pmxfReadLines(lstFile, n = maxLines)
  ## Fortran carriage control puts a stray "0" in front of some labels, so
  ## these patterns are deliberately not anchored at the start of the line -
  ## except the ID one, which must not match "EVENT ID DATA ITEM ...".
  num <- function(pat) {
    hit <- grep(pat, L, value = TRUE)
    if (!length(hit)) return(NA_integer_)
    suppressWarnings(as.integer(sub(".*?:\\s*([0-9]+).*", "\\1", hit[1])))
  }
  list(
    found        = TRUE,
    nDataRecs    = num("NO\\. OF DATA RECS IN DATA SET:"),
    nObsRecs     = num("TOT\\. NO\\. OF OBS RECS:"),
    nIndividuals = num("TOT\\. NO\\. OF INDIVIDUALS:"),
    lengthTheta  = num("LENGTH OF THETA:"),
    omegaDim     = num("OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:"),
    idItemNo     = num("^\\s*ID DATA ITEM IS DATA ITEM NO\\.:"),
    nProblems    = length(grep("^\\s*PROBLEM NO\\.:", L)),
    hasSimulation = any(grepl("SIMULATION STEP PERFORMED", L)),
    echo         = pmxfLstEcho(L)
  )
}

pmxfEmptyLst <- function() {
  list(found = FALSE, nDataRecs = NA_integer_, nObsRecs = NA_integer_,
       nIndividuals = NA_integer_, lengthTheta = NA_integer_,
       omegaDim = NA_integer_, idItemNo = NA_integer_, nProblems = NA_integer_,
       hasSimulation = NA, echo = character(0))
}

pmxfNorm <- function(x) toupper(gsub("\\s+", " ", trimws(paste(x, collapse = " "))))

## The single most important guard in the script. If someone edited $PK or
## $INPUT after the run, the .lst counters no longer describe the .mod on disk,
## and Tier A would report a failure against PMXForest that is really a stale
## artefact. The .lst echoes the control stream NM-TRAN actually read, so we
## can check. The $DATA *filename* is allowed to differ - PsN rewrites it.
pmxfEchoAgreement <- function(modLines, echoLines) {
  if (!length(echoLines)) return(list(verdict = "NO_ECHO"))
  mr <- pmxfRecords(pmxfStripComments(modLines))
  er <- pmxfRecords(pmxfStripComments(echoLines))
  inputSame <- identical(pmxfInputNames(mr), pmxfInputNames(er))
  filterOf <- function(recs) {
    txt <- paste(pmxfRecordAll(recs, "DATA"), collapse = " ")
    m <- gregexpr("(IGNORE|ACCEPT)\\s*=?\\s*(\\([^)]*\\)|[^ \t]+)", txt, ignore.case = TRUE)
    pmxfNorm(sort(unlist(regmatches(txt, m))))
  }
  filterSame <- identical(filterOf(mr), filterOf(er))
  pkSame <- identical(pmxfNorm(pmxfRecordAll(mr, "PK")), pmxfNorm(pmxfRecordAll(er, "PK")))
  verdict <- if (!inputSame) {
    "INPUT_DIFF"
  } else if (!filterSame) {
    "FILTER_DIFF"
  } else if (!pkSame) {
    "PK_DIFF"
  } else {
    "MATCH"
  }
  list(verdict = verdict, inputSame = inputSame, filterSame = filterSame, pkSame = pkSame)
}


# -- locating run artefacts ---------------------------------------------------

pmxfFindArtefacts <- function(modFile) {
  dir <- dirname(modFile)
  base <- sub("\\.[^.]*$", "", basename(modFile))
  pick <- function(exts, sub = NULL) {
    cand <- if (is.null(sub)) file.path(dir, paste0(base, exts)) else file.path(dir, sub, paste0("psn", exts))
    cand <- cand[file.exists(cand)]
    if (length(cand)) cand[1] else NA_character_
  }
  lst <- pick(c(".lst", ".res", ".out", ".LST"))
  if (is.na(lst)) lst <- pick(".lst", sub = "NM_run1")
  ext <- pick(c(".ext", ".EXT"))
  if (is.na(ext)) ext <- pick(".ext", sub = "NM_run1")
  list(mod = modFile, dir = dir, base = base, lst = lst, ext = ext)
}

pmxfDataSpec <- function(recs) {
  txt <- paste(pmxfRecordAll(recs, "DATA"), collapse = " ")
  if (!nzchar(txt)) return(list(path = NA_character_))
  rest <- trimws(sub("^\\s*\\$[A-Za-z]+", "", txt))
  path <- strsplit(rest, "\\s+")[[1]][1]
  list(path = gsub("^['\"]|['\"]$", "", path))
}

pmxfResolveData <- function(rawPath, modFile, lstFile) {
  if (is.na(rawPath)) return(list(exists = FALSE, how = "NO_DATA_RECORD"))
  tries <- list(
    c(rawPath, "ABSOLUTE"),
    c(file.path(dirname(modFile), rawPath), "REL_MOD"),
    c(file.path(dirname(modFile), basename(rawPath)), "BASENAME_MOD"),
    c(file.path(dirname(modFile), "..", basename(rawPath)), "PARENT_MOD")
  )
  if (!is.na(lstFile)) {
    tries <- c(tries, list(c(file.path(dirname(lstFile), basename(rawPath)), "BASENAME_LST")))
  }
  for (t in tries) {
    if (file.exists(t[1]) && !dir.exists(t[1])) {
      return(list(exists = TRUE, path = normalizePath(t[1]), how = t[2],
                  bytes = file.size(t[1]), mtime = file.mtime(t[1])))
    }
  }
  list(exists = FALSE, how = "NOT_FOUND")
}

## NONMEM data files come with and without a header, comma- or space-separated,
## and often with leading comment lines marked @ # C or ;.
pmxfSniffData <- function(path, nInput) {
  head20 <- pmxfReadLines(path, n = 20L)
  skip <- 0L
  while (skip < length(head20) &&
         grepl("^\\s*[@#Cc;]", head20[skip + 1L]) &&
         !grepl("^\\s*[Cc][0-9., \t-]*$", head20[skip + 1L])) {
    skip <- skip + 1L
  }
  body <- head20[(skip + 1L):length(head20)]
  body <- body[nzchar(trimws(body))]
  if (!length(body)) return(list(ok = FALSE))
  nf <- function(s, sep) length(strsplit(s, sep)[[1]])
  sep <- if (nf(body[1], ",") > 1 && nf(body[1], ",") >= nf(body[1], "[ \t]+")) "," else ""
  first <- strsplit(trimws(body[1]), if (nzchar(sep)) sep else "[ \t]+")[[1]]
  header <- any(is.na(suppressWarnings(as.numeric(first))))
  list(ok = TRUE, sep = sep, skip = skip, header = header, nCols = length(first))
}

pmxfReadData <- function(path, sniff, inputNames) {
  df <- utils::read.table(path,
    sep = if (nzchar(sniff$sep)) sniff$sep else "",
    header = sniff$header, skip = sniff$skip,
    na.strings = c("NA", ""), comment.char = "", check.names = TRUE,
    stringsAsFactors = FALSE, colClasses = NA
  )
  if (!sniff$header && length(inputNames) && ncol(df) >= length(inputNames)) {
    names(df)[seq_along(inputNames)] <- inputNames
  }
  df
}


# -- Tier A: filterByModel against the .lst -----------------------------------

pmxfObsCount <- function(kept) {
  if ("MDV" %in% names(kept)) {
    return(list(n = sum(kept$MDV == 0, na.rm = TRUE), basis = "MDV"))
  }
  if ("EVID" %in% names(kept)) {
    return(list(n = sum(kept$EVID == 0, na.rm = TRUE), basis = "EVID"))
  }
  list(n = NA_integer_, basis = "NONE")
}

pmxfTierA <- function(ctx, row) {
  skip <- function(reason) {
    row$tierAStatus <- "SKIP"
    row$tierASkipReason <- reason
    row
  }
  if (!ctx$art$lst %in% TRUE && is.na(ctx$art$lst)) return(skip("NO_LST"))
  lst <- ctx$lst
  if (is.na(lst$nDataRecs)) return(skip("LST_NO_COUNTS"))
  if (!is.na(lst$nProblems) && lst$nProblems != 1L) return(skip("MULTIPLE_PROBLEMS"))
  if (!ctx$echo$verdict %in% c("MATCH")) return(skip(paste0("LST_STALE_", ctx$echo$verdict)))

  d <- pmxfResolveData(ctx$dataSpec$path, ctx$modFile, ctx$art$lst)
  if (!isTRUE(d$exists)) return(skip(paste0("DATA_", d$how)))
  if (d$bytes > PMXF_MAX_DATA_BYTES) return(skip("DATA_TOO_LARGE"))
  if (!is.na(ctx$art$lst) && d$mtime > file.mtime(ctx$art$lst)) return(skip("DATA_CHANGED"))

  sn <- pmxfSniffData(d$path, length(ctx$inputNames))
  if (!isTRUE(sn$ok)) return(skip("DATA_UNREADABLE"))
  raw <- pmxfSafely(pmxfReadData(d$path, sn, ctx$inputNames))
  if (!raw$ok) return(skip("DATA_UNREADABLE"))
  rawDf <- raw$value
  row$tierARawRows <- nrow(rawDf)

  ft <- pmxfSafely(PMXForest::filterByModel(rawDf, ctx$modFile, useInputNames = TRUE, quiet = TRUE))
  if (!ft$ok) {
    cls <- pmxfClassifyRefusal(ft$cond, ctx$modFile)
    row$tierAStatus <- if (cls$deliberate) "PASS_REFUSED" else "ERROR"
    row$tierAErrorKind <- cls$kind
    return(row)
  }
  kept <- ft$value
  idCol <- if ("ID" %in% names(kept)) "ID" else names(kept)[lst$idItemNo %||% 1L]
  obs <- pmxfObsCount(kept)

  recPass <- nrow(kept) == lst$nDataRecs
  idPass <- if (is.na(lst$nIndividuals) || is.null(kept[[idCol]])) NA else
    length(unique(kept[[idCol]])) == lst$nIndividuals
  obsPass <- if (is.na(lst$nObsRecs) || obs$basis == "NONE") NA else obs$n == lst$nObsRecs

  row$tierARecPass <- recPass
  row$tierAIdPass <- idPass
  row$tierAObsPass <- obsPass
  row$tierARecDelta <- nrow(kept) - lst$nDataRecs
  row$tierAIdDelta <- if (isTRUE(is.na(idPass))) NA_integer_ else length(unique(kept[[idCol]])) - lst$nIndividuals
  row$tierAObsDelta <- if (isTRUE(is.na(obsPass))) NA_integer_ else obs$n - lst$nObsRecs
  row$tierAObsBasis <- obs$basis

  ## A filter that removed nothing passed a test that proved nothing.
  frac <- if (nrow(rawDf) > 0) (nrow(rawDf) - lst$nDataRecs) / nrow(rawDf) else NA_real_
  row$tierAInformative <- !is.na(frac) && frac > 0
  row$tierAFilterFraction <- if (is.na(frac)) NA_character_ else
    cut(frac, c(-Inf, 0, .01, .1, .5, Inf), labels = c("0", "<1%", "1-10%", "10-50%", ">50%"))[1]

  checks <- c(recPass, idPass, obsPass)
  row$tierAStatus <- if (any(checks %in% FALSE)) "FAIL" else "PASS"
  if (row$tierAStatus == "FAIL") {
    row$tierAFailKind <- if (all(checks %in% FALSE)) "ALL_THREE" else
      paste(c("RECS", "IDS", "OBS")[which(checks %in% FALSE)], collapse = "+")
  }
  row
}


# -- refusal taxonomy ---------------------------------------------------------
#
# Three signals, stored separately and never collapsed: the message match, the
# NULL condition call (deliberate stops use call. = FALSE), and whether the
# message names the file and line. Keeping all three means the taxonomy can be
# revised later from returned CSVs without asking anyone to re-run.

PMXF_REFUSALS <- list(
  UNSUPPORTED_PK     = "Unsupported NONMEM construct",
  NO_PK              = "^No \\$PK record found in ",
  EMPTY_PK           = "contains no statements",
  READ_BEFORE_ASSIGN = "before assigning it|before every path has assigned",
  UNKNOWN_SYMBOL     = "never assigns it and it is not in \\$INPUT|does not assign",
  NO_COV_REF         = "^No reference value could be derived from ",
  THETA_COUNT        = "THETA\\(s\\) but \\$PK references|THETA index",
  NOT_ASSIGNED       = "^Not assigned in the \\$PK block of ",
  COVREF_BAD         = "covRef",
  DATA_TEXT_VALUE    = "compares against a text value",
  DATA_BOTH_LISTS    = "cannot both appear in the \\$DATA record",
  DATA_MISSING_COL   = "which \\$INPUT does not declare",
  DATA_CHAR_COL      = "read as text rather than numbers",
  NO_INPUT           = "^No \\$INPUT record found in ",
  FILE_MISSING       = "^Model file not found: "
)

## Not a refusal despite reading like one: an internal deparse failure is a bug.
PMXF_CRASHES <- c(
  INTERNAL   = "^Internal error",
  BASE_CRASH = paste(
    "subscript out of bounds", "missing value where TRUE/FALSE",
    "object '.*' not found", "non-numeric argument", "undefined columns",
    "argument .* of length zero", "invalid 'times'",
    sep = "|"
  ),
  TIMEOUT = "reached (elapsed|CPU) time limit"
)

pmxfClassifyRefusal <- function(cond, modFile) {
  msg <- conditionMessage(cond)
  callNull <- is.null(conditionCall(cond))
  base <- basename(modFile)
  namesFile <- grepl(base, msg, fixed = TRUE)
  line <- suppressWarnings(as.integer(sub(paste0(".*", base, "[: ]+(?:line )?(\\d+).*"), "\\1", msg)))
  if (is.na(line) || !grepl("\\d", msg)) line <- NA_integer_

  for (nm in names(PMXF_CRASHES)) {
    if (grepl(PMXF_CRASHES[[nm]], msg)) {
      return(list(kind = nm, deliberate = FALSE, callNull = callNull,
                  namesFile = namesFile, line = line, msg = msg))
    }
  }
  for (nm in names(PMXF_REFUSALS)) {
    if (grepl(PMXF_REFUSALS[[nm]], msg)) {
      return(list(kind = nm, deliberate = TRUE, callNull = callNull,
                  namesFile = namesFile, line = line, msg = msg))
    }
  }
  ## Unlisted. If it looks deliberate on the other two signals, say so and flag
  ## it as a taxonomy gap rather than reporting a spurious failure.
  if (callNull && namesFile) {
    list(kind = "LIKELY_REFUSAL_UNLISTED", deliberate = TRUE, callNull = callNull,
         namesFile = namesFile, line = line, msg = msg)
  } else {
    list(kind = "UNCLASSIFIED", deliberate = FALSE, callNull = callNull,
         namesFile = namesFile, line = line, msg = msg)
  }
}


# -- Tier B: createParamFunction, structural ----------------------------------

pmxfThetasFor <- function(ctx, out) {
  n <- out$noBaseThetas
  if (!is.na(ctx$art$ext)) {
    e <- pmxfSafely(PMXForest::getExt(ctx$art$ext))
    if (e$ok) {
      fin <- e$value[e$value$ITERATION == -1000000000, , drop = FALSE]
      if (nrow(fin) && ncol(fin) >= n + 1L) {
        return(list(thetas = as.numeric(fin[1, 2:(n + 1L)]), source = "EXT"))
      }
    }
  }
  list(thetas = rep(0.5, n), source = "SYNTHETIC")
}

pmxfProbe <- function(f, thetas, df) {
  r <- pmxfSafely(f(thetas = thetas, df = df))
  if (!r$ok) return(list(ok = FALSE, finite = NA, msg = conditionMessage(r$cond)))
  v <- r$value
  list(ok = TRUE,
       finite = all(vapply(v, function(z) length(z) == 1L && is.numeric(z) && is.finite(z), NA)),
       value = v)
}

pmxfTierB <- function(ctx, row) {
  cpf <- pmxfSafely(PMXForest::createParamFunction(ctx$modFile, quiet = TRUE))
  row$tierBSeconds <- cpf$seconds
  row$tierBCovRefInferred <- sum(grepl("inferred rather than read", cpf$warnings))

  if (!cpf$ok) {
    cls <- pmxfClassifyRefusal(cpf$cond, ctx$modFile)
    row$tierBStatus <- if (cls$deliberate) "PASS_REFUSED" else "ERROR"
    row$tierBRefusalKind <- cls$kind
    row$tierBDeliberate <- cls$deliberate
    row$tierBCallNull <- cls$callNull
    row$tierBNamesFile <- cls$namesFile
    row$tierBErrLine <- cls$line
    row$tierBErrSig <- pmxfSanitise(cls$msg, ctx)
    if (cls$kind %in% c("UNCLASSIFIED", "LIKELY_REFUSAL_UNLISTED")) {
      row$flag <- TRUE
      row$flagReason <- paste(na.omit(c(row$flagReason, paste0("tierB:", cls$kind))), collapse = ";")
    }
    return(list(row = row, out = NULL))
  }

  out <- cpf$value
  row$nCov <- length(out$covRef)
  row$nParam <- length(out$primaryNames)
  row$nSecondary <- length(out$secondaryNames)
  row$nTheta <- out$noBaseThetas

  row$tierBParses <- isTRUE(pmxfSafely(parse(text = paste(out$code, collapse = "\n")))$ok)
  fEnv <- new.env(parent = globalenv())
  ev <- pmxfSafely(eval(parse(text = paste(out$code, collapse = "\n")), envir = fEnv))
  row$tierBEvals <- ev$ok && is.function(ev$value)
  if (!isTRUE(row$tierBEvals)) {
    row$tierBStatus <- "FAIL"
    return(list(row = row, out = out))
  }
  f <- ev$value
  row$tierBFormalsOk <- identical(names(formals(f)), c("thetas", "df", "..."))

  ## Free structural oracles from the .lst and the .ext, available on every
  ## model with artefacts - which matters because Tier C almost never runs.
  row$tierBThetaVsLst <- if (is.na(ctx$lst$lengthTheta)) NA_character_ else
    if (ctx$lst$lengthTheta == out$noBaseThetas) "SAME" else "DIFF"

  th <- pmxfThetasFor(ctx, out)
  row$tierBThetaSource <- th$source

  covs <- names(out$covRef)
  mkRow <- function(vals) {
    if (!length(covs)) return(data.frame(dummy = 1))
    as.data.frame(stats::setNames(as.list(vals), covs))
  }
  refRow <- mkRow(vapply(out$covRef, function(z) as.numeric(z$value), 0))
  missRow <- mkRow(rep(out$missVal, length(covs)))

  p <- pmxfProbe(f, th$thetas, refRow)
  row$tierBRefFinite <- isTRUE(p$finite)
  row$tierBReturnShapeOk <- p$ok && isTRUE(setequal(names(p$value), out$functionListName))
  row$tierBMissValOk <- isTRUE(pmxfProbe(f, th$thetas, missRow)$ok)
  row$tierBAbsentColOk <- isTRUE(pmxfProbe(f, th$thetas, data.frame(dummy = 1))$ok)
  ## These two are EXPECTED to error given the emitted one-row guard. Recorded
  ## as information, not failure - but aggregated fleet-wide they say whether
  ## NA-hardening would be worth doing.
  if (length(covs)) {
    naRow <- refRow
    naRow[[1]] <- NA_real_
    row$tierBNaSafe <- isTRUE(pmxfProbe(f, th$thetas, naRow)$ok)
    row$tierBMultiRowSafe <- isTRUE(pmxfProbe(f, th$thetas, rbind(refRow, refRow))$ok)
  }

  ok <- c(row$tierBParses, row$tierBEvals, row$tierBFormalsOk, row$tierBReturnShapeOk,
          row$tierBMissValOk, row$tierBAbsentColOk)
  synthetic <- identical(th$source, "SYNTHETIC")
  if (!synthetic) ok <- c(ok, row$tierBRefFinite)
  if (identical(row$tierBThetaVsLst, "DIFF")) ok <- c(ok, FALSE)
  row$tierBStatus <- if (all(ok %in% TRUE)) "PASS" else "FAIL"
  list(row = row, out = out)
}


# -- Tier C: verifyParamFunction ----------------------------------------------

pmxfTableRecords <- function(recs, omegaDim) {
  txts <- pmxfRecordAll(recs, "TABLE")
  if (!length(txts)) return(list())
  opts <- c("NOPRINT", "PRINT", "ONEHEADER", "NOHEADER", "NOTITLE", "NOLABEL",
            "FIRSTONLY", "LASTONLY", "FIRSTLASTONLY", "NOAPPEND", "APPEND",
            "UNCONDITIONAL", "CONDITIONAL", "OMITTED", "FORMAT", "RFORMAT")
  lapply(txts, function(t) {
    file <- if (grepl("FILE\\s*=", t, ignore.case = TRUE)) {
      gsub("^['\"]|['\"]$", "", sub(".*FILE\\s*=\\s*([^ \t]+).*", "\\1", t, ignore.case = TRUE))
    } else NA_character_
    toks <- strsplit(trimws(sub("^\\s*\\$[A-Za-z]+", "", t)), "\\s+")[[1]]
    toks <- toks[nzchar(toks) & !grepl("=", toks) & !toupper(toks) %in% opts]
    cols <- unlist(lapply(toks, function(tk) {
      m <- regmatches(tk, regexec("^ETAS?\\((\\d+):(\\d+|LAST)\\)$", toupper(tk)))[[1]]
      if (length(m) == 3) {
        hi <- if (m[3] == "LAST") omegaDim else as.integer(m[3])
        if (is.na(hi)) return(character(0))
        paste0("ETA", seq.int(as.integer(m[2]), hi))
      } else tk
    }))
    noAppend <- grepl("NOAPPEND", t, ignore.case = TRUE)
    list(file = file, cols = c(cols, if (!noAppend) c("DV", "PRED", "RES", "WRES")),
         firstOnly = grepl("FIRSTONLY", t, ignore.case = TRUE))
  })
}

pmxfSniffTable <- function(tabFile) {
  L <- pmxfReadLines(tabFile, n = 3L)
  if (length(L) < 2) return(list(ok = FALSE))
  banner <- grepl("^\\s*TABLE NO", L[1])
  hdr <- if (banner) L[2] else L[1]
  nComma <- length(strsplit(hdr, ",")[[1]])
  nWs <- length(strsplit(trimws(hdr), "[ \t]+")[[1]])
  sep <- if (nComma > 1 && nComma >= nWs) "," else ""
  nmz <- strsplit(trimws(hdr), if (nzchar(sep)) sep else "[ \t]+")[[1]]
  fmt <- paste0(if (banner) "BANNER" else "NOBANNER", "_", if (nzchar(sep)) "CSV" else "WS")
  list(ok = TRUE, banner = banner, sep = sep, names = make.names(trimws(nmz)), format = fmt)
}

pmxfTierCFeasible <- function(out, cols, ext) {
  if (is.na(ext)) return(list(feasible = FALSE, reason = "NO_EXT"))
  covs <- names(out$covRef)
  if (length(covs) && !all(covs %in% cols)) {
    return(list(feasible = FALSE, reason = "COVS_NOT_TABLED"))
  }
  can <- vapply(out$primaryNames, function(p) {
    if (paste0("TV", p) %in% cols) return(TRUE)
    p %in% cols && p %in% names(out$etaMap) && paste0("ETA", out$etaMap[[p]]) %in% cols
  }, NA)
  if (!any(can)) {
    return(list(feasible = FALSE,
                reason = if (!any(out$primaryNames %in% cols)) "PARAMS_NOT_TABLED" else "ETA_NOT_TABLED"))
  }
  list(feasible = TRUE, params = out$primaryNames[can])
}

## verifyParamFunction() assumes a "TABLE NO." banner and whitespace separation.
## Real tables are not always that, so normalise a copy in tempdir(). The shim
## being needed is itself reported - it is a finding about the package.
pmxfTableShim <- function(tabFile, sniff) {
  if (identical(sniff$format, "BANNER_WS")) return(list(path = tabFile, shim = "NONE"))
  L <- pmxfReadLines(tabFile, n = PMXF_MAX_TABLE_ROWS)
  body <- if (sniff$banner) L[-1] else L
  if (nzchar(sniff$sep)) body <- gsub(",", " ", body, fixed = TRUE)
  p <- tempfile(fileext = ".tab")
  writeLines(c("TABLE NO.  1", body), p)
  list(path = p, shim = if (nzchar(sniff$sep)) "CSV_CONVERTED" else "BANNER_ADDED")
}

pmxfTierC <- function(ctx, out, row) {
  skip <- function(reason) {
    row$tierCStatus <- "SKIP"
    row$tierCSkipReason <- reason
    row
  }
  if (is.null(out)) return(skip("TIER_B_FAILED"))
  tabs <- pmxfTableRecords(ctx$recs, ctx$lst$omegaDim)
  row$nTableRecords <- length(tabs)
  if (!length(tabs)) return(skip("NO_TABLE_RECORD"))

  present <- Filter(function(t) !is.na(t$file) && file.exists(file.path(ctx$art$dir, t$file)), tabs)
  row$nTableFilesPresent <- length(present)
  if (!length(present)) return(skip("TABLE_FILE_MISSING"))

  best <- NULL
  for (t in present) {
    path <- file.path(ctx$art$dir, t$file)
    sn <- pmxfSniffTable(path)
    if (!isTRUE(sn$ok)) next
    fe <- pmxfTierCFeasible(out, sn$names, ctx$art$ext)
    if (isTRUE(fe$feasible)) {
      best <- list(path = path, sniff = sn, fe = fe)
      break
    }
    if (is.null(best)) best <- list(path = path, sniff = sn, fe = fe)
  }
  if (is.null(best)) return(skip("TABLE_UNREADABLE"))
  row$tabFormat <- best$sniff$format
  row$tabHasParams <- any(out$primaryNames %in% best$sniff$names)
  row$tabHasCovs <- all(names(out$covRef) %in% best$sniff$names)

  if (!isTRUE(best$fe$feasible)) {
    row$tierCSuggestionAvailable <- TRUE
    ## Only parameters you would actually plot: skip the model's own TV* and
    ## MU_* intermediates, which createParamFunction() returns but nobody
    ## forest-plots, and which would otherwise make the suggestion unusable.
    want <- out$primaryNames[!grepl("^(TV|MU_|COV[0-9]*$)", out$primaryNames)]
    if (!length(want)) want <- out$primaryNames
    want <- utils::head(want, 8)
    ctx$suggest(sprintf(
      "$TABLE ID %s %s ONEHEADER NOPRINT FILE=pmxfverify.tab",
      paste(names(out$covRef), collapse = " "),
      paste0("TV", want, collapse = " ")
    ))
    return(skip(best$fe$reason))
  }

  sh <- pmxfTableShim(best$path, best$sniff)
  row$tierCShim <- sh$shim
  th <- pmxfThetasFor(ctx, out)
  if (identical(th$source, "SYNTHETIC")) return(skip("NO_FINAL_ESTIMATES"))

  v <- pmxfSafely(PMXForest::verifyParamFunction(out, sh$path, th$thetas, quiet = TRUE))
  row$tierCSeconds <- v$seconds
  if (!v$ok) {
    row$tierCStatus <- "ERROR"
    row$tierCSkipReason <- pmxfSanitise(conditionMessage(v$cond), ctx)
    return(row)
  }
  ck <- attr(v$value, "checks")
  row$tierCnParams <- nrow(ck)
  row$tierCnPass <- sum(ck$PASS %in% TRUE)
  row$tierCnFail <- sum(ck$PASS %in% FALSE)
  row$tierCnNA <- sum(is.na(ck$PASS))
  row$tierCMaxRelDiff <- suppressWarnings(max(ck$MAXRELDIFF, na.rm = TRUE))
  row$tierCStatus <- if (row$tierCnFail > 0) "FAIL" else if (row$tierCnNA > 0) "PARTIAL" else "PASS"
  row
}


# -- Tier D: differential against a hand-written paramFunction ----------------

pmxfFindHandwritten <- function(modFile, handwritten, manifest) {
  base <- sub("\\.[^.]*$", "", basename(modFile))
  if (!is.null(manifest) && nrow(manifest)) {
    hit <- manifest[normalizePath(manifest$modFile, mustWork = FALSE) ==
                      normalizePath(modFile, mustWork = FALSE), , drop = FALSE]
    if (nrow(hit)) return(list(file = hit$funFile[1], funName = hit$funName[1] %||% NA, source = "MANIFEST"))
  }
  dirs <- unique(c(dirname(modFile), handwritten))
  cand <- as.vector(outer(dirs, c(
    paste0(base, ".paramFunction.R"), paste0("paramFunction_", base, ".R"), "paramFunction.R"
  ), file.path))
  cand <- cand[file.exists(cand)]
  if (length(cand)) list(file = cand[1], funName = NA, source = "SIDECAR") else NULL
}

pmxfTierD <- function(ctx, out, row) {
  skip <- function(r) {
    row$tierDStatus <- "SKIP"
    row$tierDSource <- r
    row
  }
  if (is.null(out)) return(skip("TIER_B_FAILED"))
  hw <- pmxfFindHandwritten(ctx$modFile, ctx$handwritten, ctx$manifest)
  if (is.null(hw)) return(skip("NO_HANDWRITTEN"))

  env <- new.env(parent = globalenv())
  ld <- pmxfSafely(sys.source(hw$file, envir = env))
  if (!ld$ok) return(skip("LOAD_FAILED"))
  fns <- Filter(function(n) is.function(get(n, envir = env)), ls(env))
  nm <- if (!is.na(hw$funName) && hw$funName %in% fns) hw$funName else
    if ("paramFunction" %in% fns) "paramFunction" else if (length(fns) == 1L) fns[1] else NA
  if (is.na(nm)) return(skip("AMBIGUOUS_FUNCTION"))
  hand <- get(nm, envir = env)

  fmls <- names(formals(hand))
  conv <- if (identical(fmls[1:2], c("thetas", "df"))) "NEW" else
    if (any(c("basethetas", "dfrow") %in% fmls)) "LEGACY" else "UNKNOWN"
  row$tierDConvention <- conv
  row$tierDConventionRisk <- conv == "LEGACY"
  if (conv == "UNKNOWN") return(skip("SIGNATURE"))

  gen <- eval(parse(text = paste(out$code, collapse = "\n")), envir = new.env(parent = globalenv()))
  th <- pmxfThetasFor(ctx, out)
  covs <- names(out$covRef)
  refRow <- if (length(covs)) {
    as.data.frame(stats::setNames(as.list(vapply(out$covRef, function(z) as.numeric(z$value), 0)), covs))
  } else data.frame(dummy = 1)

  grid <- list(refRow)
  for (cv in covs) {
    for (mult in c(0.5, 0.75, 1.25, 2)) {
      g <- refRow
      g[[cv]] <- g[[cv]] * mult
      grid[[length(grid) + 1L]] <- g
    }
    g <- refRow
    g[[cv]] <- out$missVal
    grid[[length(grid) + 1L]] <- g
  }
  row$tierDGridRows <- length(grid)

  callHand <- function(dfRow) {
    if (conv == "NEW") hand(thetas = th$thetas, df = dfRow) else
      hand(basethetas = th$thetas, covthetas = numeric(0), dfrow = dfRow, etas = rep(0, 10))
  }
  nDiff <- 0L
  maxRel <- 0
  firstDiff <- NA_integer_
  common <- NULL
  for (i in seq_along(grid)) {
    a <- pmxfSafely(gen(thetas = th$thetas, df = grid[[i]]))
    b <- pmxfSafely(callHand(grid[[i]]))
    if (!a$ok || !b$ok) next
    av <- unlist(a$value)
    bv <- unlist(b$value)
    common <- intersect(names(av), names(bv))
    if (!length(common)) {
      row$tierDStatus <- "SHAPE"
      return(row)
    }
    rel <- abs(av[common] - bv[common]) / pmax(abs(bv[common]), 1e-12)
    rel[is.na(rel) & is.na(av[common]) & is.na(bv[common])] <- 0
    if (any(rel > 1e-8, na.rm = TRUE)) {
      nDiff <- nDiff + 1L
      if (is.na(firstDiff)) firstDiff <- i
    }
    maxRel <- max(maxRel, max(rel, na.rm = TRUE))
  }
  row$tierDCommonParams <- length(common %||% character(0))
  row$tierDnDiffer <- nDiff
  row$tierDMaxRelDiff <- maxRel
  row$tierDFirstDiffRow <- firstDiff
  row$tierDStatus <- if (nDiff > 0) "DIFF" else "PASS"
  row
}


# -- redaction ----------------------------------------------------------------

## The package quotes the offending source line in some errors, and control
## stream lines contain covariate names. Subject counts and covariate names are
## the two things that identify a study to anyone who knows the programme.
pmxfSanitise <- function(msg, ctx) {
  if (is.null(msg) || !length(msg)) return(NA_character_)
  msg <- paste(msg, collapse = " | ")
  msg <- gsub("\\s+", " ", msg)
  if (!ctx$full) {
    msg <- gsub(basename(ctx$modFile), "<model>", msg, fixed = TRUE)
    msg <- gsub("(/[^ '\"]+)+", "<path>", msg)
    msg <- gsub("'[^']+'", "'<sym>'", msg)
  }
  substr(msg, 1, 300)
}

PMXF_SENSITIVE <- c(
  "modelFile", "covNames", "paramNames", "tierARawRows", "tierAIdDelta",
  "tierAObsDelta", "tierARecDelta"
)

pmxfApplyRedaction <- function(row, full) {
  if (full) return(row)
  for (nm in intersect(PMXF_SENSITIVE, names(row))) row[[nm]] <- NA
  ## Deltas are reduced to their sign: zero carries the whole test result.
  for (nm in c("tierARecDelta", "tierAIdDelta", "tierAObsDelta")) row[[nm]] <- NA
  row
}

## Belt and braces. A new error message in a future release can leak a path at
## any time and no schema review would catch it.
pmxfOutboundAudit <- function(df) {
  n <- 0L
  for (j in seq_along(df)) {
    if (!is.character(df[[j]])) next
    bad <- !is.na(df[[j]]) &
      (grepl("^/", df[[j]]) | grepl("^[A-Za-z]:\\\\", df[[j]]) |
         grepl("@", df[[j]]) | nchar(df[[j]]) > 300)
    if (any(bad)) {
      df[[j]][bad] <- "<redacted>"
      n <- n + sum(bad)
    }
  }
  list(df = df, blanked = n)
}


# -- the report row -----------------------------------------------------------

pmxfNewRow <- function() {
  as.data.frame(c(
    list(schemaVersion = 1L, runId = NA_character_, runStamp = NA_character_,
         scriptVersion = PMXF_SCRIPT_VERSION, pmxforestVersion = NA_character_,
         pmxforestPriorVersion = NA_character_, rVersion = NA_character_,
         platform = NA_character_, redactLevel = NA_character_,
         modelId = NA_character_, pathHash = NA_character_, modelFile = NA_character_,
         pkHash = NA_character_, dupOfModelId = NA_character_,
         modEol = NA_character_, modLines = NA_integer_, advan = NA_character_,
         nInput = NA_integer_, nTheta = NA_integer_, nCov = NA_integer_,
         nParam = NA_integer_, nSecondary = NA_integer_, covNames = NA_character_,
         isMuReferenced = NA, isSimulation = NA, hasLst = NA, hasExt = NA,
         lstEchoMatch = NA_character_, nTableRecords = NA_integer_,
         nTableFilesPresent = NA_integer_, tabFormat = NA_character_,
         tabHasParams = NA, tabHasCovs = NA,
         tierAStatus = NA_character_, tierASkipReason = NA_character_,
         tierARecPass = NA, tierAIdPass = NA, tierAObsPass = NA,
         tierARecDelta = NA_integer_, tierAIdDelta = NA_integer_,
         tierAObsDelta = NA_integer_, tierARawRows = NA_integer_,
         tierAObsBasis = NA_character_, tierAInformative = NA,
         tierAFilterFraction = NA_character_, tierAFailKind = NA_character_,
         tierAErrorKind = NA_character_,
         tierBStatus = NA_character_, tierBRefusalKind = NA_character_,
         tierBDeliberate = NA, tierBCallNull = NA, tierBNamesFile = NA,
         tierBErrLine = NA_integer_, tierBErrSig = NA_character_,
         tierBParses = NA, tierBEvals = NA, tierBFormalsOk = NA,
         tierBReturnShapeOk = NA, tierBRefFinite = NA, tierBMissValOk = NA,
         tierBAbsentColOk = NA, tierBNaSafe = NA, tierBMultiRowSafe = NA,
         tierBThetaSource = NA_character_, tierBThetaVsLst = NA_character_,
         tierBCovRefInferred = NA_integer_, tierBSeconds = NA_real_,
         tierCStatus = NA_character_, tierCSkipReason = NA_character_,
         tierCShim = NA_character_, tierCnParams = NA_integer_,
         tierCnPass = NA_integer_, tierCnFail = NA_integer_, tierCnNA = NA_integer_,
         tierCMaxRelDiff = NA_real_, tierCSuggestionAvailable = NA,
         tierCSeconds = NA_real_,
         tierDStatus = NA_character_, tierDSource = NA_character_,
         tierDConvention = NA_character_, tierDConventionRisk = NA,
         tierDGridRows = NA_integer_, tierDCommonParams = NA_integer_,
         tierDnDiffer = NA_integer_, tierDMaxRelDiff = NA_real_,
         tierDFirstDiffRow = NA_integer_,
         overall = NA_character_, flag = FALSE, flagReason = NA_character_,
         totalSeconds = NA_real_)
  ), stringsAsFactors = FALSE)
}


# -- one model ----------------------------------------------------------------

pmxfOneModel <- function(modFile, idx, opts, prov) {
  t0 <- proc.time()[["elapsed"]]
  row <- pmxfNewRow()
  for (nm in names(prov)) row[[nm]] <- prov[[nm]]
  row$modelId <- sprintf("M%03d", idx)
  row$pathHash <- pmxfHash(normalizePath(modFile, mustWork = FALSE), opts$salt)
  row$modelFile <- basename(modFile)

  lines <- pmxfReadLines(modFile)
  row$modEol <- pmxfEolKind(modFile)
  row$modLines <- length(lines)
  clean <- pmxfStripComments(lines)
  recs <- pmxfRecords(clean)
  pkTxt <- pmxfRecordAll(recs, "PK")
  row$pkHash <- if (length(pkTxt)) pmxfHash(pmxfNorm(pkTxt)) else NA_character_
  row$isMuReferenced <- any(grepl("\\bMU_[0-9]+", clean))
  row$isSimulation <- any(recs$record == "SIMULATION")
  row$advan <- {
    s <- pmxfRecordAll(recs, "SUBROUTINE")
    if (length(s)) toupper(sub(".*?(ADVAN\\s*=?\\s*[0-9]+).*", "\\1", s[1])) else NA_character_
  }
  inputNames <- pmxfInputNames(recs)
  row$nInput <- length(inputNames)

  art <- pmxfFindArtefacts(modFile)
  row$hasLst <- !is.na(art$lst)
  row$hasExt <- !is.na(art$ext)
  lst <- if (!is.na(art$lst)) pmxfParseLst(art$lst) else pmxfEmptyLst()
  echo <- if (length(lst$echo)) pmxfEchoAgreement(lines, lst$echo) else list(verdict = "NO_ECHO")
  row$lstEchoMatch <- echo$verdict

  ctx <- list(modFile = modFile, recs = recs, inputNames = inputNames, art = art,
              lst = lst, echo = echo, dataSpec = pmxfDataSpec(recs),
              full = opts$full, handwritten = opts$handwritten,
              manifest = opts$manifest, suggest = opts$suggest)

  row <- pmxfTierA(ctx, row)
  tb <- pmxfTierB(ctx, row)
  row <- tb$row
  row <- pmxfTierC(ctx, tb$out, row)
  if (opts$tierD) row <- pmxfTierD(ctx, tb$out, row)
  if (!is.null(tb$out)) row$covNames <- paste(names(tb$out$covRef), collapse = " ")

  row$overall <- if (identical(row$tierBStatus, "ERROR") || identical(row$tierAStatus, "ERROR")) {
    "ERROR"
  } else if (any(c(row$tierAStatus, row$tierBStatus, row$tierCStatus, row$tierDStatus) %in%
                 c("FAIL", "DIFF"))) {
    "FAIL"
  } else if (identical(row$tierBStatus, "PASS_REFUSED")) {
    "REFUSED"
  } else {
    "PASS"
  }
  if (row$overall %in% c("ERROR", "FAIL")) row$flag <- TRUE
  row$totalSeconds <- round(proc.time()[["elapsed"]] - t0, 2)
  row
}


# -- self-test ----------------------------------------------------------------

## The script's own parsers are the newest and least-tested code in the
## pipeline. Before judging anyone's models, it must reproduce a result that is
## known to be right. If it cannot, the harness is wrong, not the corpus.
pmxfSelfTest <- function() {
  d <- system.file("extdata", "SimVal", package = "PMXForest")
  if (!nzchar(d)) return(list(ok = FALSE, why = "bundled SimVal data not found"))
  mod <- file.path(d, "run7.mod")
  lst <- file.path(d, "run7.lst")
  dat <- file.path(d, "DAT-1-MI-PMX-2.csv")
  if (!all(file.exists(mod, lst, dat))) return(list(ok = FALSE, why = "SimVal incomplete"))
  p <- pmxfParseLst(lst)
  raw <- utils::read.csv(dat)
  kept <- PMXForest::filterByModel(raw, mod, quiet = TRUE)
  obs <- pmxfObsCount(kept)
  got <- c(nrow(kept), length(unique(kept$ID)), obs$n)
  want <- c(p$nDataRecs, p$nIndividuals, p$nObsRecs)
  ok <- identical(as.integer(got), as.integer(want)) &&
    identical(as.integer(want), c(33885L, 754L, 6284L))
  list(ok = ok, got = got, want = want,
       why = if (ok) "" else sprintf("expected 33885/754/6284, .lst gave %s, filterByModel gave %s",
                                     paste(want, collapse = "/"), paste(got, collapse = "/")))
}


# -- main ---------------------------------------------------------------------

pmxfDiscoverModels <- function(paths) {
  out <- character(0)
  for (p in paths) {
    p <- path.expand(p)
    if (dir.exists(p)) {
      out <- c(out, list.files(p, pattern = "\\.(mod|ctl|con)$", recursive = TRUE,
                               full.names = TRUE, ignore.case = TRUE))
    } else if (file.exists(p)) {
      out <- c(out, p)
    } else {
      warning("path not found, skipped: ", p, call. = FALSE)
    }
  }
  sort(unique(out))
}

verifyPMXForest <- function(paths,
                            out = NULL,
                            full = FALSE,
                            install = FALSE,
                            handwritten = NULL,
                            manifest = NULL,
                            tierD = TRUE,
                            maxModels = Inf) {
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
  outDir <- out %||% file.path(getwd(), paste0("pmxf-stress-", stamp))
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

  pmxfRule("PMXForest stress test")
  priorVer <- tryCatch(as.character(utils::packageVersion("PMXForest")), error = function(e) NA_character_)
  priorLib <- tryCatch(dirname(find.package("PMXForest")), error = function(e) NA_character_)
  pmxfMsg("Installed PMXForest : ", priorVer %||% "none", "  in  ", priorLib %||% "-")
  pmxfMsg("Wanted for this test: ", PMXF_WANT_VERSION)

  if (install) {
    pmxfMsg("")
    pmxfMsg("Installing will OVERWRITE PMXForest ", priorVer, " in ", priorLib, ".")
    pmxfMsg("That is a change on disk, not just in this session.")
    if (!requireNamespace("PMXRenv", quietly = TRUE)) {
      stop("PMXRenv is not available; install PMXForest ", PMXF_WANT_VERSION, " yourself.", call. = FALSE)
    }
    PMXRenv::activate.unqualified.packages()
    PMXRenv::install.unqualified.packages("PMXForest", repoName = "development")
    pmxfMsg("To go back: reinstall PMXForest ", priorVer, " into ", priorLib,
            ", then start a FRESH R session (activate.unqualified.packages() ",
            "changes .libPaths() for this session only).")
  }

  have <- as.character(utils::packageVersion("PMXForest"))
  if (!identical(have, PMXF_WANT_VERSION)) {
    pmxfMsg("")
    pmxfMsg("PMXForest ", have, " is loaded but this test wants ", PMXF_WANT_VERSION, ".")
    pmxfMsg("Run these two lines, then start a fresh R session and try again:")
    pmxfMsg('  PMXRenv::activate.unqualified.packages()')
    pmxfMsg('  PMXRenv::install.unqualified.packages("PMXForest", repoName = "development")')
    stop("wrong PMXForest version", call. = FALSE)
  }

  pmxfRule("Self-test")
  st <- pmxfSelfTest()
  if (!st$ok) {
    pmxfMsg("FAILED: ", st$why)
    stop("the harness itself is wrong; not scanning your models", call. = FALSE)
  }
  pmxfMsg("OK - reproduced 33885 records / 754 individuals / 6284 observations ",
          "on the bundled example.")

  models <- pmxfDiscoverModels(paths)
  if (length(models) > maxModels) models <- models[seq_len(maxModels)]
  pmxfRule("Corpus")
  pmxfMsg("Control streams found: ", length(models))
  if (!length(models)) return(invisible(pmxfNewRow()[0, ]))

  if (is.character(manifest)) manifest <- utils::read.csv(manifest, stringsAsFactors = FALSE)
  suggestions <- new.env(parent = emptyenv())
  suggestions$v <- character(0)

  opts <- list(full = full, salt = pmxfHash(stamp), handwritten = handwritten,
               manifest = manifest, tierD = tierD,
               suggest = function(s) suggestions$v <- c(suggestions$v, s))
  prov <- list(runId = stamp, runStamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
               pmxforestVersion = have, pmxforestPriorVersion = priorVer,
               rVersion = as.character(getRversion()), platform = R.version$platform,
               redactLevel = if (full) "full" else "safe")

  csvPath <- file.path(outDir, paste0("pmxf-stress-", stamp, "-report.csv"))
  rows <- vector("list", length(models))
  for (i in seq_along(models)) {
    cat(sprintf("\r[%d/%d] %-50s", i, length(models),
                substr(basename(models[i]), 1, 50)))
    r <- tryCatch(pmxfOneModel(models[i], i, opts, prov),
      error = function(e) {
        rr <- pmxfNewRow()
        for (nm in names(prov)) rr[[nm]] <- prov[[nm]]
        rr$modelId <- sprintf("M%03d", i)
        rr$overall <- "ERROR"
        rr$flag <- TRUE
        rr$flagReason <- "harness"
        rr$tierBErrSig <- substr(conditionMessage(e), 1, 300)
        rr
      }
    )
    rows[[i]] <- r
    ## Journaled: a hard crash loses one model, not the run.
    utils::write.table(pmxfApplyRedaction(r, full), csvPath, sep = ",",
      row.names = FALSE, col.names = (i == 1L), append = (i > 1L),
      qmethod = "double", na = ""
    )
  }
  cat("\r", strrep(" ", 70), "\r", sep = "")
  report <- do.call(rbind, rows)

  ## Rewrite the CSV once through the outbound audit.
  red <- pmxfApplyRedaction(report, full)
  aud <- pmxfOutboundAudit(red)
  utils::write.csv(aud$df, csvPath, row.names = FALSE, na = "")
  rdsPath <- file.path(outDir, paste0("pmxf-stress-", stamp, "-local.rds"))
  saveRDS(list(report = report, suggestions = suggestions$v, models = models), rdsPath)

  pmxfSummarise(report, suggestions$v, csvPath, rdsPath, aud$blanked, full)
  invisible(report)
}

pmxfTally <- function(x) {
  t <- table(factor(x[!is.na(x)]))
  if (!length(t)) return("  (none)")
  paste0("  ", format(names(t), width = 22), " ", as.integer(t), collapse = "\n")
}

pmxfSummarise <- function(report, suggestions, csvPath, rdsPath, blanked, full) {
  nDistinct <- length(unique(report$pkHash[!is.na(report$pkHash)]))
  pmxfRule("Results")
  pmxfMsg("Models: ", nrow(report), "   distinct $PK blocks: ", nDistinct)
  for (tier in c("A", "B", "C", "D")) {
    col <- paste0("tier", tier, "Status")
    pmxfMsg("\nTier ", tier, ":")
    cat(pmxfTally(report[[col]]), "\n")
  }

  infA <- report[report$tierAInformative %in% TRUE, ]
  pmxfRule("Headline")
  pmxfMsg("Tier A pass rate (informative tests only): ",
          sum(infA$tierAStatus == "PASS"), "/", nrow(infA))
  pmxfMsg("Tier B unclassified errors               : ",
          sum(report$tierBRefusalKind %in% c("UNCLASSIFIED", "LIKELY_REFUSAL_UNLISTED")))
  pmxfMsg("Tier B crashes                           : ",
          sum(report$tierBStatus %in% "ERROR"))
  pmxfMsg("Tier C shimmed (table format unreadable) : ",
          sum(!report$tierCShim %in% c(NA, "NONE")))
  pmxfMsg("Tier C failures                          : ", sum(report$tierCnFail %in% TRUE))
  pmxfMsg("Tier D differences                       : ", sum(report$tierDStatus %in% "DIFF"))

  flagged <- report[report$flag %in% TRUE, ]
  if (nrow(flagged)) {
    pmxfRule("Needs a look")
    for (i in seq_len(nrow(flagged))) {
      pmxfMsg("  ", flagged$modelId[i], "  ", flagged$overall[i], "  ",
              flagged$tierBRefusalKind[i] %||% "", "  ", flagged$tierBErrSig[i] %||% "")
    }
  }

  if (length(suggestions)) {
    pmxfRule("To make more models verifiable next time you run them")
    pmxfMsg("Add a $TABLE like this (console only - not written to the shared file):")
    for (s in unique(suggestions)[seq_len(min(5, length(unique(suggestions))))]) pmxfMsg("  ", s)
  }

  pmxfRule("Files")
  pmxfMsg("SEND BACK : ", csvPath)
  pmxfMsg("KEEP LOCAL: ", rdsPath, "  (paths, covariate names, raw messages)")
  if (blanked) pmxfMsg(blanked, " cell(s) blanked by the outbound filter.")
  if (!full) {
    pmxfMsg("\nThe file to send back is redacted: no paths, no covariate names, ",
            "no subject counts.")
    pmxfMsg("First rows of what you would be sending:")
    print(utils::head(utils::read.csv(csvPath)[, c("modelId", "pkHash", "advan",
      "tierAStatus", "tierBStatus", "tierCStatus", "overall")], 3))
  }
}


# -- command line -------------------------------------------------------------

if (!interactive() && identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) {
    flags <- grepl("^--", args)
    opt <- args[flags]
    verifyPMXForest(
      paths = args[!flags],
      full = "--full" %in% opt,
      install = "--install" %in% opt,
      tierD = !"--no-tier-d" %in% opt
    )
  }
}
