## checkFilterByModel.R ---------------------------------------------------------
##
## Run PMXForest::verifyFilterByModel() over a list of real models and write a
## single report you can paste back.
##
## For each model it compares what filterByModel() keeps against the record,
## subject and observation counts NONMEM printed in that run's .lst. Nothing is
## re-run; it reads the control stream, the .lst and the $DATA file.
##
##   Rscript checkFilterByModel.R
##   Rscript checkFilterByModel.R --out ~/pmxf-filter-check.txt
##   Rscript checkFilterByModel.R --max-mb 200      # skip larger data files
##
## Writes two files next to each other:
##   <out>.txt   a compact report, meant to be pasted whole
##   <out>.csv   the same thing per-model, if you would rather attach it
##
## Nothing is written anywhere near the run directories.
## ---------------------------------------------------------------------------- ##

## Model paths are NOT baked in. This file ships inside the package, so it must
## not carry anyone's project layout. Give the models on the command line, or
## put one path per line in a file and pass --models.
##
##   Rscript checkFilterByModel.R a/run1.mod b/run2.mod
##   Rscript checkFilterByModel.R --models mymodels.txt
MODELS <- character(0)

args <- commandArgs(trailingOnly = TRUE)
argOf <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
OUT <- argOf("--out", file.path(getwd(), "pmxf-filter-check"))
MAXMB <- as.numeric(argOf("--max-mb", "500"))

flagged <- c("--out", "--max-mb", "--models")
taken <- sort(unique(c(match(flagged, args), match(flagged, args) + 1L)))
loose <- args[setdiff(seq_along(args), taken[!is.na(taken)])]
MODELS <- c(MODELS, loose[!grepl("^--", loose)])
modelsFile <- argOf("--models", NA_character_)
if (!is.na(modelsFile) && file.exists(modelsFile)) {
  ln <- trimws(readLines(modelsFile, warn = FALSE))
  MODELS <- c(MODELS, ln[nzchar(ln) & !startsWith(ln, "#")])
}
if (length(MODELS) == 0) {
  cat("No models given.\n\n",
    "  Rscript checkFilterByModel.R path/to/run1.mod path/to/run2.mod\n",
    "  Rscript checkFilterByModel.R --models mymodels.txt\n\n",
    "Each model needs its .lst and its $DATA file beside it.\n",
    sep = ""
  )
  quit(status = 1)
}

suppressPackageStartupMessages(library(PMXForest))

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- one model ------------------------------------------------------------

checkOne <- function(mod) {
  r <- list(
    model = basename(mod), status = NA_character_, note = "",
    records = NA, subjects = NA, observations = NA,
    lstRecords = NA, lstSubjects = NA, lstObservations = NA,
    rawRows = NA, removed = NA, informative = NA, obsBasis = NA_character_,
    seconds = NA_real_, warnings = ""
  )
  t0 <- proc.time()[["elapsed"]]

  if (!file.exists(mod)) {
    r$status <- "NO_MODEL"
    return(r)
  }

  ## Size-gate the data file before opening it. Reading half of one would make
  ## the counts meaningless, so this skips rather than truncates.
  dat <- tryCatch(nmDataPath(mod), error = function(e) NA_character_)
  if (!is.na(dat) && file.exists(dat)) {
    mb <- file.size(dat) / 1e6
    r$note <- sprintf("data %.0f MB", mb)
    if (mb > MAXMB) {
      r$status <- "SKIP_DATA_TOO_LARGE"
      return(r)
    }
  }

  warns <- character(0)
  res <- withCallingHandlers(
    tryCatch(verifyFilterByModel(mod, quiet = TRUE), error = function(e) e),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  r$seconds <- round(proc.time()[["elapsed"]] - t0, 1)
  r$warnings <- paste(gsub("\\s+", " ", warns), collapse = " | ")

  if (inherits(res, "error")) {
    r$status <- "ERROR"
    r$note <- gsub("\\s+", " ", conditionMessage(res))
    return(r)
  }

  d <- attr(res, "checks")
  r$status <- if (isTRUE(unclass(res)[1])) "PASS" else "FAIL"
  r$records <- d$FILTERED[d$CHECK == "RECORDS"]
  r$subjects <- d$FILTERED[d$CHECK == "SUBJECTS"]
  r$observations <- d$FILTERED[d$CHECK == "OBSERVATIONS"]
  r$lstRecords <- d$NONMEM[d$CHECK == "RECORDS"]
  r$lstSubjects <- d$NONMEM[d$CHECK == "SUBJECTS"]
  r$lstObservations <- d$NONMEM[d$CHECK == "OBSERVATIONS"]
  r$removed <- attr(res, "removed")
  r$informative <- attr(res, "informative")
  r$obsBasis <- attr(res, "obsBasis")
  ## From the object, not reconstructed: filtered + removed only reconciles
  ## when the check passes, and a failing model is where this report matters.
  ## An older PMXForest has no such attribute, so say so rather than print NULL.
  r$rawRows <- attr(res, "rawRows") %||% NA_integer_
  r
}

## The $DATA path, resolved the way verifyFilterByModel() resolves it.
nmDataPath <- function(mod) {
  L <- sub("\r$", "", readLines(mod, warn = FALSE))
  i <- grep("^\\s*\\$DAT", L)
  if (!length(i)) {
    return(NA_character_)
  }
  j <- grep("^\\s*\\$", L)
  j <- j[j > i[1]]
  txt <- paste(L[i[1]:(if (length(j)) j[1] - 1L else length(L))], collapse = " ")
  raw <- strsplit(trimws(sub("^\\s*\\$[A-Za-z]+", "", txt)), "\\s+")[[1]][1]
  raw <- gsub("^['\"]|['\"]$", "", raw)
  cand <- c(raw, file.path(dirname(mod), raw),
    file.path(dirname(mod), basename(raw))
  )
  cand <- cand[file.exists(cand) & !dir.exists(cand)]
  if (length(cand)) normalizePath(cand[1]) else NA_character_
}

## ---- run ------------------------------------------------------------------

cat("PMXForest ", as.character(packageVersion("PMXForest")),
  "  |  ", R.version.string, "\n\n",
  sep = ""
)

rows <- list()
for (i in seq_along(MODELS)) {
  cat(sprintf("[%d/%d] %-18s ", i, length(MODELS), basename(MODELS[i])))
  utils::flush.console()
  rows[[i]] <- checkOne(MODELS[i])
  cat(rows[[i]]$status, "\n")
}
report <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))

## ---- the report you paste back --------------------------------------------

txt <- c(
  "=== verifyFilterByModel() over real models ===",
  paste0("PMXForest ", as.character(packageVersion("PMXForest")),
    " | ", R.version.string, " | ", format(Sys.time(), "%Y-%m-%d %H:%M")
  ),
  ""
)
for (r in rows) {
  txt <- c(txt, sprintf("%-16s %s%s", r$model, r$status,
    if (nzchar(r$note)) paste0("   [", r$note, "]") else ""
  ))
  if (r$status %in% c("PASS", "FAIL")) {
    txt <- c(txt, sprintf(
      "    %-13s filtered %10s   nonmem %10s   %s",
      c("RECORDS", "SUBJECTS", "OBSERVATIONS"),
      format(c(r$records, r$subjects, r$observations)),
      format(c(r$lstRecords, r$lstSubjects, r$lstObservations)),
      ifelse(c(r$records, r$subjects, r$observations) ==
        c(r$lstRecords, r$lstSubjects, r$lstObservations), "ok", "MISMATCH")
    ))
    txt <- c(txt, sprintf(
      "    raw rows %s, removed %s%s, observations counted from %s, %ss",
      format(r$rawRows), format(r$removed),
      if (isTRUE(r$informative)) "" else "  <- filter removed nothing",
      r$obsBasis, format(r$seconds)
    ))
  }
  if (nzchar(r$warnings)) {
    txt <- c(txt, paste0("    warning: ", substr(r$warnings, 1, 400)))
  }
  txt <- c(txt, "")
}
txt <- c(txt,
  "=== summary ===",
  paste0("  ", names(table(report$status)), ": ", as.integer(table(report$status))),
  paste0("  informative (filter removed something): ",
    sum(report$informative %in% TRUE), " of ",
    sum(report$status %in% c("PASS", "FAIL"))
  )
)

writeLines(txt, paste0(OUT, ".txt"))
utils::write.csv(report, paste0(OUT, ".csv"), row.names = FALSE)

cat("\n")
cat(txt, sep = "\n")
cat("\n\nWritten:\n  ", OUT, ".txt   <- paste this back\n  ",
  OUT, ".csv   <- or attach this\n",
  sep = ""
)
