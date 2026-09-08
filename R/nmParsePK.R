## Internal machinery for translating a NONMEM $PK block into R source.
##
## The pipeline is: read and normalise the control stream -> extract a record ->
## lex -> parse into statements -> classify symbols -> derive covariate
## references -> emit R. Everything here is internal; `createParamFunction()` is
## the entry point.
##
## The parser is deliberately narrow. It accepts assignments, IF statements and
## closed-form arithmetic, and refuses everything else with a file:line message,
## because a silently mistranslated $PK block produces a plausible-looking but
## wrong Forest plot.

## ---------------------------------------------------------------------------
## Reading and record extraction
## ---------------------------------------------------------------------------

#' Read a control stream and split code from comments
#'
#' Returns a data frame with one row per physical line: `lineno`, `code` (the
#' line with its comment and trailing whitespace removed) and `comment` (the
#' text after the first `;`, or `""`).
#'
#' @noRd
nmReadModel <- function(modFile) {
  if (!file.exists(modFile)) stop("Model file not found: ", modFile)
  raw <- readLines(modFile, warn = FALSE)

  # A ";" always starts a comment in NONMEM abbreviated code; there are no
  # string literals to protect it from.
  pos     <- regexpr(";", raw, fixed = TRUE)
  code    <- ifelse(pos > 0, substr(raw, 1, pos - 1), raw)
  comment <- ifelse(pos > 0, substr(raw, pos + 1, nchar(raw)), "")

  data.frame(
    lineno  = seq_along(raw),
    code    = trimws(code, which = "right"),
    comment = trimws(comment),
    stringsAsFactors = FALSE
  )
}

#' Extract the lines belonging to one control-stream record
#'
#' `record` is a regular expression matched against the start of the line, e.g.
#' `"\\$PK"` or `"\\$INP(U(T)?)?"`. Matching is case-insensitive. The record ends
#' at the next line starting with `$`. Returns the rows of `mod` for that record
#' with the record name itself stripped from the first line, or a zero-row frame
#' when the record is absent.
#'
#' @noRd
nmRecord <- function(mod, record) {
  starts <- grep(paste0("^\\s*", record), mod$code, ignore.case = TRUE)
  if (length(starts) == 0) return(mod[0, , drop = FALSE])

  start <- starts[1]
  rest  <- grep("^\\s*\\$", mod$code)
  rest  <- rest[rest > start]
  stop_ <- if (length(rest) == 0) nrow(mod) else rest[1] - 1

  out <- mod[start:stop_, , drop = FALSE]
  out$code[1] <- sub(paste0("^\\s*", record), "", out$code[1], ignore.case = TRUE)
  out
}

#' Column names declared in $INPUT
#'
#' Handles the `$INP` / `$INPU` / `$INPUT` abbreviations, `SYNONYM=REAL` pairs,
#' `=DROP` / `=SKIP` columns and continuation across lines.
#'
#' @noRd
nmInputNames <- function(mod) {
  rec <- nmRecord(mod, "\\$INP(U(T)?)?\\b")
  if (nrow(rec) == 0) return(character(0))

  items <- unlist(strsplit(trimws(paste(rec$code, collapse = " ")), "[[:space:],]+"))
  items <- items[nzchar(items)]
  if (length(items) == 0) return(character(0))

  # DROP/SKIP columns are not read by NONMEM, so they are not covariates.
  items <- items[!grepl("(^|=)(DROP|SKIP)$", items, ignore.case = TRUE)]
  # "SYNONYM=REAL" declares one column reachable under either name; keep both.
  unique(unlist(strsplit(items, "=", fixed = TRUE)))
}

#' Column names declared by $INPUT, one per position
#'
#' NONMEM reads a data file positionally: the header line is skipped and
#' `$INPUT` names the columns by position, so the names in the file need not
#' agree with the names the model uses. This returns one entry per position,
#' keeping `DROP`/`SKIP` columns because they still occupy one, together with
#' the alternate names introduced by `SYNONYM=REAL` pairs.
#'
#' Returns `list(names = <character, one per position>, aliases = <named
#' character, alternate name -> primary name>)`.
#'
#' @noRd
nmInputPositions <- function(mod) {
  rec <- nmRecord(mod, "\\$INP(U(T)?)?\\b")
  if (nrow(rec) == 0) return(list(names = character(0), aliases = character(0)))

  items <- unlist(strsplit(trimws(paste(rec$code, collapse = " ")), "[[:space:],]+"))
  items <- items[nzchar(items)]

  nms     <- character(length(items))
  aliases <- character(0)
  for (i in seq_along(items)) {
    parts <- strsplit(items[i], "=", fixed = TRUE)[[1]]
    nms[i] <- parts[1]
    if (length(parts) > 1 && !grepl("^(DROP|SKIP)$", parts[2], ignore.case = TRUE)) {
      # SYNONYM=REAL: either name refers to this column.
      aliases[parts[2]] <- parts[1]
    }
  }
  list(names = nms, aliases = aliases)
}

#' Number of THETAs declared by the $THETA records
#'
#' Counts actual THETAs rather than `$THETA` lines: a `(low,init,up)` triplet is
#' one THETA, `FIX`/`FIXED` flags are dropped, and an `xN` repeat count expands.
#'
#' @noRd
nmCountThetas <- function(mod) {
  starts <- grep("^\\s*\\$THE(T(A)?)?\\b", mod$code, ignore.case = TRUE)
  if (length(starts) == 0) return(0L)

  allStarts <- grep("^\\s*\\$", mod$code)
  n <- 0L
  for (s in starts) {
    nxt   <- allStarts[allStarts > s]
    stop_ <- if (length(nxt) == 0) nrow(mod) else nxt[1] - 1
    txt   <- paste(mod$code[s:stop_], collapse = " ")
    txt   <- sub("^\\s*\\$THE(T(A)?)?\\b", "", txt, ignore.case = TRUE)
    n     <- n + nmCountThetaRecord(txt)
  }
  n
}

#' @noRd
nmCountThetaRecord <- function(txt) {
  # `code` is normally comment-free already; strip defensively so a stray ";"
  # can never contribute a spurious theta.
  txt <- sub(";.*$", "", txt)
  txt <- gsub("\\b(FIXED|FIX)\\b", " ", txt, ignore.case = TRUE)
  # A parenthesised group is one THETA unless it carries an explicit xN count.
  n <- 0L
  while (grepl("\\(", txt)) {
    m <- regexpr("\\([^()]*\\)\\s*(x\\s*[0-9]+)?", txt, ignore.case = TRUE)
    if (m < 0) break
    grp <- substr(txt, m, m + attr(m, "match.length") - 1)
    rep <- sub(".*[xX]\\s*([0-9]+)\\s*$", "\\1", grp)
    n   <- n + if (grepl("[xX]\\s*[0-9]+\\s*$", grp)) as.integer(rep) else 1L
    txt <- paste0(substr(txt, 1, m - 1), " ", substr(txt, m + attr(m, "match.length"), nchar(txt)))
  }
  # Remaining bare values, each optionally with an xN repeat count.
  toks <- unlist(strsplit(trimws(txt), "\\s+"))
  toks <- toks[nzchar(toks)]
  i <- 1L
  while (i <= length(toks)) {
    if (grepl("^[xX][0-9]+$", toks[i])) {
      n <- n + as.integer(sub("^[xX]", "", toks[i])) - 1L
    } else if (grepl("^[-+]?([0-9]*\\.)?[0-9]+([EeDd][-+]?[0-9]+)?$", toks[i])) {
      n <- n + 1L
    }
    i <- i + 1L
  }
  as.integer(n)
}

## ---------------------------------------------------------------------------
## Lexer
## ---------------------------------------------------------------------------

## Dot operators are matched before numbers so that ".EQ." is never mistaken for
## a fractional part, and longest-first so ".EQN." beats ".EQ.".
nmDotOps <- c(
  ".EQN." = "==", ".NEN." = "!=",
  ".AND." = "&",  ".NOT." = "!", ".OR." = "|",
  ".EQ."  = "==", ".NE."  = "!=", ".GE." = ">=",
  ".GT."  = ">",  ".LE."  = "<=", ".LT." = "<"
)

nmFunctions <- c(
  EXP = "exp", LOG = "log", LOG10 = "log10", SQRT = "sqrt", ABS = "abs",
  MIN = "min", MAX = "max", INT = "trunc",
  SIN = "sin", COS = "cos", TAN = "tan", ATAN = "atan", ASIN = "asin",
  ACOS = "acos", GAMLN = "lgamma"
)

#' Tokenise one line of NONMEM abbreviated code
#'
#' @noRd
nmLex <- function(text, lineno, modFile) {
  toks <- list()
  i    <- 1L
  n    <- nchar(text)

  bad <- function(ch) {
    stop("Unsupported NONMEM construct in $PK at ", basename(modFile), ":", lineno,
         " - unrecognised character '", ch, "' in:\n  ", trimws(text),
         "\ncreateParamFunction() handles assignments, IF statements and ",
         "closed-form arithmetic only.", call. = FALSE)
  }

  while (i <= n) {
    ch <- substr(text, i, i)

    if (grepl("\\s", ch)) { i <- i + 1L; next }

    # Dot operator, checked before a number so ".EQ." is unambiguous.
    matched <- FALSE
    if (ch == ".") {
      for (op in names(nmDotOps)) {
        if (toupper(substr(text, i, i + nchar(op) - 1L)) == op) {
          toks[[length(toks) + 1L]] <- list(type = "op", value = nmDotOps[[op]])
          i <- i + nchar(op)
          matched <- TRUE
          break
        }
      }
    }
    if (matched) next

    # Number, including Fortran D exponents.
    if (grepl("[0-9]", ch) || (ch == "." && grepl("[0-9]", substr(text, i + 1L, i + 1L)))) {
      m <- regexpr("^([0-9]*\\.)?[0-9]+([EeDd][-+]?[0-9]+)?", substr(text, i, n))
      lit <- substr(text, i, i + attr(m, "match.length") - 1L)
      toks[[length(toks) + 1L]] <- list(
        type = "num", value = as.numeric(sub("[Dd]", "E", lit))
      )
      i <- i + attr(m, "match.length")
      next
    }

    # Identifier.
    if (grepl("[A-Za-z_]", ch)) {
      m <- regexpr("^[A-Za-z_][A-Za-z0-9_]*", substr(text, i, n))
      toks[[length(toks) + 1L]] <- list(
        type = "sym", value = substr(text, i, i + attr(m, "match.length") - 1L)
      )
      i <- i + attr(m, "match.length")
      next
    }

    # Operators and punctuation.
    two <- substr(text, i, i + 1L)
    if (two == "**") {
      toks[[length(toks) + 1L]] <- list(type = "op", value = "^"); i <- i + 2L; next
    }
    if (two %in% c("==", "/=", "<=", ">=")) {
      toks[[length(toks) + 1L]] <- list(type = "op", value = sub("/=", "!=", two))
      i <- i + 2L; next
    }
    if (ch %in% c("+", "-", "*", "/", "^", "<", ">")) {
      toks[[length(toks) + 1L]] <- list(type = "op", value = ch); i <- i + 1L; next
    }
    if (ch == "=") { toks[[length(toks) + 1L]] <- list(type = "assign"); i <- i + 1L; next }
    if (ch == "(") { toks[[length(toks) + 1L]] <- list(type = "lparen"); i <- i + 1L; next }
    if (ch == ")") { toks[[length(toks) + 1L]] <- list(type = "rparen"); i <- i + 1L; next }
    if (ch == ",") { toks[[length(toks) + 1L]] <- list(type = "comma");  i <- i + 1L; next }

    bad(ch)
  }
  toks
}

## ---------------------------------------------------------------------------
## Expression parser (recursive descent, precedence climbing)
## ---------------------------------------------------------------------------

## Binary operator precedence, low to high. `^` is right-associative, matching
## both Fortran's `**` and R's `^`.
nmBinPrec <- c("|" = 1, "&" = 2,
               "==" = 4, "!=" = 4, "<" = 4, ">" = 4, "<=" = 4, ">=" = 4,
               "+" = 5, "-" = 5, "*" = 6, "/" = 6, "%%" = 6, "^" = 8)

#' Parser state: a token list plus a cursor
#' @noRd
nmParser <- function(toks, lineno, modFile) {
  env <- new.env(parent = emptyenv())
  env$toks <- toks
  env$pos  <- 1L
  env$lineno <- lineno
  env$modFile <- modFile
  env
}

#' @noRd
nmPeek <- function(p) if (p$pos <= length(p$toks)) p$toks[[p$pos]] else NULL

#' @noRd
nmNext <- function(p) { t <- nmPeek(p); p$pos <- p$pos + 1L; t }

#' @noRd
nmFail <- function(p, msg) {
  stop("Unsupported NONMEM construct in $PK at ", basename(p$modFile), ":",
       p$lineno, " - ", msg,
       "\ncreateParamFunction() handles assignments, IF statements and ",
       "closed-form arithmetic only.", call. = FALSE)
}

#' Parse an expression with precedence climbing
#' @noRd
nmParseExpr <- function(p, minPrec = 0) {
  lhs <- nmParseUnary(p)
  repeat {
    t <- nmPeek(p)
    if (is.null(t) || t$type != "op" || !(t$value %in% names(nmBinPrec))) break
    prec <- nmBinPrec[[t$value]]
    if (prec < minPrec) break
    nmNext(p)
    # `^` is right-associative; everything else binds left.
    rhs <- nmParseExpr(p, if (t$value == "^") prec else prec + 1)
    lhs <- list(type = "binop", op = t$value, lhs = lhs, rhs = rhs)
  }
  lhs
}

#' @noRd
nmParseUnary <- function(p) {
  t <- nmPeek(p)
  if (!is.null(t) && t$type == "op" && t$value %in% c("-", "+", "!")) {
    nmNext(p)
    # Unary minus binds looser than `^`, as in both Fortran and R: -2**2 == -4.
    arg <- nmParseExpr(p, if (t$value == "!") 3 else 7)
    if (t$value == "+") return(arg)
    return(list(type = "unop", op = t$value, arg = arg))
  }
  nmParseAtom(p)
}

#' @noRd
nmParseAtom <- function(p) {
  t <- nmNext(p)
  if (is.null(t)) nmFail(p, "unexpected end of expression")

  if (t$type == "num") return(list(type = "num", value = t$value))

  if (t$type == "lparen") {
    e <- nmParseExpr(p)
    cl <- nmNext(p)
    if (is.null(cl) || cl$type != "rparen") nmFail(p, "unbalanced parentheses")
    return(e)
  }

  if (t$type == "sym") {
    nm <- toupper(t$value)
    nxt <- nmPeek(p)

    if (!is.null(nxt) && nxt$type == "lparen") {
      nmNext(p)
      args <- list()
      if (!is.null(nmPeek(p)) && nmPeek(p)$type != "rparen") {
        repeat {
          args[[length(args) + 1L]] <- nmParseExpr(p)
          nx <- nmPeek(p)
          if (!is.null(nx) && nx$type == "comma") { nmNext(p); next }
          break
        }
      }
      cl <- nmNext(p)
      if (is.null(cl) || cl$type != "rparen") nmFail(p, "unbalanced parentheses")

      if (nm == "THETA") {
        if (length(args) != 1 || args[[1]]$type != "num") {
          nmFail(p, "THETA() index must be a literal integer")
        }
        return(list(type = "theta", index = as.integer(args[[1]]$value)))
      }
      if (nm == "ETA") {
        if (length(args) != 1 || args[[1]]$type != "num") {
          nmFail(p, "ETA() index must be a literal integer")
        }
        return(list(type = "eta", index = as.integer(args[[1]]$value)))
      }
      if (nm %in% c("EPS", "ERR")) {
        nmFail(p, paste0(nm, "() cannot appear in $PK"))
      }
      if (nm == "A") {
        nmFail(p, "A() refers to a compartment amount and needs an ODE solution")
      }
      if (nm == "MOD") {
        # R spells the remainder as an infix operator; emitting it as a call
        # ("%%(a, b)") would produce source that does not parse.
        if (length(args) != 2) nmFail(p, "MOD() takes exactly two arguments")
        return(list(type = "binop", op = "%%", lhs = args[[1]], rhs = args[[2]]))
      }
      if (nm %in% names(nmFunctions)) {
        return(list(type = "call", fn = nmFunctions[[nm]], args = args))
      }
      nmFail(p, paste0("unsupported function '", t$value, "()'"))
    }

    return(list(type = "sym", name = t$value))
  }

  nmFail(p, paste0("unexpected token of type '", t$type, "'"))
}

## ---------------------------------------------------------------------------
## Statement parser
## ---------------------------------------------------------------------------

#' Parse the code lines of a record into a statement list
#'
#' Statements are `assign` (lhs, rhs, lineno, comment) and `if` (cond, then,
#' elifs, else_, lineno). Blocks nest.
#'
#' @noRd
nmParseStatements <- function(rec, modFile) {
  keep <- nzchar(trimws(rec$code))
  rec  <- rec[keep, , drop = FALSE]

  verbatim <- grep('^\\s*"', rec$code)
  if (length(verbatim) > 0) {
    stop("Unsupported NONMEM construct in $PK at ", basename(modFile), ":",
         rec$lineno[verbatim[1]], " - verbatim FORTRAN code.", call. = FALSE)
  }

  state <- new.env(parent = emptyenv())
  state$i <- 1L
  res <- nmParseBlock(rec, state, modFile, terminators = character(0))
  if (state$i <= nrow(rec)) {
    stop("Unsupported NONMEM construct in $PK at ", basename(modFile), ":",
         rec$lineno[state$i], " - unexpected '", trimws(rec$code[state$i]), "'.",
         call. = FALSE)
  }
  res
}

#' @noRd
nmParseBlock <- function(rec, state, modFile, terminators) {
  stmts <- list()
  while (state$i <= nrow(rec)) {
    line   <- rec$code[state$i]
    lineno <- rec$lineno[state$i]
    up     <- toupper(trimws(line))

    if (any(vapply(terminators, function(tm) grepl(tm, up), logical(1)))) break

    # Block IF: "IF (cond) THEN"
    if (grepl("^IF\\s*\\(.*\\)\\s*THEN$", up)) {
      state$i <- state$i + 1L
      cond <- nmParseCondition(line, lineno, modFile)
      terms <- c("^ELSE\\s+IF\\b", "^ELSE$", "^END\\s*IF$")
      thenStmts <- nmParseBlock(rec, state, modFile, terms)

      elifs <- list()
      elseStmts <- NULL
      repeat {
        if (state$i > nrow(rec)) {
          stop("Unsupported NONMEM construct in $PK at ", basename(modFile), ":",
               lineno, " - IF block is never closed by END IF.", call. = FALSE)
        }
        cur   <- rec$code[state$i]
        curUp <- toupper(trimws(cur))
        if (grepl("^ELSE\\s+IF\\b", curUp)) {
          state$i <- state$i + 1L
          elifs[[length(elifs) + 1L]] <- list(
            cond  = nmParseCondition(sub("(?i)^\\s*ELSE\\s+", "", cur, perl = TRUE),
                                     rec$lineno[state$i - 1L], modFile),
            stmts = nmParseBlock(rec, state, modFile, terms)
          )
          next
        }
        if (grepl("^ELSE$", curUp)) {
          state$i <- state$i + 1L
          elseStmts <- nmParseBlock(rec, state, modFile, terms)
          next
        }
        if (grepl("^END\\s*IF$", curUp)) { state$i <- state$i + 1L; break }
        stop("Unsupported NONMEM construct in $PK at ", basename(modFile), ":",
             rec$lineno[state$i], " - '", trimws(cur), "' inside an IF block.",
             call. = FALSE)
      }

      stmts[[length(stmts) + 1L]] <- list(
        type = "if", cond = cond, then = thenStmts, elifs = elifs,
        else_ = elseStmts, lineno = lineno
      )
      next
    }

    # One-line IF: "IF (cond) VAR = expr"
    if (grepl("^IF\\s*\\(", up)) {
      close <- nmMatchParen(line, modFile, lineno)
      cond  <- nmParseCondition(substr(line, 1, close), lineno, modFile)
      body  <- trimws(substr(line, close + 1L, nchar(line)))
      if (!nzchar(body)) {
        stop("Unsupported NONMEM construct in $PK at ", basename(modFile), ":",
             lineno, " - IF without THEN and without a statement.", call. = FALSE)
      }
      inner <- nmParseAssign(body, lineno, rec$comment[state$i], modFile)
      stmts[[length(stmts) + 1L]] <- list(
        type = "if", cond = cond, then = list(inner), elifs = list(),
        else_ = NULL, lineno = lineno, oneline = TRUE
      )
      state$i <- state$i + 1L
      next
    }

    if (grepl("^(DO|WHILE|CALL|EXIT|GOTO|GO\\s+TO|RETURN)\\b", up)) {
      stop("Unsupported NONMEM construct in $PK at ", basename(modFile), ":",
           lineno, " - '", trimws(line), "'.", call. = FALSE)
    }

    stmts[[length(stmts) + 1L]] <-
      nmParseAssign(line, lineno, rec$comment[state$i], modFile)
    state$i <- state$i + 1L
  }
  stmts
}

#' Position of the parenthesis closing the one opened after IF
#' @noRd
nmMatchParen <- function(line, modFile, lineno) {
  chars <- strsplit(line, "")[[1]]
  depth <- 0L
  for (k in seq_along(chars)) {
    if (chars[k] == "(") depth <- depth + 1L
    if (chars[k] == ")") {
      depth <- depth - 1L
      if (depth == 0L) return(k)
    }
  }
  stop("Unsupported NONMEM construct in $PK at ", basename(modFile), ":", lineno,
       " - unbalanced parentheses in '", trimws(line), "'.", call. = FALSE)
}

#' @noRd
nmParseCondition <- function(line, lineno, modFile) {
  txt <- sub("(?i)^\\s*IF\\s*", "", line, perl = TRUE)
  txt <- sub("(?i)\\s*THEN\\s*$", "", txt, perl = TRUE)
  p <- nmParser(nmLex(txt, lineno, modFile), lineno, modFile)
  e <- nmParseExpr(p)
  if (p$pos <= length(p$toks)) nmFail(p, "trailing tokens in IF condition")
  e
}

#' @noRd
nmParseAssign <- function(line, lineno, comment, modFile) {
  toks <- nmLex(line, lineno, modFile)
  eq   <- which(vapply(toks, function(t) t$type == "assign", logical(1)))
  if (length(eq) == 0) {
    stop("Unsupported NONMEM construct in $PK at ", basename(modFile), ":",
         lineno, " - '", trimws(line), "' is not an assignment.", call. = FALSE)
  }
  eq <- eq[1]
  if (eq != 2L || toks[[1]]$type != "sym") {
    stop("Unsupported NONMEM construct in $PK at ", basename(modFile), ":",
         lineno, " - only assignment to a plain variable is supported, got '",
         trimws(line), "'.", call. = FALSE)
  }
  p <- nmParser(toks[(eq + 1L):length(toks)], lineno, modFile)
  rhs <- nmParseExpr(p)
  if (p$pos <= length(p$toks)) nmFail(p, "trailing tokens after the assignment")

  list(type = "assign", lhs = toks[[1]]$value, rhs = rhs,
       lineno = lineno, comment = comment)
}

## ---------------------------------------------------------------------------
## Deparse to R
## ---------------------------------------------------------------------------

## Precedence of the emitted R operators, used to decide where parentheses are
## needed. Matches R's own table.
nmRPrec <- c("|" = 1, "&" = 2, "!" = 3,
             "==" = 4, "!=" = 4, "<" = 4, ">" = 4, "<=" = 4, ">=" = 4,
             "+" = 5, "-" = 5, "*" = 6, "/" = 6, "%%" = 6, "u-" = 7, "^" = 8)

#' @noRd
nmPrecOf <- function(node) {
  switch(node$type,
    binop = nmRPrec[[node$op]],
    unop  = if (node$op == "-") nmRPrec[["u-"]] else nmRPrec[["!"]],
    99
  )
}

#' Render a parsed `$PK` expression node as R source
#'
#' @description Turns one expression node from the `statements` tree of
#'   [nmParsePK()] into a string of R code. `THETA(n)` becomes `<thetaVar>[n]`
#'   and every `ETA(n)` becomes the literal text `etaValue` (`"0"` gives typical
#'   values). Symbols and covariate names are emitted unchanged, `**` becomes
#'   `^`, and parentheses are added only where operator precedence needs them.
#'
#' @param node An expression node - a `rhs` or `cond` from an [nmParsePK()]
#'   statement, or a sub-node of one.
#' @param thetaVar Name of the vector that `THETA(n)` indexes into. Default
#'   `"thetas"`.
#' @param etaValue Text substituted for every `ETA(n)`. Default `"0"`. Pass, for
#'   example, `"eta[3]"` to keep the random effect.
#'
#' @return A length-one character string of R source.
#'
#' @seealso [nmParsePK()].
#'
#' @export
#'
#' @examples
#' modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")
#' p <- nmParsePK(modFile)
#' # render the right-hand side of the first assignment
#' first <- Find(function(s) s$type == "assign", p$statements)
#' nmDeparse(first$rhs, thetaVar = "thetas")
nmDeparse <- function(node, thetaVar = "thetas", etaValue = "0") {
  wrap <- function(child, parentPrec, side = c("left", "right")) {
    side <- match.arg(side)
    txt  <- nmDeparse(child, thetaVar, etaValue)
    cp   <- nmPrecOf(child)
    need <- cp < parentPrec
    # Left-associative operators need parentheses on the right at equal
    # precedence (a - (b - c)); `^` is right-associative, so the reverse.
    if (!need && cp == parentPrec && child$type == "binop") {
      need <- if (parentPrec == nmRPrec[["^"]]) side == "left" else side == "right"
    }
    if (need) paste0("(", txt, ")") else txt
  }

  switch(node$type,
    num   = nmFormatNum(node$value),
    sym   = node$name,
    theta = paste0(thetaVar, "[", node$index, "]"),
    eta   = etaValue,
    call  = paste0(node$fn, "(",
                   paste(vapply(node$args, nmDeparse, character(1),
                                thetaVar, etaValue), collapse = ", "), ")"),
    unop  = paste0(node$op, wrap(node$arg, nmPrecOf(node), "right")),
    binop = {
      prec <- nmRPrec[[node$op]]
      sep  <- if (node$op == "^") "" else " "
      paste0(wrap(node$lhs, prec, "left"), sep, node$op, sep,
             wrap(node$rhs, prec, "right"))
    },
    stop("Internal error: cannot deparse node type '", node$type, "'.")
  )
}

#' Fold the constants introduced by setting ETA() to 0
#'
#' Substituting `ETA(n) -> 0` leaves artefacts such as `TVCL * exp(0)`. Folding
#' them away makes the emitted source read like the control stream, which is the
#' point of generating it. Only identities that hold for every real value are
#' applied, so this cannot change the arithmetic.
#'
#' @noRd
nmSimplify <- function(node) {
  isNum <- function(n, v) n$type == "num" && !is.na(n$value) && n$value == v

  # Typical values: substitute ETA() in the tree rather than at deparse time, so
  # the folding below can see the resulting constants.
  if (node$type == "eta") return(list(type = "num", value = 0))

  if (node$type == "call") {
    node$args <- lapply(node$args, nmSimplify)
    if (node$fn == "exp" && length(node$args) == 1L && isNum(node$args[[1]], 0)) {
      return(list(type = "num", value = 1))
    }
    if (node$fn == "log" && length(node$args) == 1L && isNum(node$args[[1]], 1)) {
      return(list(type = "num", value = 0))
    }
    return(node)
  }

  if (node$type == "unop") {
    node$arg <- nmSimplify(node$arg)
    if (node$op == "-" && node$arg$type == "num") {
      return(list(type = "num", value = -node$arg$value))
    }
    return(node)
  }

  if (node$type == "binop") {
    node$lhs <- nmSimplify(node$lhs)
    node$rhs <- nmSimplify(node$rhs)
    op <- node$op

    if (op == "*") {
      if (isNum(node$lhs, 1)) return(node$rhs)
      if (isNum(node$rhs, 1)) return(node$lhs)
    }
    if (op == "+") {
      if (isNum(node$lhs, 0)) return(node$rhs)
      if (isNum(node$rhs, 0)) return(node$lhs)
    }
    if (op == "-" && isNum(node$rhs, 0)) return(node$lhs)
    if (op == "/" && isNum(node$rhs, 1)) return(node$lhs)
    if (op == "^" && isNum(node$rhs, 1)) return(node$lhs)
    return(node)
  }

  node
}

#' Apply `nmSimplify()` across a statement list
#' @noRd
nmSimplifyStmts <- function(stmts) {
  lapply(stmts, function(s) {
    if (s$type == "assign") {
      # Recorded before folding, so the emitted source can still note where an
      # ETA() was dropped even though `exp(0)` has been simplified away.
      s$hadEta <- nmHasEta(s$rhs)
      s$rhs    <- nmSimplify(s$rhs)
    } else {
      s$cond  <- nmSimplify(s$cond)
      s$then  <- nmSimplifyStmts(s$then)
      s$elifs <- lapply(s$elifs, function(e) {
        e$cond  <- nmSimplify(e$cond)
        e$stmts <- nmSimplifyStmts(e$stmts)
        e
      })
      if (!is.null(s$else_)) s$else_ <- nmSimplifyStmts(s$else_)
    }
    s
  })
}

#' Format a numeric literal as the shortest string that reads back exactly
#'
#' @description Keeps ordinary values readable (`75`, `0.75`) while allowing
#'   scientific notation where writing the digits out would be absurd
#'   (`1e-12`). Exported for emitters built on [nmParsePK()] that need to write
#'   covariate reference values into generated source.
#'
#' @param x A single numeric value (or `NA`).
#'
#' @return A length-one character string that `as.numeric()` maps back to `x`.
#'
#' @seealso [nmParsePK()].
#'
#' @export
#'
#' @examples
#' nmFormatNum(75)
#' nmFormatNum(1e-12)
nmFormatNum <- function(x) {
  if (is.na(x)) return("NA")
  for (d in 1:17) {
    s <- format(x, digits = d)
    if (identical(as.numeric(s), x)) return(s)
  }
  format(x, digits = 17)
}

## ---------------------------------------------------------------------------
## Symbol classification
## ---------------------------------------------------------------------------

#' Walk a statement list, collecting assigned and referenced symbols
#'
#' Returns `list(assigned = <chr, in first-assignment order>, used = <chr>)`.
#'
#' @noRd
nmSymbols <- function(stmts) {
  assigned <- character(0)
  used     <- character(0)

  walkExpr <- function(node) {
    switch(node$type,
      sym   = used <<- c(used, node$name),
      call  = lapply(node$args, walkExpr),
      unop  = walkExpr(node$arg),
      binop = { walkExpr(node$lhs); walkExpr(node$rhs) },
      NULL
    )
    invisible(NULL)
  }
  walkStmts <- function(ss) {
    for (s in ss) {
      if (s$type == "assign") {
        assigned <<- c(assigned, s$lhs)
        walkExpr(s$rhs)
      } else {
        walkExpr(s$cond)
        walkStmts(s$then)
        for (e in s$elifs) { walkExpr(e$cond); walkStmts(e$stmts) }
        if (!is.null(s$else_)) walkStmts(s$else_)
      }
    }
  }
  walkStmts(stmts)

  list(assigned = unique(assigned), used = unique(used))
}

#' Map parameters onto the ETA that carries their between-subject variability
#'
#' Detects the exponential-IIV idiom `P = <expr> * EXP(ETA(n))` at the top level
#' of an assignment and returns a named integer vector `P -> n`. NONMEM tables
#' write individual values, so this is what lets [verifyParamFunction()] divide
#' the tabled value by `exp(ETA)` and recover the typical value.
#'
#' Only this one pattern is recognised; parameters written any other way are
#' absent from the result and are reported as unverifiable.
#'
#' @noRd
nmEtaMap <- function(stmts) {
  map <- integer(0)

  etaInExp <- function(node) {
    # <expr> * EXP(ETA(n)) or EXP(ETA(n)) * <expr>
    if (node$type != "binop" || node$op != "*") return(NA_integer_)
    for (side in list(node$lhs, node$rhs)) {
      if (side$type == "call" && side$fn == "exp" && length(side$args) == 1L &&
          side$args[[1]]$type == "eta") {
        return(side$args[[1]]$index)
      }
    }
    NA_integer_
  }

  for (s in stmts) {
    if (s$type != "assign") next
    idx <- etaInExp(s$rhs)
    if (!is.na(idx)) map[s$lhs] <- idx
  }
  map
}

#' Highest THETA index referenced anywhere in a statement list
#' @noRd
nmMaxTheta <- function(stmts) {
  mx <- 0L
  walkExpr <- function(node) {
    switch(node$type,
      theta = mx <<- max(mx, node$index),
      call  = lapply(node$args, walkExpr),
      unop  = walkExpr(node$arg),
      binop = { walkExpr(node$lhs); walkExpr(node$rhs) },
      NULL
    )
    invisible(NULL)
  }
  walkStmts <- function(ss) {
    for (s in ss) {
      if (s$type == "assign") walkExpr(s$rhs)
      else {
        walkExpr(s$cond); walkStmts(s$then)
        for (e in s$elifs) { walkExpr(e$cond); walkStmts(e$stmts) }
        if (!is.null(s$else_)) walkStmts(s$else_)
      }
    }
  }
  walkStmts(stmts)
  mx
}
