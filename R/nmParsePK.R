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
  pos <- regexpr(";", raw, fixed = TRUE)
  code <- ifelse(pos > 0, substr(raw, 1, pos - 1), raw)
  comment <- ifelse(pos > 0, substr(raw, pos + 1, nchar(raw)), "")

  data.frame(
    lineno = seq_along(raw),
    code = trimws(code, which = "right"),
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
  if (length(starts) == 0) {
    return(mod[0, , drop = FALSE])
  }

  start <- starts[1]
  rest <- grep("^\\s*\\$", mod$code)
  rest <- rest[rest > start]
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
  if (nrow(rec) == 0) {
    return(character(0))
  }

  items <- unlist(strsplit(trimws(paste(rec$code, collapse = " ")), "[[:space:],]+"))
  items <- items[nzchar(items)]
  if (length(items) == 0) {
    return(character(0))
  }

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
#' character, alternate name -> primary name>, dropped = <logical, one per
#' position, TRUE for a DROP/SKIP column>)`.
#'
#' @noRd
nmInputPositions <- function(mod) {
  rec <- nmRecord(mod, "\\$INP(U(T)?)?\\b")
  if (nrow(rec) == 0) {
    return(list(
      names = character(0), aliases = character(0),
      dropped = logical(0)
    ))
  }

  items <- unlist(strsplit(trimws(paste(rec$code, collapse = " ")), "[[:space:],]+"))
  items <- items[nzchar(items)]

  nms <- character(length(items))
  aliases <- character(0)
  for (i in seq_along(items)) {
    parts <- strsplit(items[i], "=", fixed = TRUE)[[1]]
    nms[i] <- parts[1]
    if (length(parts) > 1 && !grepl("^(DROP|SKIP)$", parts[2], ignore.case = TRUE)) {
      # SYNONYM=REAL: either name refers to this column.
      aliases[parts[2]] <- parts[1]
    }
  }
  list(
    names = nms, aliases = aliases,
    dropped = grepl("(^|=)(DROP|SKIP)$", items, ignore.case = TRUE)
  )
}

#' Number of THETAs declared by the $THETA records
#'
#' Counts actual THETAs rather than `$THETA` lines: a `(low,init,up)` triplet is
#' one THETA, `FIX`/`FIXED` flags are dropped, and an `xN` repeat count expands.
#'
#' @noRd
nmCountThetas <- function(mod) {
  starts <- grep("^\\s*\\$THE(T(A)?)?\\b", mod$code, ignore.case = TRUE)
  if (length(starts) == 0) {
    return(0L)
  }

  allStarts <- grep("^\\s*\\$", mod$code)
  n <- 0L
  for (s in starts) {
    nxt <- allStarts[allStarts > s]
    stop_ <- if (length(nxt) == 0) nrow(mod) else nxt[1] - 1
    txt <- paste(mod$code[s:stop_], collapse = " ")
    txt <- sub("^\\s*\\$THE(T(A)?)?\\b", "", txt, ignore.case = TRUE)
    n <- n + nmCountThetaRecord(txt)
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
    n <- n + if (grepl("[xX]\\s*[0-9]+\\s*$", grp)) as.integer(rep) else 1L
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
  ".AND." = "&", ".NOT." = "!", ".OR." = "|",
  ".EQ." = "==", ".NE." = "!=", ".GE." = ">=",
  ".GT." = ">", ".LE." = "<=", ".LT." = "<"
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
nmLex <- function(text, lineno, modFile, block = "$PK") {
  toks <- list()
  i <- 1L
  n <- nchar(text)

  bad <- function(ch) {
    stop("Unsupported NONMEM construct in ", block, " at ", basename(modFile), ":", lineno,
      " - unrecognised character '", ch, "' in:\n  ", trimws(text),
      "\ncreateParamFunction() handles assignments, IF statements and ",
      "closed-form arithmetic only.",
      call. = FALSE
    )
  }

  while (i <= n) {
    ch <- substr(text, i, i)

    if (grepl("\\s", ch)) {
      i <- i + 1L
      next
    }

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
      toks[[length(toks) + 1L]] <- list(type = "op", value = "^")
      i <- i + 2L
      next
    }
    if (two %in% c("==", "/=", "<=", ">=")) {
      toks[[length(toks) + 1L]] <- list(type = "op", value = sub("/=", "!=", two))
      i <- i + 2L
      next
    }
    if (ch %in% c("+", "-", "*", "/", "^", "<", ">")) {
      toks[[length(toks) + 1L]] <- list(type = "op", value = ch)
      i <- i + 1L
      next
    }
    if (ch == "=") {
      toks[[length(toks) + 1L]] <- list(type = "assign")
      i <- i + 1L
      next
    }
    if (ch == "(") {
      toks[[length(toks) + 1L]] <- list(type = "lparen")
      i <- i + 1L
      next
    }
    if (ch == ")") {
      toks[[length(toks) + 1L]] <- list(type = "rparen")
      i <- i + 1L
      next
    }
    if (ch == ",") {
      toks[[length(toks) + 1L]] <- list(type = "comma")
      i <- i + 1L
      next
    }

    bad(ch)
  }
  toks
}

## ---------------------------------------------------------------------------
## Expression parser (recursive descent, precedence climbing)
## ---------------------------------------------------------------------------

## Binary operator precedence, low to high. `^` is right-associative, matching
## both Fortran's `**` and R's `^`.
nmBinPrec <- c(
  "|" = 1, "&" = 2,
  "==" = 4, "!=" = 4, "<" = 4, ">" = 4, "<=" = 4, ">=" = 4,
  "+" = 5, "-" = 5, "*" = 6, "/" = 6, "^" = 8
)

#' Parser state: a token list plus a cursor
#' @noRd
nmParser <- function(toks, lineno, modFile, block = "$PK") {
  env <- new.env(parent = emptyenv())
  env$toks <- toks
  env$pos <- 1L
  env$lineno <- lineno
  env$modFile <- modFile
  env$block <- block
  env
}

#' @noRd
nmPeek <- function(p) if (p$pos <= length(p$toks)) p$toks[[p$pos]] else NULL

#' @noRd
nmNext <- function(p) {
  t <- nmPeek(p)
  p$pos <- p$pos + 1L
  t
}

#' @noRd
nmFail <- function(p, msg) {
  stop("Unsupported NONMEM construct in ", p$block, " at ", basename(p$modFile), ":",
    p$lineno, " - ", msg,
    "\ncreateParamFunction() handles assignments, IF statements and ",
    "closed-form arithmetic only.",
    call. = FALSE
  )
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
    if (t$value == "+") {
      return(arg)
    }
    return(list(type = "unop", op = t$value, arg = arg))
  }
  nmParseAtom(p)
}

#' @noRd
nmParseAtom <- function(p) {
  t <- nmNext(p)
  if (is.null(t)) nmFail(p, "unexpected end of expression")

  if (t$type == "num") {
    return(list(type = "num", value = t$value))
  }

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
          if (!is.null(nx) && nx$type == "comma") {
            nmNext(p)
            next
          }
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
        ## $PRED builds Y from the residual error in the same block as the
        ## parameters, so the block cannot parse without these. $PK has no
        ## business containing one - NM-TRAN does not allow it either - so
        ## there the refusal stands.
        ##
        ## ERR() is ambiguous by NM-TRAN's own definition: EPS(n) for
        ## population data, ETA(n) for single-subject. Both fold to 0, so
        ## folding is safe either way, but it can never earn an etaMap entry.
        if (!identical(p$block, "$PRED")) {
          nmFail(p, paste0(nm, "() cannot appear in ", p$block))
        }
        if (length(args) != 1 || args[[1]]$type != "num") {
          nmFail(p, paste0(nm, "() index must be a literal integer"))
        }
        return(list(type = "eps", index = as.integer(args[[1]]$value)))
      }
      if (nm == "A") {
        nmFail(p, "A() refers to a compartment amount and needs an ODE solution")
      }
      if (nm == "MOD") {
        if (length(args) != 2) nmFail(p, "MOD() takes exactly two arguments")
        # Not `%%`. Two reasons, either of which alone is a wrong answer:
        #   * Fortran MOD() truncates towards zero, R's `%%` floors, so they
        #     disagree whenever the first argument is negative -
        #     MOD(-7, 3) is -1 but -7 %% 3 is 2.
        #   * `%any%` binds tighter than `*` and `/` in R but MOD() is a call,
        #     so `MOD(A*B, C)` would have emitted `A * B %% C` = A * (B %% C).
        # Desugaring to `a - b * trunc(a / b)` is the faithful translation and,
        # being built from nodes the deparser already knows, it gets its
        # parentheses from the ordinary precedence rules. $PK expressions are
        # pure, so evaluating `a` and `b` twice is safe.
        return(list(
          type = "binop", op = "-", lhs = args[[1]],
          rhs = list(
            type = "binop", op = "*", lhs = args[[2]],
            rhs = list(
              type = "call", fn = "trunc",
              args = list(list(
                type = "binop",
                op = "/",
                lhs = args[[1]],
                rhs = args[[2]]
              ))
            )
          )
        ))
      }
      if (nm %in% names(nmFunctions)) {
        return(list(type = "call", fn = nmFunctions[[nm]], args = args))
      }
      ## OMEGA(i,j) and SIGMA(i,j) are NONMEM's own variance matrices. They are
      ## constants at the final estimates, so they are recorded here and folded
      ## to literals once the .ext has been read.
      if (nm %in% c("OMEGA", "SIGMA")) {
        idx <- vapply(args, function(a) {
          if (is.list(a) && identical(a$type, "num")) a$value else NA_real_
        }, numeric(1))
        if (length(idx) != 2L || anyNA(idx)) {
          nmFail(p, paste0(
            nm, "() takes two literal indices, as in ", nm, "(2,2)"
          ))
        }
        return(list(
          type = "nmmatrix", mat = nm,
          i = as.integer(idx[1]), j = as.integer(idx[2])
        ))
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

## The record a model defines its parameters in.
##
## $PK for a PREDPP model, $PRED for one that codes the prediction by hand.
## They are the same abbreviated-code language; what differs is that $PRED
## computes the prediction and the residual error in the same block, so it also
## contains Y and EPS()/ERR().
##
## @noRd
nmParamBlock <- function(mod, modFile) {
  for (b in c("$PK", "$PRED")) {
    rec <- nmRecord(mod, paste0("\\", b, "\\b"))
    if (nrow(rec) > 0) {
      return(list(rec = rec, block = b))
    }
  }
  stop("No $PK or $PRED record found in ", basename(modFile),
    ". createParamFunction() reads the block a model defines its parameters ",
    "in, and this model has neither.",
    call. = FALSE
  )
}

#' Parse the code lines of a record into a statement list
#'
#' Statements are `assign` (lhs, rhs, lineno, comment) and `if` (cond, then,
#' elifs, else_, lineno). Blocks nest.
#'
#' @noRd
nmParseStatements <- function(rec, modFile, ignoreVerbatim = FALSE,
                              block = "$PK") {
  keep <- nzchar(trimws(rec$code))
  rec <- rec[keep, , drop = FALSE]

  verbatim <- grep('^\\s*"', rec$code)
  if (length(verbatim) > 0) {
    ## Verbatim FORTRAN can define variables the rest of the block reads, and
    ## nothing here can tell that apart from a solver directive that touches
    ## no parameter. Refusing is the default for that reason; a caller who has
    ## read the block can say it is inert.
    if (!ignoreVerbatim) {
      stop("Unsupported NONMEM construct in ", block, " at ", basename(modFile), ":",
        rec$lineno[verbatim[1]], " - verbatim FORTRAN code.",
        "\nIf it defines nothing ", block, " reads - a solver directive, say - pass ",
        "ignoreVerbatim = TRUE.",
        call. = FALSE
      )
    }
    rec <- rec[-verbatim, , drop = FALSE]
  }

  state <- new.env(parent = emptyenv())
  state$i <- 1L
  res <- nmParseBlock(rec, state, modFile,
    terminators = character(0), block = block
  )
  if (state$i <= nrow(rec)) {
    stop("Unsupported NONMEM construct in ", block, " at ", basename(modFile), ":",
      rec$lineno[state$i], " - unexpected '", trimws(rec$code[state$i]), "'.",
      call. = FALSE
    )
  }
  res
}

#' @noRd
nmParseBlock <- function(rec, state, modFile, terminators, block = "$PK") {
  stmts <- list()
  while (state$i <= nrow(rec)) {
    line <- rec$code[state$i]
    lineno <- rec$lineno[state$i]
    up <- toupper(trimws(line))

    if (any(vapply(terminators, function(tm) grepl(tm, up), logical(1)))) break

    # Block IF: "IF (cond) THEN"
    if (grepl("^IF\\s*\\(.*\\)\\s*THEN$", up)) {
      state$i <- state$i + 1L
      cond <- nmParseCondition(line, lineno, modFile, block)
      terms <- c("^ELSE\\s+IF\\b", "^ELSE$", "^END\\s*IF$")
      thenStmts <- nmParseBlock(rec, state, modFile, terms, block)

      elifs <- list()
      elseStmts <- NULL
      repeat {
        if (state$i > nrow(rec)) {
          stop("Unsupported NONMEM construct in ", block, " at ", basename(modFile), ":",
            lineno, " - IF block is never closed by END IF.",
            call. = FALSE
          )
        }
        cur <- rec$code[state$i]
        curUp <- toupper(trimws(cur))
        if (grepl("^ELSE\\s+IF\\b", curUp)) {
          state$i <- state$i + 1L
          elifs[[length(elifs) + 1L]] <- list(
            cond = nmParseCondition(
              sub("(?i)^\\s*ELSE\\s+", "", cur, perl = TRUE),
              rec$lineno[state$i - 1L], modFile, block
            ),
            stmts = nmParseBlock(rec, state, modFile, terms, block)
          )
          next
        }
        if (grepl("^ELSE$", curUp)) {
          state$i <- state$i + 1L
          elseStmts <- nmParseBlock(rec, state, modFile, terms, block)
          next
        }
        if (grepl("^END\\s*IF$", curUp)) {
          state$i <- state$i + 1L
          break
        }
        stop("Unsupported NONMEM construct in ", block, " at ", basename(modFile), ":",
          rec$lineno[state$i], " - '", trimws(cur), "' inside an IF block.",
          call. = FALSE
        )
      }

      stmts[[length(stmts) + 1L]] <- list(
        type = "if", cond = cond, then = thenStmts, elifs = elifs,
        else_ = elseStmts, lineno = lineno
      )
      next
    }

    # One-line IF: "IF (cond) VAR = expr"
    if (grepl("^IF\\s*\\(", up)) {
      close <- nmMatchParen(line, modFile, lineno, block)
      cond <- nmParseCondition(substr(line, 1, close), lineno, modFile, block)
      body <- trimws(substr(line, close + 1L, nchar(line)))
      if (!nzchar(body)) {
        stop("Unsupported NONMEM construct in ", block, " at ", basename(modFile), ":",
          lineno, " - IF without THEN and without a statement.",
          call. = FALSE
        )
      }
      inner <- nmParseAssign(
        body, lineno, rec$comment[state$i], modFile,
        block
      )
      stmts[[length(stmts) + 1L]] <- list(
        type = "if", cond = cond, then = list(inner), elifs = list(),
        else_ = NULL, lineno = lineno, oneline = TRUE
      )
      state$i <- state$i + 1L
      next
    }

    if (grepl("^(DO|WHILE|CALL|EXIT|GOTO|GO\\s+TO|RETURN)\\b", up)) {
      stop("Unsupported NONMEM construct in ", block, " at ", basename(modFile), ":",
        lineno, " - '", trimws(line), "'.",
        call. = FALSE
      )
    }

    stmts[[length(stmts) + 1L]] <-
      nmParseAssign(line, lineno, rec$comment[state$i], modFile, block)
    state$i <- state$i + 1L
  }
  stmts
}

#' Position of the parenthesis closing the one opened after IF
#' @noRd
nmMatchParen <- function(line, modFile, lineno, block = "$PK") {
  chars <- strsplit(line, "")[[1]]
  depth <- 0L
  for (k in seq_along(chars)) {
    if (chars[k] == "(") depth <- depth + 1L
    if (chars[k] == ")") {
      depth <- depth - 1L
      if (depth == 0L) {
        return(k)
      }
    }
  }
  stop("Unsupported NONMEM construct in ", block, " at ", basename(modFile), ":", lineno,
    " - unbalanced parentheses in '", trimws(line), "'.",
    call. = FALSE
  )
}

#' @noRd
nmParseCondition <- function(line, lineno, modFile, block = "$PK") {
  txt <- sub("(?i)^\\s*IF\\s*", "", line, perl = TRUE)
  txt <- sub("(?i)\\s*THEN\\s*$", "", txt, perl = TRUE)
  p <- nmParser(nmLex(txt, lineno, modFile, block), lineno, modFile, block)
  e <- nmParseExpr(p)
  if (p$pos <= length(p$toks)) nmFail(p, "trailing tokens in IF condition")
  e
}

#' @noRd
nmParseAssign <- function(line, lineno, comment, modFile, block = "$PK") {
  toks <- nmLex(line, lineno, modFile, block)
  eq <- which(vapply(toks, function(t) t$type == "assign", logical(1)))
  if (length(eq) == 0) {
    stop("Unsupported NONMEM construct in ", block, " at ", basename(modFile), ":",
      lineno, " - '", trimws(line), "' is not an assignment.",
      call. = FALSE
    )
  }
  eq <- eq[1]
  if (eq != 2L || toks[[1]]$type != "sym") {
    stop("Unsupported NONMEM construct in ", block, " at ", basename(modFile), ":",
      lineno, " - only assignment to a plain variable is supported, got '",
      trimws(line), "'.",
      call. = FALSE
    )
  }
  p <- nmParser(toks[(eq + 1L):length(toks)], lineno, modFile, block)
  rhs <- nmParseExpr(p)
  if (p$pos <= length(p$toks)) nmFail(p, "trailing tokens after the assignment")

  list(
    type = "assign", lhs = toks[[1]]$value, rhs = rhs,
    lineno = lineno, comment = comment
  )
}

## ---------------------------------------------------------------------------
## Deparse to R
## ---------------------------------------------------------------------------

## Precedence of the emitted R operators, used to decide where parentheses are
## needed. Matches R's own table.
nmRPrec <- c(
  "|" = 1, "&" = 2, "!" = 3,
  "==" = 4, "!=" = 4, "<" = 4, ">" = 4, "<=" = 4, ">=" = 4,
  "+" = 5, "-" = 5, "*" = 6, "/" = 6, "u-" = 7, "^" = 8
)

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
#' @param epsValue Text substituted for every `EPS(n)`/`ERR(n)`, which only a
#'   `$PRED` block can contain. Default `"0"`, so the residual error drops out
#'   of a typical-value function exactly as `ETA(n)` does.
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
nmDeparse <- function(node, thetaVar = "thetas", etaValue = "0",
                      epsValue = "0") {
  wrap <- function(child, parentPrec, side = c("left", "right")) {
    side <- match.arg(side)
    txt <- nmDeparse(child, thetaVar, etaValue, epsValue)
    cp <- nmPrecOf(child)
    need <- cp < parentPrec
    # Left-associative operators need parentheses on the right at equal
    # precedence (a - (b - c)); `^` is right-associative, so the reverse.
    if (!need && cp == parentPrec && child$type == "binop") {
      need <- if (parentPrec == nmRPrec[["^"]]) side == "left" else side == "right"
    }
    if (need) paste0("(", txt, ")") else txt
  }

  switch(node$type,
    num = nmFormatNum(node$value),
    sym = node$name,
    theta = paste0(thetaVar, "[", node$index, "]"),
    eta = etaValue,
    eps = epsValue,
    call = paste0(
      node$fn, "(",
      paste(vapply(
        node$args, nmDeparse, character(1),
        thetaVar, etaValue, epsValue
      ), collapse = ", "), ")"
    ),
    unop = paste0(node$op, wrap(node$arg, nmPrecOf(node), "right")),
    binop = {
      prec <- nmRPrec[[node$op]]
      sep <- if (node$op == "^") "" else " "
      paste0(
        wrap(node$lhs, prec, "left"), sep, node$op, sep,
        wrap(node$rhs, prec, "right")
      )
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
  if (node$type == "eta" || node$type == "eps") {
    return(list(type = "num", value = 0))
  }

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
      if (isNum(node$lhs, 1)) {
        return(node$rhs)
      }
      if (isNum(node$rhs, 1)) {
        return(node$lhs)
      }
      ## A folded ETA()/EPS() leaves a literal 0, and the IOV idiom multiplies
      ## it by an occasion indicator: EXP(ETA(1) + ETA(2)*(1-OCC) + ETA(3)*OCC)
      ## would otherwise emit exp(0 * (1 - OCC) + 0 * OCC) instead of
      ## collapsing, and nmTypicalMap() would not see that TVCL is the typical
      ## value. The term is structurally zero whatever the occasion is.
      ##
      ## This makes 0 * NA fold to 0 rather than propagate NA. For a typical
      ## value that is the right answer: the term contributes nothing.
      if (isNum(node$lhs, 0) || isNum(node$rhs, 0)) {
        return(list(type = "num", value = 0))
      }
    }
    if (op == "+") {
      if (isNum(node$lhs, 0)) {
        return(node$rhs)
      }
      if (isNum(node$rhs, 0)) {
        return(node$lhs)
      }
    }
    if (op == "-" && isNum(node$rhs, 0)) {
      return(node$lhs)
    }
    if (op == "/" && isNum(node$rhs, 1)) {
      return(node$lhs)
    }
    if (op == "^" && isNum(node$rhs, 1)) {
      return(node$lhs)
    }
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
      s$rhs <- nmSimplify(s$rhs)
    } else {
      s$cond <- nmSimplify(s$cond)
      s$then <- nmSimplifyStmts(s$then)
      s$elifs <- lapply(s$elifs, function(e) {
        e$cond <- nmSimplify(e$cond)
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
  if (is.na(x)) {
    return("NA")
  }
  for (d in 1:17) {
    s <- format(x, digits = d)
    if (identical(as.numeric(s), x)) {
      return(s)
    }
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
  used <- character(0)

  walkExpr <- function(node) {
    switch(node$type,
      sym = used <<- c(used, node$name),
      call = lapply(node$args, walkExpr),
      unop = walkExpr(node$arg),
      binop = {
        walkExpr(node$lhs)
        walkExpr(node$rhs)
      },
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
        for (e in s$elifs) {
          walkExpr(e$cond)
          walkStmts(e$stmts)
        }
        if (!is.null(s$else_)) walkStmts(s$else_)
      }
    }
  }
  walkStmts(stmts)

  list(assigned = unique(assigned), used = unique(used))
}

#' Symbols read by one `$PK` expression
#' @noRd
nmExprSyms <- function(node) {
  out <- character(0)
  walk <- function(n) {
    switch(n$type,
      sym = out <<- c(out, n$name),
      call = lapply(n$args, walk),
      unop = walk(n$arg),
      binop = {
        walk(n$lhs)
        walk(n$rhs)
      },
      NULL
    )
    invisible(NULL)
  }
  walk(node)
  unique(out)
}

#' The first symbol read before any path has bound it
#'
#' A reaching-assignment walk. `bound` is what holds a value on entry; after an
#' `IF` construct a name counts as bound if *any* branch assigns it.
#'
#' The weaker "any" rather than "every" is deliberate. PsN's scm writes
#' exhaustive one-line branches with no `ELSE` -
#' `IF(FORM.EQ.1) FRELFORM = 1` / `IF(FORM.EQ.0) FRELFORM = (1 + THETA(12))`,
#' as run7.mod itself does - which no static analysis can tell apart from a
#' genuinely non-exhaustive branch. Requiring every path to assign would reject
#' the package's own reference model, so a name assigned only under a condition
#' is accepted here; if no branch fires at run time the generated function
#' fails loudly with "object not found", which is the right outcome anyway.
#'
#' This is what tells a covariate apart from a local. NONMEM populates the
#' `$INPUT` data items before `$PK` runs, so a name read before `$PK` assigns
#' it is being read from the data record - whether or not `$PK` later
#' reassigns it. `IF(WT.EQ.-99) WT = 75` reads WT in the condition, so WT is
#' correctly found this way; classifying on "is it assigned anywhere" instead
#' made that idiom - rule 1 of [nmCovRef()] - unreachable.
#'
#' Returns `list(name, lineno, everAssigned)` for the first offending read, or
#' `NULL` when every read is safe.
#'
#' @noRd
nmFirstUnboundUse <- function(stmts, bound, everAssigned) {
  hit <- NULL

  checkExpr <- function(node, lineno, bnd) {
    if (!is.null(hit) || is.null(node)) {
      return(invisible(NULL))
    }
    for (nm in nmExprSyms(node)) {
      if (!nm %in% bnd) {
        hit <<- list(
          name = nm, lineno = lineno,
          everAssigned = nm %in% everAssigned
        )
        return(invisible(NULL))
      }
    }
    invisible(NULL)
  }

  ## Returns the set of names bound after running `ss`.
  run <- function(ss, bnd) {
    for (s in ss) {
      if (!is.null(hit)) {
        return(bnd)
      }
      if (s$type == "assign") {
        checkExpr(s$rhs, s$lineno, bnd)
        bnd <- union(bnd, s$lhs)
      } else {
        lineno <- if (is.null(s$lineno)) NA_integer_ else s$lineno
        checkExpr(s$cond, lineno, bnd)
        outs <- list(run(s$then, bnd))
        for (e in s$elifs) {
          checkExpr(e$cond, lineno, bnd)
          outs[[length(outs) + 1L]] <- run(e$stmts, bnd)
        }
        if (!is.null(s$else_)) outs[[length(outs) + 1L]] <- run(s$else_, bnd)
        bnd <- Reduce(union, outs)
      }
    }
    bnd
  }

  run(stmts, bound)
  hit
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
## Which symbol holds each parameter's typical value.
##
## "TV" in front of a parameter name is a convention, not a rule: run7 has
## `TVMAT = THETA(6); MAT = TVMAT * EXP(ETA(5))`, where TVMAT really is MAT's
## typical value, and immediately afterwards `TVD1 = THETA(7);
## D1 = MAT*(1-TVD1)`, where TVD1 is a dimensionless fraction that D1 is
## computed from. Only an assignment of the form `P = <sym> * EXP(ETA(n))` -
## or a bare `P = <sym>` for a parameter without IIV - tells us that <sym> is
## P's typical value, so that is what this records. Anything else is left out
## rather than guessed at.
##
## @noRd
nmTypicalMap <- function(stmts) {
  map <- character(0)

  ## Called on the *folded* tree, where every ETA()/EPS() is already 0 and the
  ## constants around them have collapsed. That reduces the question to "does
  ## this parameter reduce to a single symbol?", which is exactly what a
  ## typical value is - and it covers idioms an exp(ETA(n)) pattern match
  ## cannot, notably IOV: TVCL*EXP(ETA(1) + ETA(2)*(1-OCC) + ETA(3)*OCC) folds
  ## to TVCL, so TVCL is recorded and verifyParamFunction() can compare against
  ## a TVCL column. A parameter computed from others - D1 = MAT*(1-TVD1) -
  ## folds to a binop and gets no entry, which is the right answer.
  tvSym <- function(node) {
    if (node$type == "sym") {
      return(node$name)
    }
    NA_character_
  }

  for (s in stmts) {
    if (s$type != "assign") next
    nm <- tvSym(s$rhs)
    if (!is.na(nm)) map[s$lhs] <- nm else map <- map[names(map) != s$lhs]
  }
  map
}

nmEtaMap <- function(stmts) {
  map <- integer(0)

  ## Count every ETA() in a subtree, so a separable one can be told from an
  ## eta that is scaled or otherwise entangled.
  etasIn <- function(n) {
    if (!is.list(n) || is.null(n$type)) {
      return(integer(0))
    }
    if (identical(n$type, "eta")) {
      return(n$index)
    }
    out <- integer(0)
    ## `arg` too: the parser builds unary minus as a `unop` node, and an eta
    ## under it would otherwise be invisible to the count - so the guard would
    ## pass and a non-separable expression be claimed.
    for (f in c("lhs", "rhs", "cond", "arg")) {
      if (!is.null(n[[f]])) out <- c(out, etasIn(n[[f]]))
    }
    if (!is.null(n$args)) for (a in n$args) out <- c(out, etasIn(a))
    out
  }

  ## Is there an EPS()/ERR() anywhere in this subtree?
  epsIn <- function(n) {
    if (!is.list(n) || is.null(n$type)) {
      return(FALSE)
    }
    if (identical(n$type, "eps")) {
      return(TRUE)
    }
    for (f in c("lhs", "rhs", "cond", "arg")) {
      if (!is.null(n[[f]]) && epsIn(n[[f]])) {
        return(TRUE)
      }
    }
    if (!is.null(n$args)) {
      for (a in n$args) {
        if (epsIn(a)) {
          return(TRUE)
        }
      }
    }
    FALSE
  }

  ## Names whose value carries randomness, transitively. An eta can reach an
  ## exponent through a symbol - the IOV idiom assigns ETA() to a name in an
  ## occasion block and adds that name inside EXP() - and a syntactic count of
  ## ETA() nodes cannot see it. Without this, EXP(ETA(1) + IOV) looked like the
  ## MU-referenced idiom and the entry was claimed separable, after which
  ## dividing the tabled value by exp(eta1) leaves exp(IOV) behind.
  etaDerived <- local({
    binds <- list()
    walk <- function(ss) {
      for (st in ss) {
        if (identical(st$type, "assign")) {
          binds[[length(binds) + 1L]] <<- list(lhs = st$lhs, rhs = st$rhs)
        } else {
          walk(st$then)
          for (e in st$elifs) walk(e$stmts)
          if (!is.null(st$else_)) walk(st$else_)
        }
      }
    }
    walk(stmts)
    seeds <- character(0)
    repeat {
      before <- seeds
      for (b in binds) {
        if (b$lhs %in% seeds) next
        if (length(etasIn(b$rhs)) > 0 ||
          any(nmExprSyms(b$rhs) %in% seeds)) {
          seeds <- c(seeds, b$lhs)
        }
      }
      if (identical(seeds, before)) break
    }
    seeds
  })

  ## Top-level additive terms. Only `+`: under `-` the eta enters with the
  ## wrong sign and dividing by exp(ETA) would not undo it.
  addTerms <- function(n) {
    if (is.list(n) && identical(n$type, "binop") && identical(n$op, "+")) {
      c(addTerms(n$lhs), addTerms(n$rhs))
    } else {
      list(n)
    }
  }

  ## exp(...) whose argument carries exactly one ETA, as a bare additive term.
  ## That covers both idioms in use:
  ##   CL = TVCL * EXP(ETA(3))          classic
  ##   CL = EXP(MU_6 + ETA(6))          MU-referenced, as IMP and SAEM write it
  ## and both are separable for the same reason - EXP(a + eta) is
  ## EXP(a) * EXP(eta) - so the tabled individual value divided by exp(eta)
  ## gives the typical value either way.
  fromExp <- function(e) {
    if (!(is.list(e) && identical(e$type, "call") && identical(e$fn, "exp") &&
      length(e$args) == 1L)) {
      return(NA_integer_)
    }
    if (length(etasIn(e$args[[1]])) != 1L) {
      return(NA_integer_)
    }
    ## An EPS() in the exponent is a second source of randomness, and unlike
    ## the IOV-through-a-symbol case it is not a symbol, so the sibling check
    ## below cannot see it. EXP(ETA(1) + EPS(1)) counts one ETA() and would
    ## otherwise be claimed - after which dividing by exp(eta1) leaves
    ## exp(EPS(1)) behind.
    if (epsIn(e$args[[1]])) {
      return(NA_integer_)
    }
    terms <- addTerms(e$args[[1]])
    bare <- Filter(function(t) identical(t$type, "eta"), terms)
    if (length(bare) != 1L) {
      return(NA_integer_)
    }
    ## No sibling of the bare eta may itself carry randomness.
    others <- Filter(function(t) !identical(t$type, "eta"), terms)
    for (t in others) {
      if (any(nmExprSyms(t) %in% etaDerived)) {
        return(NA_integer_)
      }
    }
    bare[[1]]$index
  }

  etaInExp <- function(node) {
    cands <- list(node)
    n <- node
    while (is.list(n) && identical(n$type, "binop") && identical(n$op, "*")) {
      cands <- c(cands, list(n$lhs, n$rhs))
      n <- n$lhs
    }
    for (c0 in cands) {
      idx <- fromExp(c0)
      if (!is.na(idx)) {
        return(idx)
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
      call = lapply(node$args, walkExpr),
      unop = walkExpr(node$arg),
      binop = {
        walkExpr(node$lhs)
        walkExpr(node$rhs)
      },
      NULL
    )
    invisible(NULL)
  }
  walkStmts <- function(ss) {
    for (s in ss) {
      if (s$type == "assign") {
        walkExpr(s$rhs)
      } else {
        walkExpr(s$cond)
        walkStmts(s$then)
        for (e in s$elifs) {
          walkExpr(e$cond)
          walkStmts(e$stmts)
        }
        if (!is.null(s$else_)) walkStmts(s$else_)
      }
    }
  }
  walkStmts(stmts)
  mx
}


## Fold OMEGA(i,j) / SIGMA(i,j) to the value in the .ext file.
##
## They are constants once the model is estimated, so a parameter function can
## carry the number. Without an .ext there is nothing to resolve them from, and
## guessing is not on the table.
##
## @noRd
nmResolveMatrixRefs <- function(stmts, ext, modFile, block = "$PK") {
  needed <- list()

  walk <- function(node) {
    if (!is.list(node) || is.null(node$type)) {
      return(node)
    }
    if (identical(node$type, "nmmatrix")) {
      needed[[length(needed) + 1L]] <<- node
      return(node)
    }
    for (f in c("lhs", "rhs", "cond")) {
      if (!is.null(node[[f]])) node[[f]] <- walk(node[[f]])
    }
    if (!is.null(node$args)) node$args <- lapply(node$args, walk)
    node
  }
  scan <- function(ss) {
    for (st in ss) {
      if (identical(st$type, "assign")) {
        walk(st$rhs)
      } else {
        walk(st$cond)
        scan(st$then)
        for (e in st$elifs) {
          walk(e$cond)
          scan(e$stmts)
        }
        if (!is.null(st$else_)) scan(st$else_)
      }
    }
  }
  scan(stmts)
  if (length(needed) == 0L) {
    return(stmts)
  }

  if (is.null(ext)) {
    nm <- unique(vapply(needed, function(n) {
      sprintf("%s(%d,%d)", n$mat, n$i, n$j)
    }, ""))
    stop("The ", block, " block of ", basename(modFile), " reads ",
      paste(nm, collapse = ", "),
      ", whose value is in the model's .ext file. Supply extFile.",
      call. = FALSE
    )
  }

  final <- ext[ext$ITERATION == -1000000000, , drop = FALSE]
  if (nrow(final) == 0L) {
    stop("No final estimates in the .ext file, so ",
      "OMEGA()/SIGMA() cannot be resolved.",
      call. = FALSE
    )
  }
  lookup <- function(node) {
    ## The .ext carries the lower triangle, and R mangles "OMEGA(2,2)" to
    ## "OMEGA.2.2." on read; the matrix is symmetric, so try both orders.
    cand <- c(
      sprintf("%s.%d.%d.", node$mat, node$i, node$j),
      sprintf("%s.%d.%d.", node$mat, node$j, node$i)
    )
    hit <- cand[cand %in% names(final)]
    if (length(hit) == 0L) {
      stop(sprintf(
        "%s(%d,%d) is not in the .ext file of %s.",
        node$mat, node$i, node$j, basename(modFile)
      ), call. = FALSE)
    }
    as.numeric(final[[hit[1]]][1])
  }
  subst <- function(node) {
    if (!is.list(node) || is.null(node$type)) {
      return(node)
    }
    if (identical(node$type, "nmmatrix")) {
      return(list(type = "num", value = lookup(node)))
    }
    for (f in c("lhs", "rhs", "cond")) {
      if (!is.null(node[[f]])) node[[f]] <- subst(node[[f]])
    }
    if (!is.null(node$args)) node$args <- lapply(node$args, subst)
    node
  }
  fix <- function(ss) {
    lapply(ss, function(st) {
      if (identical(st$type, "assign")) {
        st$rhs <- subst(st$rhs)
      } else {
        st$cond <- subst(st$cond)
        st$then <- fix(st$then)
        st$elifs <- lapply(st$elifs, function(e) {
          e$cond <- subst(e$cond)
          e$stmts <- fix(e$stmts)
          e
        })
        if (!is.null(st$else_)) st$else_ <- fix(st$else_)
      }
      st
    })
  }
  fix(stmts)
}

## Keep only the statements the requested parameters depend on.
##
## $PK is all parameter definitions, so emitting the whole block is merely
## verbose - the extra entries are correct, just unwanted. Pruning is still
## worth doing: it shortens the source the caller has to read against the
## control stream, and it shrinks the covariate set to those the requested
## parameters actually use, so no reference has to be justified for a covariate
## that cannot reach the answer.
##
## Dependencies in $PK only ever run backwards: a statement can be affected by
## earlier statements and never by later ones. One reverse pass therefore
## suffices. A variable is never taken out of the needed set once it is in it,
## so a re-assignment chain - `TVCL = THETA(4)*CLCOV1` and then
## `TVCL = CLCOV*TVCL` - keeps both of its links.
##
## An IF block is kept whole when anything inside it is needed. Pruning within
## a block would leave empty branches for no gain; carrying a few extra
## assignments is the cheaper mistake.
##
## @noRd
nmPruneToParameters <- function(stmts, parameters) {
  need <- parameters

  blockAssigns <- function(st) {
    if (identical(st$type, "assign")) {
      return(st$lhs)
    }
    out <- unlist(lapply(st$then, blockAssigns))
    for (e in st$elifs) out <- c(out, unlist(lapply(e$stmts, blockAssigns)))
    if (!is.null(st$else_)) out <- c(out, unlist(lapply(st$else_, blockAssigns)))
    out
  }
  blockSyms <- function(st) {
    if (identical(st$type, "assign")) {
      return(nmExprSyms(st$rhs))
    }
    out <- c(nmExprSyms(st$cond), unlist(lapply(st$then, blockSyms)))
    for (e in st$elifs) {
      out <- c(out, nmExprSyms(e$cond), unlist(lapply(e$stmts, blockSyms)))
    }
    if (!is.null(st$else_)) out <- c(out, unlist(lapply(st$else_, blockSyms)))
    out
  }

  keptRev <- list()
  for (i in rev(seq_along(stmts))) {
    st <- stmts[[i]]
    if (!any(blockAssigns(st) %in% need)) next
    need <- unique(c(need, blockSyms(st)))
    keptRev[[length(keptRev) + 1L]] <- st
  }
  rev(keptRev)
}

## Every name assigned anywhere in a statement tree.
## @noRd
nmAssignedNames <- function(stmts) {
  out <- character(0)
  walk <- function(ss) {
    for (st in ss) {
      if (identical(st$type, "assign")) {
        out <<- c(out, st$lhs)
      } else {
        walk(st$then)
        for (e in st$elifs) walk(e$stmts)
        if (!is.null(st$else_)) walk(st$else_)
      }
    }
  }
  walk(stmts)
  unique(out)
}

## Every symbol read anywhere in a statement (conditions included).
## @noRd
nmStmtSyms <- function(st) {
  if (identical(st$type, "assign")) {
    return(nmExprSyms(st$rhs))
  }
  out <- nmExprSyms(st$cond)
  for (x in st$then) out <- c(out, nmStmtSyms(x))
  for (e in st$elifs) {
    out <- c(out, nmExprSyms(e$cond))
    for (x in e$stmts) out <- c(out, nmStmtSyms(x))
  }
  if (!is.null(st$else_)) for (x in st$else_) out <- c(out, nmStmtSyms(x))
  unique(out)
}
