## Unit tests for the internal NONMEM parsing machinery. The public behaviour is
## covered in test-createParamFunction.R; these pin the pieces individually.

## Parse one expression string and render it as R.
r_of <- function(txt) {
  p <- nmParser(nmLex(txt, 1L, "unit.mod"), 1L, "unit.mod")
  nmDeparse(nmParseExpr(p))
}

## Parse a whole $PK body given as a character vector.
stmts_of <- function(lines, file = "unit.mod") {
  mod <- data.frame(lineno = seq_along(lines), code = lines,
                    comment = "", stringsAsFactors = FALSE)
  nmParseStatements(mod, file)
}

test_that("the lexer handles numbers, symbols and dot operators", {
  expect_equal(r_of("1.0000E+00"), "1")
  expect_equal(r_of("1E-12"), "1e-12")
  expect_equal(r_of("0.75"), "0.75")
  expect_equal(r_of("THETA(4)"), "thetas[4]")
  expect_equal(r_of("WT"), "WT")
  # A Fortran D exponent is a number, not a dot operator.
  expect_equal(r_of("1.5D+02"), "150")
})

test_that("dot operators translate to their R equivalents", {
  expect_equal(r_of("SEX.EQ.1"), "SEX == 1")
  expect_equal(r_of("SEX.NE.1"), "SEX != 1")
  expect_equal(r_of("A.GE.1"),   "A >= 1")
  expect_equal(r_of("A.LE.1"),   "A <= 1")
  expect_equal(r_of("A.GT.1"),   "A > 1")
  expect_equal(r_of("A.LT.1"),   "A < 1")
  # .EQN. must win over .EQ. despite sharing a prefix.
  expect_equal(r_of("A.EQN.2"),  "A == 2")
  expect_equal(r_of("A.EQ.1.AND.B.EQ.2"), "A == 1 & B == 2")
  expect_equal(r_of("A.EQ.1.OR.B.EQ.2"),  "A == 1 | B == 2")
})

test_that("arithmetic keeps NONMEM's meaning", {
  expect_equal(r_of("(WT/75)**THETA(2)"), "(WT / 75)^thetas[2]")
  expect_equal(r_of("A+B*C"),   "A + B * C")
  expect_equal(r_of("(A+B)*C"), "(A + B) * C")
  expect_equal(r_of("A-B-C"),   "A - B - C")
  expect_equal(r_of("A-(B-C)"), "A - (B - C)")
  expect_equal(r_of("A/B/C"),   "A / B / C")
  expect_equal(r_of("A/(B/C)"), "A / (B / C)")
  # ** is right-associative in Fortran, as ^ is in R. The values must differ:
  # with 2,2,2 both associativities give 16, so the assertion could not fail.
  expect_equal(eval(parse(text = chartr("ABC", "232", r_of("A**B**C")))), 2^(3^2))
  # unary minus binds looser than a power: -2**2 is -4
  expect_equal(eval(parse(text = r_of("-2**2"))), -4)
})

test_that("MOD() keeps Fortran's meaning, not R's `%%`", {
  # `%%` is wrong twice over, so neither half may be reintroduced.

  # 1. Precedence. MOD() is a call, but `%any%` binds tighter than * and / in
  #    R, so emitting `%%` turned MOD(A*B, C) into A * (B %% C).
  expect_equal(eval(parse(text = chartr("ABC", "732", r_of("MOD(A*B,C)")))), 1)
  expect_equal(eval(parse(text = chartr("ABC", "732", r_of("MOD(A/B,C)")))),
               (7 / 3) - 2 * trunc((7 / 3) / 2))
  # and MOD() nested inside an outer operator stays intact
  expect_equal(eval(parse(text = chartr("ABC", "732", r_of("C*MOD(A,B)")))), 2)

  # 2. Sign. Fortran MOD() truncates towards zero, `%%` floors, so they differ
  #    for a negative first argument: MOD(-7, 3) is -1, but -7 %% 3 is 2.
  expect_equal(eval(parse(text = r_of("MOD(-7,3)"))), -1)
  expect_equal(eval(parse(text = r_of("MOD(7,-3)"))), 1)
  expect_equal(eval(parse(text = r_of("MOD(-7,-3)"))), -1)

  expect_false(grepl("%%", r_of("MOD(A,B)"), fixed = TRUE))
  expect_error(r_of("MOD(A)"), "exactly two arguments")
})

test_that("intrinsic functions map to R", {
  expect_equal(r_of("EXP(X)"),   "exp(X)")
  expect_equal(r_of("LOG(X)"),   "log(X)")
  expect_equal(r_of("SQRT(X)"),  "sqrt(X)")
  expect_equal(r_of("LOG10(X)"), "log10(X)")
  expect_equal(r_of("MAX(A,B)"), "max(A, B)")
})

test_that("constructs needing an ODE solution or unknown syntax are refused", {
  expect_error(r_of("A(2)*1000"), "compartment amount")
  expect_error(r_of("EPS(1)"),    "cannot appear in \\$PK")
  expect_error(r_of("FOO(X)"),    "unsupported function")
  expect_error(r_of("A @ B"),     "unrecognised character")
  expect_error(r_of("(A+B"),      "unbalanced parentheses")
})

test_that("assignments, one-line IFs and IF blocks all parse", {
  s <- stmts_of(c("CL = THETA(1)"))
  expect_length(s, 1)
  expect_equal(s[[1]]$type, "assign")
  expect_equal(s[[1]]$lhs, "CL")

  s <- stmts_of(c("IF(SEX.EQ.1) CL = THETA(1)"))
  expect_equal(s[[1]]$type, "if")
  expect_true(s[[1]]$oneline)
  expect_length(s[[1]]$then, 1)

  s <- stmts_of(c("IF (SEX.EQ.1) THEN", "  CL = THETA(1)", "ELSE",
                  "  CL = THETA(2)", "END IF"))
  expect_equal(s[[1]]$type, "if")
  expect_null(s[[1]]$oneline)
  expect_length(s[[1]]$then, 1)
  expect_length(s[[1]]$else_, 1)

  s <- stmts_of(c("IF (A.EQ.1) THEN", "  X = 1", "ELSE IF (A.EQ.2) THEN",
                  "  X = 2", "ELSE", "  X = 3", "ENDIF"))
  expect_length(s[[1]]$elifs, 1)
  expect_length(s[[1]]$else_, 1)
})

test_that("malformed or unsupported statements are refused with a line number", {
  expect_error(stmts_of(c("CL = THETA(1)", "DO WHILE (X.LT.2)")),
               "unit\\.mod:2")
  expect_error(stmts_of(c("CALL SUBR(X)")), "CALL")
  expect_error(stmts_of(c("IF (A.EQ.1) THEN", "  X = 1")), "never closed")
  expect_error(stmts_of(c('"  CALL FOO')), "verbatim FORTRAN")
  expect_error(stmts_of(c("JUST A PHRASE")), "not an assignment")
})

test_that("ETA() folds to typical values and artefacts are simplified away", {
  s <- nmSimplifyStmts(stmts_of(c("CL = TVCL*EXP(ETA(3))")))
  expect_equal(nmDeparse(s[[1]]$rhs), "TVCL")
  expect_true(s[[1]]$hadEta)

  s <- nmSimplifyStmts(stmts_of(c("Y = X+ETA(1)")))
  expect_equal(nmDeparse(s[[1]]$rhs), "X")

  # Folding must not touch arithmetic that is not an identity
  s <- nmSimplifyStmts(stmts_of(c("Y = X*2")))
  expect_equal(nmDeparse(s[[1]]$rhs), "X * 2")
})

test_that("symbols are split into assigned and used", {
  s <- stmts_of(c("CLWT = (WT/75)**THETA(2)", "TVCL = THETA(4)*CLWT"))
  sym <- nmSymbols(s)
  expect_equal(sym$assigned, c("CLWT", "TVCL"))
  expect_true(all(c("WT", "CLWT") %in% sym$used))
  expect_false("TVCL" %in% sym$used)
})

test_that("the highest THETA index is found, including inside IF blocks", {
  s <- stmts_of(c("A = THETA(1)", "IF (X.EQ.1) B = THETA(9)"))
  expect_equal(nmMaxTheta(s), 9)
})

test_that("the exponential-IIV pattern is recorded, other forms are not", {
  s <- stmts_of(c("CL = TVCL*EXP(ETA(3))", "V = EXP(ETA(4))*TVV",
                  "MAT = TVMAT+ETA(5)"))
  m <- nmEtaMap(s)
  expect_equal(m[["CL"]], 3L)
  expect_equal(m[["V"]],  4L)
  expect_false("MAT" %in% names(m))
})

test_that("$THETA records are counted by theta, not by line", {
  cnt <- function(x) {
    nmCountThetas(data.frame(lineno = seq_along(x), code = x, comment = "",
                             stringsAsFactors = FALSE))
  }
  expect_equal(cnt("$THETA 1 2 3"), 3)                      # not 1
  expect_equal(cnt("$THETA (0,11.9)"), 1)
  expect_equal(cnt(c("$THETA (0,1) (0,2)", "$THETA 3")), 3)
  expect_equal(cnt("$THETA 1 FIX"), 1)
  expect_equal(cnt("$THETA (0,1)x3"), 3)
  expect_equal(cnt("$THETA 0.1 ; a comment with 999 in it"), 1)
  expect_equal(cnt("$OMEGA 1"), 0)
})

test_that("$INPUT names are read, with DROP columns and synonyms handled", {
  inp <- function(x) {
    nmInputNames(data.frame(lineno = seq_along(x), code = x, comment = "",
                            stringsAsFactors = FALSE))
  }
  expect_equal(inp("$INPUT ID DV WT"), c("ID", "DV", "WT"))
  expect_equal(inp("$INP ID DV"), c("ID", "DV"))
  expect_false("JUNK" %in% inp("$INPUT ID JUNK=DROP WT"))
  expect_true(all(c("CONC", "DV") %in% inp("$INPUT ID CONC=DV")))
  # continuation onto a second line
  expect_equal(inp(c("$INPUT ID DV", "       WT AGE")), c("ID", "DV", "WT", "AGE"))
})

test_that("comments are separated from code but kept available", {
  f <- withr::local_tempfile(fileext = ".mod")
  writeLines(c("$PK", "CL = THETA(1) ; the clearance"), f)
  mod <- nmReadModel(f)
  expect_equal(trimws(mod$code[2]), "CL = THETA(1)")
  expect_equal(mod$comment[2], "the clearance")
})

test_that("an explicit missing-value branch is used as the reference (rule 1)", {
  s <- stmts_of(c("IF(WT.EQ.-99) WT = 70", "CLWT = WT*THETA(1)"))
  r <- nmCovRef(s, "WT", missVal = -99)
  expect_equal(r$WT$value, 70)
  expect_true(r$WT$confident)
  expect_match(r$WT$source, "explicit missing-value handling")
})

test_that("reference rules are tried in order of confidence", {
  # An explicit branch beats a normalisation constant.
  s <- stmts_of(c("IF(WT.EQ.-99) WT = 70", "CLWT = (WT/75)**THETA(1)"))
  expect_equal(nmCovRef(s, "WT", -99)$WT$value, 70)

  # With no explicit branch, the normalisation constant is used.
  s <- stmts_of(c("CLWT = (WT/75)**THETA(1)"))
  expect_equal(nmCovRef(s, "WT", -99)$WT$value, 75)

  # A subtractive centering constant works too.
  s <- stmts_of(c("RF = (AGE-50)*THETA(1)"))
  expect_equal(nmCovRef(s, "AGE", -99)$AGE$value, 50)

  # A covariate entering linearly has no derivable reference.
  s <- stmts_of(c("EFF = THETA(1)*EXPO"))
  expect_length(nmCovRef(s, "EXPO", -99), 0)
})

test_that("the identity-value branch is recognised without a marker (rule 2b)", {
  s <- stmts_of(c("IF(FORM.EQ.1) FRELFORM = 1",
                  "IF(FORM.EQ.0) FRELFORM = 1+THETA(1)"))
  r <- nmCovRef(s, "FORM", -99)
  expect_equal(r$FORM$value, 1)
  expect_true(r$FORM$confident)
  expect_match(r$FORM$source, "identity value")

  # An additive model has 0 as its identity.
  s <- stmts_of(c("IF(SMOK.EQ.0) EFF = 0", "IF(SMOK.EQ.1) EFF = THETA(1)"))
  expect_equal(nmCovRef(s, "SMOK", -99)$SMOK$value, 0)
})

test_that("equality tests are found whichever way round they are written", {
  expect_equal(nmEqualityTest(nmParseExpr(nmParser(nmLex("SEX.EQ.2", 1L, "u"), 1L, "u")), "SEX"), 2)
  expect_equal(nmEqualityTest(nmParseExpr(nmParser(nmLex("2.EQ.SEX", 1L, "u"), 1L, "u")), "SEX"), 2)
  # and through an .AND. chain, from either side
  cond <- nmParseExpr(nmParser(nmLex("STUDY.EQ.1.AND.SEX.EQ.2", 1L, "u"), 1L, "u"))
  expect_equal(nmEqualityTest(cond, "SEX"), 2)
  expect_equal(nmEqualityTest(cond, "STUDY"), 1)
  expect_null(nmEqualityTest(cond, "AGE"))
  # a non-equality comparison is not a reference
  expect_null(nmEqualityTest(nmParseExpr(nmParser(nmLex("WT.GT.70", 1L, "u"), 1L, "u")), "WT"))
  expect_null(nmEqualityTest(NULL, "WT"))
})

test_that("the normalisation constant is found inside a nested expression", {
  s <- stmts_of(c("CL = THETA(1)*EXP(LOG(WT/75)*THETA(2))"))
  expect_equal(nmCovRef(s, "WT", -99)$WT$value, 75)
  s <- stmts_of(c("CL = -(WT/80)*THETA(1)"))
  expect_equal(nmCovRef(s, "WT", -99)$WT$value, 80)
})

test_that("references are derived from inside IF blocks too", {
  s <- stmts_of(c("IF (STUDY.EQ.1) THEN", "  CLWT = (WT/75)**THETA(1)",
                  "ELSE", "  CLWT = 1", "END IF"))
  expect_equal(nmCovRef(s, "WT", -99)$WT$value, 75)
})

test_that("nmRecord finds a record, its end, and copes with absence", {
  mod <- nmReadModel(system.file("extdata", "SimVal/run7.mod", package = "PMXForest"))
  pk <- nmRecord(mod, "\\$PK\\b")
  expect_gt(nrow(pk), 10)
  # the record stops before the next $ record
  expect_false(any(grepl("^\\s*\\$ERROR", pk$code)))
  # matching is case-insensitive
  expect_equal(nrow(nmRecord(mod, "\\$pk\\b")), nrow(pk))
  # an absent record gives zero rows rather than an error
  expect_equal(nrow(nmRecord(mod, "\\$NOSUCH\\b")), 0)
  expect_equal(nmInputNames(nmRecord(mod, "\\$NOSUCH\\b")), character(0))
})

test_that("a missing model file is reported clearly", {
  expect_error(nmReadModel("no-such-file.mod"), "Model file not found")
})

test_that("emitted IF blocks keep their else and else-if branches", {
  s <- nmSimplifyStmts(stmts_of(c(
    "IF (A.EQ.1) THEN", "  X = 1", "ELSE IF (A.EQ.2) THEN", "  X = 2",
    "ELSE", "  X = 3", "END IF"
  )))
  code <- nmEmitStmts(s, "  ", "u.mod")
  expect_true(any(grepl("^\\s*if \\(A == 1\\) \\{", code)))
  expect_true(any(grepl("^\\s*\\} else if \\(A == 2\\) \\{", code)))
  expect_true(any(grepl("^\\s*\\} else \\{", code)))
  # the branch bodies round-trip through R
  f <- eval(parse(text = c("function(A) {", code, "X }")))
  expect_equal(f(1), 1); expect_equal(f(2), 2); expect_equal(f(3), 3)
})

test_that("nmHasEta sees an ETA at any depth", {
  has <- function(txt) nmHasEta(stmts_of(paste0("Y = ", txt))[[1]]$rhs)
  expect_true(has("EXP(ETA(1))"))
  expect_true(has("A*B+EXP(ETA(1))"))
  expect_true(has("-ETA(1)"))
  expect_false(has("A*B+THETA(1)"))
})

test_that("trailing tokens and malformed statements are refused", {
  expect_error(stmts_of("IF (A.EQ.1 2) X = 1"), "unbalanced|trailing")
  expect_error(stmts_of("X = 1 2"), "trailing tokens after the assignment")
  expect_error(stmts_of("IF (A.EQ.1) X = 1 2"), "trailing tokens after the assignment")
  expect_error(stmts_of("IF (A.EQ.1)"), "without THEN")
  expect_error(stmts_of("X(1) = 2"), "assignment to a plain variable")
  # an unsupported construct is refused wherever it appears, including in a block
  expect_error(stmts_of(c("IF (A.EQ.1) THEN", "  DO I=1,2", "END IF")),
               "unit\\.mod:2")
})

test_that("constant folding covers every identity it claims", {
  fold <- function(txt) nmDeparse(nmSimplify(
    nmParseExpr(nmParser(nmLex(txt, 1L, "u"), 1L, "u"))))
  expect_equal(fold("0+X"),        "X")
  expect_equal(fold("X+0"),        "X")
  expect_equal(fold("X-0"),        "X")
  expect_equal(fold("1*X"),        "X")
  expect_equal(fold("X*1"),        "X")
  expect_equal(fold("X/1"),        "X")
  expect_equal(fold("X**1"),       "X")
  expect_equal(fold("LOG(1)"),     "0")
  expect_equal(fold("EXP(0)"),     "1")
  # A negated literal folds to a number. "-2" alone cannot show this: it
  # deparses as "-2" whether or not the fold ran. These need the fold, because
  # the identity tests match a `num` node and not a negated one.
  expect_equal(fold("EXP(-0)"),    "1")
  expect_equal(fold("X+-0"),       "X")
  expect_equal(fold("-2"),         "-2")
  # nothing that is not an identity is touched
  expect_equal(fold("X/2"),        "X / 2")
  expect_equal(fold("X**2"),       "X^2")
  expect_equal(fold("LOG(X)"),     "log(X)")
  expect_equal(fold("-X"),         "-X")
})

test_that("statement walkers descend into every branch of an IF block", {
  s <- stmts_of(c("IF (A.EQ.1) THEN", "  X = -THETA(1)",
                  "ELSE IF (B.EQ.2) THEN", "  X = THETA(7)",
                  "ELSE", "  X = THETA(3)*C", "END IF"))
  sym <- nmSymbols(s)
  expect_setequal(sym$assigned, "X")
  # symbols from the condition, the else-if and the else branch are all seen
  expect_true(all(c("A", "B", "C") %in% sym$used))
  expect_equal(nmMaxTheta(s), 7)
})

test_that("NA and extreme numeric literals format sensibly", {
  expect_equal(nmFormatNum(NA_real_), "NA")
  expect_equal(nmFormatNum(75), "75")
  expect_equal(nmFormatNum(0.75), "0.75")
  expect_equal(as.numeric(nmFormatNum(1 / 3)), 1 / 3)
  expect_equal(as.numeric(nmFormatNum(1e-12)), 1e-12)
})

test_that("deparse renders ETA and unary nodes, and rejects an unknown node", {
  eta <- list(type = "eta", index = 2L)
  expect_equal(nmDeparse(eta), "0")
  expect_equal(nmDeparse(eta, etaValue = "etas[2]"), "etas[2]")
  expect_equal(nmDeparse(list(type = "unop", op = "!",
                              arg = list(type = "sym", name = "A"))), "!A")
  expect_error(nmDeparse(list(type = "nonsense")), "cannot deparse node type")
})
