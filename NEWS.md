# PMXForest 1.3.0

This release adds a way to generate the parameter function from the control
stream instead of writing it by hand, convenience functions for the other
`getForestDFSCM()` inputs, and a set of fixes - two of which change numbers in
plots you have already made. Start with **Changes to existing behaviour**.

## Changes to existing behaviour

* **Forest plots built from SIR files will change.** `getSamples()`'s
  SIR-detection check looked for a `samples_order` column that PsN does not
  produce (the column is `sample_order`), so SIR files were routed through the
  bootstrap path. That path filters on `ofv != 0`, and therefore returned the
  full SIR *proposal* distribution rather than the importance-resampled
  parameter vectors (`resamples == 1`). It now returns the resampled vectors,
  with the final estimates prepended as the first row, as for every other input
  type. For a well-converged SIR run the confidence intervals shift by a few
  percent, not systematically in one direction; for a poorly initialised run the
  change can be larger. Bootstrap and `.cov` inputs are unaffected - those files
  carry neither column, so the check was always false for them. `getSamples()`
  also gains `quiet` (default `FALSE`), which prints a message when a SIR file
  is detected.

* **The statistics table shows more digits on the relative scale.** With
  neither `sigdigits` nor `decimals` set - the default - the table now uses 2
  decimal places on the relative scale, and 2 significant digits on the absolute
  scale as before. Relative covariate effects cluster around 1, where 2
  significant digits discarded the second decimal, so a row that read `1.0` now
  reads `1.03`. Only the printed numbers change; points, intervals and axis
  labels are unaffected. Calls passing `sigdigits` explicitly are unchanged.
  `forestPlot()` gains `decimals` alongside `sigdigits`; supply one or the other.

* **An explicit `catRef` (or `refLevels`) now sets the reference *value* as well
  as the one-hot column names.** Previously `refLevels` only chose which level
  was dropped when encoding, while the reference state was always the most
  common level in the data - so `refLevels = list(RACE = 1)` produced columns
  keyed on level 1 but a background sitting at the modal level. The two are now
  one setting. Calls that pass neither are unaffected: the encoding still drops
  the lowest level and the reference still comes from the mode.

* **`refLevels` is deprecated in favour of `catRef`.** `getCovStats()`,
  `setupDfCovs()`, `setupDfRefRow()` and `setupCovExpressionsList()` all take
  `catRef`, symmetric with `contRef`. `refLevels` still works and forwards to
  `catRef` with a deprecation warning; supplying both is an error.

* **`getCovStats()` sorts binary covariate levels.** They were returned in order
  of appearance in the data, so row order in a Forest plot depended on which
  subject came first. Now always `c(0, 1)`.

* **`getCovStats()` refuses a covariate with no non-missing value** instead of
  returning an empty result. `createInputForestData()` emits no rows for such a
  covariate, so it silently vanished from the plot. `setupDfCovs()` inherits the
  error. `refValues()` already behaved this way.

* **The statistics panel label gains its separating space.** `forestPlot()`'s
  `statisticsLabel` defaulted to `"Statistics:"` and is prepended verbatim, so
  every panel was titled `Statistics:CL (L/h)`. The default is now
  `"Statistics: "`. A label you pass yourself is still used exactly as given, so
  include the space you want.

## Generating the parameter function from the control stream

The parameter function used to be written by hand, duplicating algebra that
already exists in `$PK` - the easiest place in the workflow to introduce an
error no test would catch.

* **`createParamFunction()`** translates a `$PK` block into R source for a
  `function(thetas, df, ...)`. It returns **text**, for you to read and check;
  nothing is evaluated on your behalf. Every `ETA(n)` is set to 0, so the
  function returns typical values. All missing-covariate handling is hoisted
  into one annotated preamble, each reference value taken from the control
  stream itself by four rules in decreasing order of confidence: explicit
  `IF(WT.EQ.-99)` handling; the branch PsN's scm marks `; Most common`; the
  branch assigning the identity value; the normalisation constant in `(WT/75)`
  or `(AGE-50)`. A weaker fifth rule proposes the level no `IF()` tests and
  warns. A covariate matching no rule is an error, and can be supplied through
  `covRef`.

  The parser accepts assignments, `IF` constructs and closed-form arithmetic,
  and refuses what it cannot translate faithfully - `$DES`, compartment amounts
  `A(n)`, verbatim FORTRAN, `DO`, `CALL` - naming the file and line rather than
  guessing. It also refuses a symbol read before anything assigns it: NONMEM
  neither initialises `$PK` variables nor clears them between data records, so
  such a model reads whatever the previous subject left behind, and a parameter
  function evaluated one row at a time cannot reproduce that.

* **`secondary`** attaches quantities `$PK` does not contain - AUC, `Cmax`, an
  event probability. Each entry is a line of R code
  (`secondary = list(AUC = "df$DOSE / CL")`) or the path to an `.R` file of
  arbitrary code, including a `deSolve` or `mrgsolve` simulation, whose text is
  **inlined** so the result stays self-contained. Entries are spliced inside
  `local({ ... })`, so each sees `thetas`, `df` and every structural parameter
  by name while its own temporaries do not leak, and are emitted in order so a
  later one may use an earlier one. An entry may also be
  `list(source = <string>, dose = 100, tau = 12)`, emitting the named constants
  immediately before the code, so a reusable secondary file can be parametrised
  at the call site. Names are appended to `functionListName` and recorded in
  `secondaryNames`, with `primaryNames` holding the `$PK` parameters alone.

* **`verifyParamFunction()`** evaluates the generated function over a NONMEM
  `$TABLE` and compares it with the values NONMEM wrote, turning "do I trust
  this translation?" into a pass or a fail. Table output is written after
  `IGNORE`/`ACCEPT`, so taking both parameters and covariates from the table
  removes any need to reproduce the `$DATA` filtering. Where `$PK` shows
  exponential IIV the tabled value is divided by `exp(ETA(n))` to recover the
  typical value. It returns a single `TRUE`/`FALSE`, usable directly in an `if`,
  with the per-parameter table in `attr(., "checks")`. Only the exponential-IIV
  idiom `P = <expr> * EXP(ETA(n))` is recognised, so a MU-referenced model
  yields an empty `etaMap` and nothing to compare.

* **The parser is reusable.** `nmParsePK()` returns the parsed `$PK` tree with
  `ETA()` intact, the covariates and their references, the THETA count and the
  `etaMap`, without emitting source. `nmDeparse()`, `nmFormatNum()` and
  `nmResolveSecondary()` are exported alongside it, so other packages can build
  their own emitters on the same tested parser; PMXFrem uses it for FREM
  parameter functions.

## Building the other inputs from the data

* **`setupDfCovs()`** composes `getCovStats()` and `createInputForestData()`
  into one call: continuous covariates become their 5th and 95th percentiles,
  categorical ones one row per level. `conditionalCovs` names the covariates the
  others are conditioned on - fed state rather than fasted, patients rather than
  healthy volunteers. They get their own rows and sit at their reference on every
  other row, rather than at `missVal`. `useMissVal` toggles whether inactive
  cells hold `missVal` or a computed baseline. Statistics are computed
  on deduplicated data - one record per `idVar` - so longitudinal data cannot
  skew them.

* **`setupDfRefRow()`** builds the reference row `getForestDF*()` expects,
  mapping computed baseline values onto the `dfCovs` geometry. `singleRef = TRUE`
  gives one row; `singleRef = FALSE` gives a full overlay preserving background
  missingness.

* **`setupCovExpressionsList()`** is the empirical counterpart of
  `setupDfCovs()`, producing the `covExpressionsList` and matching
  `cdfCovsNames` that `getForestDFemp()` consumes. Continuous covariates split
  at the `probs` tails or the median (`contSplit = "median"`); multi-level
  categoricals emit one row per level. Every generated expression is checked
  against the data and the call stops if one selects fewer than `minSubjects`
  subjects (default 10). It applies the same level and quantile rules as
  `getCovStats()` but implements them separately, because the two return
  different shapes - a change to one does not follow into the other.

* **`filterByModel()`** applies the `IGNORE` and `ACCEPT` statements in a
  control stream's `$DATA` record, returning the rows the model actually used.
  Summarising the raw analysis file instead describes a population the model
  never saw: on the bundled `run7`, the file holds 964 subjects and the model
  reads 754, and CRCL's 5th percentile moves from 74.9 to 77.7.

  Columns are matched **by position, not by name**. NONMEM skips the header and
  takes its names from `$INPUT`, so the two can disagree - in the bundled data
  file they do, `$INPUT`'s `DV` being the file's `LNDV` - and matching by name
  would silently read the wrong column. Columns beyond `$INPUT` are ignored,
  `DROP` columns still occupy a position, and either half of a `SYNONYM=REAL`
  pair may be used. `useInputNames = TRUE` returns the data renamed as NONMEM
  sees it. The single-character form (`IGNORE=C`) is a rule about the raw record
  rather than the data and cannot be applied; it is skipped with a warning,
  silently for the conventional `@` and `#`.

* **`oneHotEncode()`** adds `<covariate><sep><level>` dummy columns following
  the convention shared with NONMEM FREM data sets and
  `PMXFrem::createFREMModel()`. It is idempotent on data that already carries
  them. `getForestDFSCM()` and `getForestDFemp()` gain `oneHot` / `oneHotSep`,
  so raw multi-level columns are encoded before the parameter function runs and
  one function written against the dummies serves both workflows. `oneHot`
  defaults to `NULL`; output is unchanged without it.

## Keeping the reference row and the parameter function in agreement

`contRef` and `catRef` in `setupDfCovs()` and `setupDfRefRow()` now accept
either a single setting for every covariate, or a named list with an entry per
covariate and an optional `default`. `contRef` takes a number, `"mean"`,
`"median"` or `"model"`; `catRef` a level, `"mode"`, `"lowest"` or `"model"` -
for example `contRef = list(WT = 75, AGE = "mean", default = "median")`.
`"model"` reads the reference from the control stream through the new `model`
argument, which takes a `.mod` path or the list `createParamFunction()` returns.

This addresses a mismatch that was easy to miss. A parameter function must
return something when a covariate is inactive, and it falls back to the model's
own reference - the normalisation constant in `(WT/75)`, or the level in the
branch marked `; Most common`. `setupDfRefRow()` meanwhile took the median and
the modal level from the data. Where those differed, every row in which that
covariate was inactive was displaced from 1: on `run7`, by
`(75/85.4)^0.75 = 0.907` for CL, `(75/85.4) = 0.878` for V and
`1/(1-0.145) = 1.17` for Frel. Taking both from `$PK` removes the
disagreement. Defaults are unchanged, so existing plots are reproduced exactly.

Note that `"model"` as a bare scalar applies to every covariate, so one
covariate absent from `$PK` fails the call; use the named-list form. A
one-hot coded categorical needs a literal level rather than `"model"`, because
`$PK` names the dummies (`GENO1`, `GENO3`, `GENO4`) and never the covariate they
encode.

## Plots

* **`addStamp()`** appends a caption recording where and when a figure was
  produced - the working directory name, and inside a `knitr` chunk the input
  file and chunk label, plus the creation time - so a plot pasted into a report
  can be traced back to the code that made it. Segments that are unknown are
  omitted rather than left as empty path elements. It is a native
  reimplementation of the stamp `PhRame::add_stamp()` applies, written against
  `ggplot2` alone so the public packages can use it; `PhRame` is internal, so
  depending on it would break for external users. Only the "return the annotated
  plot" behaviour is reproduced - saving and printing are left to the caller.

## Bug fixes

* **A covariate could be read from the wrong column.** The missing-value
  preamble of a generated parameter function used `df$WT`, and `$`
  partial-matches on a data frame - so a data set carrying `WTKG` but no `WT`
  silently used `WTKG` as the covariate instead of falling back to the
  reference. Name families like this are ordinary; `run7`'s own `$INPUT` has
  `NCI`/`NCIL` and `RACEL`/`RACEL1`. The generated code now uses `df[["WT"]]`.

* **`MOD()` translated to the wrong arithmetic**, in two independent ways.
  `%%` binds tighter than `*` and `/` in R while `MOD()` is a call, so
  `MOD(A*B, C)` emitted `A * B %% C`; and Fortran `MOD()` truncates towards
  zero where R's `%%` floors, so they disagree for a negative first argument.
  `MOD(a, b)` now emits `a - b * trunc(a / b)`, which is correct for either
  sign and needs no special parenthesisation.

* **The documented rule 1 could never fire.** Covariates were taken to be
  `$INPUT` columns *never assigned* in `$PK`, but the primary reference rule
  reads `IF(WT.EQ.-99) WT = 75`, which assigns `WT` - so such a covariate was
  excluded before the rule was consulted, no preamble was emitted, and the
  generated function referred to an unbound variable. Covariates are now those
  `$PK` reads before assigning, which is what NONMEM does: data items are
  populated before `$PK` runs.

* **`covRef` was spliced in unchecked.** A misspelled name was silently ignored
  and left the derived reference in place, a non-numeric value produced source
  that parsed but failed when called, and a length-2 value died inside a
  formatting helper naming neither the argument nor the covariate. Names must
  now be covariates of the model and values a single finite number.

* **`getCovStats()` leaked `NA` into level counting and quantiles.**
  `x != missVal` is `NA` where `x` is `NA`, and indexing rows by a logical `NA`
  keeps an all-`NA` row rather than dropping it. A binary covariate with a
  genuine `NA` therefore silently changed shape - the `NA` counted as a third
  level, returning a nested one-hot list where the documentation promises a
  sorted vector - and a continuous one reached `quantile()`, whose `na.rm` is
  `FALSE`, and failed with an error naming neither the covariate nor the cause.
  `setupDfCovs()` inherited both.

* **A per-covariate length error was reported unhelpfully.**
  `contRef = c(70, 80)` was rejected clearly, but the same mistake written
  `contRef = list(WT = c(70, 80))` passed the check and failed much later with
  "replacement has 2 rows, data has 1". The length check now covers the named
  entry and `default` too.

* **Rounded quantiles could collapse onto one another.**
  `setupCovExpressionsList()` rounds its cut points with `signif()`, and when
  two distinct quantiles round to the same value the documented contrast
  inverts: the two rows then cover every subject instead of leaving those
  between them in neither. `minSubjects` cannot catch it, because both rows are
  non-empty. It now warns.

* **`verifyParamFunction()` with an explicitly named parameter:** naming a
  `secondary` parameter crashed with "subscript out of bounds" instead of
  warning that it could not be reconciled against the table. Naming a parameter
  the function does not return at all likewise crashed. Both now warn and report
  `PASS = NA`. Found through `covr::package_coverage()`, which exercises the
  path; the routine test suite did not, because the test helper had the same
  bug and masked it.

* **`nmResolveSecondary()` missed an unnamed constant.** An entry such as
  `list(source = "x", 100)` should have been rejected as unnamed, but the extra
  elements were selected with `v[setdiff(names(v), "source")]`, and indexing a
  list by `""` never matches in R - so the constant vanished before the name
  check ran and a confusing error surfaced further down. Rewritten as boolean
  indexing.

* **Windows / PSOCK parallelisation** in `getForestDFSCM()` /
  `getForestDFemp()`: with `ncores > 1` on a PSOCK platform the `foreach` loop
  could fail with "object not found", because its static global detection does
  not reliably follow the internal closure's free variables. The local
  environment is now exported explicitly, mirroring the fix in
  `PMXFrem::getExplainedVar()`. The cluster is also torn down through
  `on.exit()`, so it is released if the function errors mid-run. Fork-based
  parallelism and single-core runs are unaffected; results are identical.

* **`getForestDFSCM()` with a single covariate:** a `dfCovs` with one covariate
  column failed, because `dfCovs[i, ]` dropped to a vector. The remaining
  accesses now use `drop = FALSE`. (Original ticket dates to 1.0.5/1.0.6.)

* **`tibble` inputs** to `getForestDFSCM()` / `getForestDFemp()` failed with an
  unclear `vctrs` error, because the internal code relies on base-R `[`
  dropping a single-column selection to a vector. `dfCovs`, `dfRefRow` and
  `dfData` are now coerced with `as.data.frame()` on entry.

* **Reversed relative confidence intervals:** the `Q*_REL_REFFUNC` and
  `Q*_REL_REFFINAL` columns had their limits swapped when the parameter
  function returned a negative reference value. The relative quantiles are now
  computed from the ratio directly rather than by dividing the absolute
  quantiles by a possibly negative reference.

* **`1:length(x)` in loop bounds:** `1:0` is `c(1, 0)`, so a loop over an empty
  vector ran twice with out-of-range indices instead of not running. Replaced
  with `seq_along()` / `seq_len()` throughout.

* **Programmatic evaluation:** `getCovStats()`'s `idVar` is now evaluated with
  `rlang::sym()` rather than `rlang::ensym()`, so the function can be wrapped
  and called programmatically without scoping errors.

* **Unexported helpers:** `setupDfCovs()` and `setupDfRefRow()` were added
  without `NAMESPACE` entries or help pages, and so were not reachable as
  `PMXForest::setupDfCovs()`. Both are now exported and documented.

## Documentation

* The vignettes are a three-tier set: a Quick start, an end-to-end Walkthrough,
  and four deep dives (Forest plot inputs, R-coded models, secondary parameters,
  time-to-event models). They build with `rmarkdown::html_document` rather than
  `bookdown`.

* **The Quick start and the Walkthrough now teach the convenience functions.**
  Both generate the parameter function rather than writing one out, and the
  Walkthrough opens with `filterByModel()`, uses `setupDfCovs()` and
  `setupDfRefRow()` with references taken from the model, generates its
  empirical expressions with `setupCovExpressionsList()`, and closes with
  `addStamp()`. It prints `POINT` beside `REFFUNC` to show that with a
  model-derived reference the rows for covariates the model does not affect land
  at exactly 1.

* **The Forest plot inputs deep dive is the explicit layer.** Its introduction
  maps each convenience function to the primitives underneath, so it can be read
  in both directions. New sections cover what a parameter function must do -
  the contract, and why the `-99` guard is needed even when the data has no
  missing values - and what the generator refuses and why.

* Every exported function has a runnable `@examples` section against the
  bundled `SimVal` output, and eleven now carry a `@seealso` pointing at the
  vignette that puts them in context. The previously undocumented `forestPlot()`
  arguments `setSignEff`, `size` and `xlim` are documented.

* The README uses the current API and `system.file()` paths, so its example runs
  from an installed package rather than only from the source tree.

## Internal

* **Faster result assembly** in `getForestDFSCM()` / `getForestDFemp()`. Both
  built their per-parameter-vector result with a per-cell `data.frame()` and a
  growing `bind_rows()` loop. They now fill typed column vectors and construct
  the data frame once. Output is unchanged; the assembly step alone is roughly
  50-80x faster on a representative shape (10.7 s to ~0.15 s). Effect on total
  runtime depends on how expensive `functionList` is.

* Declared `rlang` in `Imports`, and `withr` in `Suggests`. Both were already
  used and neither was listed.

* Dropped the `table1` dependency. Its only use was `signif_pad()`, now a small
  base-R helper verified to produce identical output over a wide fuzz range.

* `setupForestPlotData()` is no longer exported; it is an implementation detail
  of `forestPlot()`, which is unaffected.

* Test coverage is 96.4% (measured for this release), with a `make coverage`
  target that fails below a 95% floor. The remaining gaps are the parallel
  (`ncores > 1`) branches and a few unreachable defensive guards. Added a `.lintr` configuration, and the package is
  now clean under `lintr` and `styler`.

# PMXForest 1.2.15

* Fixed the "small bootstrap" logic (.csv with n provided) to correctly prepend the base estimates from the .ext file, ensuring the output is always n+1 rows.
  * Refactored the data.frame input logic to act as a pure, generic MVRNORM sampler.
  * Added an is.numeric safety check to protect cov() when a data.frame is provided.
  * Updated documentation and unit tests to explicitly define the data.frame input as an exception to the n+1 rule (returning exactly n rows).

# PMXForest 1.2.14

* Added protective coding for the case when there is a missmatch between .ext and .cov.

# PMXForest 1.2.13
 
* Updated unit tests to improve coverage.
 
# PMXForest 1.2.12

* Update the DESCRIPTION fiel to add Pharmetheus as the copyright holder.

# PMXForest 1.2.11

* Fixed a bug in the handling of missing covariates in getCovStats().
* Updated the documentation for `noBaseThetas` in getForestDFSCM.
* Added unit tests for createInputForestData()

# PMXForest 1.2.10

* Making sure that IDs with missing data are only filtered out for the missing covariate concerned instead of across all covariates

# PMXForest 1.2.9

* Repeating the 1.2.8 release since something went wrong in the packaging.

# PMXForest 1.2.8

* Clarified documentations in some functions.
* Added an `xlim` argument to `forestPlot`.
* Changed the rounding function in `getCovStats`from `round` to `signif`.
* Clarified the documentation for `getCovStats()` so it is clear it requires a data.frame and not a file.
* Moved getForestDFFREM to PMXFrem
* Added argument parameterLabelsPrefix that will prepend parameterLabelsPrefix to the facet labels for the parameter panels.

# PMXForest 1.2.7

* Added a `NEWS.md` file to track changes to the package.

# PMXForest 1.2.6

* Fixed a bug related to the order of facet panels and facet labels when the user
  provides `parameters` to `forestPlot` that has a different order than `functionListName`.
* `statisticsLabels`should now be a string instead of a vector (see blow).

There is a potentially serious bug in PMXForest in versions prior to and including 1.2.4. It occurs when the user requests certain parameters using the parameters argument to the forestPlot function and when the order of the requested parameters is different from the order the parameters have in the return statement from the paramFunction. In this case there will be a missmatch between the order of the panels in the forest plot and the facet labels.:

    paramFunction <- function(thetas, df,…} {
    …
    return(c(CL,Frel,AUC))
    }
    
    functionListName <- c("CL","Frel","AUC")
    
    forestPlot(dfres,parameters=c(“Frel”,”CL”))

In this hypothetical example, the facet order will be Frel and CL but the face labels will have the same order as in `functionListName`, i.e. CL, Frel.

This has been fixed in PMXForest 1.2.6. Note that this fix necessitated a change to the `statisticsLabels` argument to the `forestPlot` function. Previously `statisticsLabels` expected a vector of facet labels for the statistics facets, e.g. `paste(“Statistics:”,parameters)`. In PMXForest 1.2.6, `statisticsLabels`  should be a string that will be prepended to the `parameterLabels` (default is "Statistics:").  `parameterLabels` is by default is the same as `parameters`. Here are a few examples:

    forestPlot(dfres) 
Facet order: CL, Frel, AUC
Facet labels for the plot panels: CL, Frel, AUC
Facet labels for the statistics panels:  Statistics: CL, Statistics: Frel, Statistics: AUC 

    forestPlots(dfres,statisticsLabels=“Stat:”)
Facet order: CL, Frel, AUC
Facet labels for the plot panels: CL, Frel, AUC
Facet labels for the statistics panels:  Stat: CL, Stat: Frel, Stat: AUC 

    forestPlots(dfres,statisticsLabels=“Stat:”,parameters=c(“Frel”,”CL”),parameterLabels=c(“F”,”Clearance”))
(The order of parameterLabels must match the order of parameters.)
Facet order: Frel, CL
Facet labels for the plot panels: F, Clearance
Facet labels for the statistics panels:  Stat: F, Stat: Clearance

