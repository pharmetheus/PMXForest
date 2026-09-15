# PMXForest 1.3.0

Every bug and issue reported against PMXForest has been addressed in this
release. Beyond that, 1.3.0 is about removing the hand-written steps from the
Forest plot workflow: the parameter function can now be generated from the
control stream, the covariate table and the reference row can be built from the
data, and the whole API is documented with runnable examples and a three-tier
set of vignettes.

Two fixes change numbers in plots you have already made - SIR-based uncertainty
and the statistics table. Both are described under **Changes to existing
behaviour**; read that section before regenerating a figure you have shipped.

## Full documentation

The package had no vignettes and patchy help pages. It now has both.

* **A three-tier vignette set.** [Quick start](https://rpkgs-docs.pmx.one/PMXForest/articles/Part1-quick-start.html)
  gets a plot out of a NONMEM model in five steps;
  [Walkthrough](https://rpkgs-docs.pmx.one/PMXForest/articles/Part2-walkthrough.html) builds a publication figure
  end to end; and five deep dives go under the surface -
  [Preparing Forest plot inputs](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-forest-plot-inputs.html),
  [How the NONMEM model is parsed](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-nonmem-parsing.html),
  [Secondary parameters](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-secondary-parameters.html),
  [R-coded models](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-r-coded-models.html) and
  [Time-to-event models](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-tte-models.html). Quick start opens
  with a "Which vignette do I want?" table.

* **The top two vignettes teach the convenience functions**, and the Forest plot
  inputs deep dive is the explicit layer underneath: its introduction maps each
  convenience function to the primitives it composes, so it can be read in both
  directions. New sections there cover the contract a parameter function must
  satisfy - and why the `-99` guard is needed even when the data has no missing
  values - and what the generator refuses and why. The Walkthrough prints
  `POINT` beside `REFFUNC` to show that with a model-derived reference, the rows
  for covariates the model does not affect land at exactly 1.

* The vignettes build with `rmarkdown::html_document` rather than `bookdown`.

* **Runnable examples everywhere.** Every user-facing exported function has an
  `@examples` section that runs against the bundled `SimVal` output, so the help
  page is something you can execute rather than only read. Eleven carry a
  `@seealso` pointing at the vignette that puts them in context.

* **Cross-references and links work.** Help pages render their code references
  as working links to the functions they name, instead of showing the markup.
  The previously undocumented `forestPlot()` arguments `setSignEff`, `size` and
  `xlim` are documented.

## Convenience functions

Four of the five inputs `getForestDFSCM()` needs used to be assembled by hand.
Each was a place to make a silent mistake: a covariate table that did not match
the parameter function, a reference row taken from the data while the function
fell back to the model, algebra retyped out of `$PK`. These functions build
those inputs from the model and the data instead.

They are a starting point, not a straitjacket. Where the output does not match
what a particular plot needs, the intended workflow is to take what the function
returns and modify it programmatically - shown in
[Preparing Forest plot inputs](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-forest-plot-inputs.html).

### `setupDfCovs()`

`setupDfCovs()` composes `getCovStats()` and `createInputForestData()` into
one call: one row of
`dfCovs` per Forest plot row, continuous covariates at their 5th and 95th
percentiles, categorical ones one row per level.

* `conditionalCovs` names the covariates the others are conditioned on - fed
  state rather than fasted, patients rather than healthy volunteers. They get
  their own rows and sit at their reference on every other row, rather than at
  `missVal`.
* Statistics are computed on deduplicated data - one record per `idVar` - so
  longitudinal data cannot skew them.
* `useMissVal` toggles whether inactive cells hold `missVal` or a computed
  baseline.
* `contRef` and `catRef` accept either one setting for every covariate or a
  named list with an entry per covariate and an optional `default`. `contRef`
  takes a number, `"mean"`, `"median"` or `"model"`; `catRef` a level, `"mode"`,
  `"lowest"` or `"model"` - so `contRef = list(WT = 75, AGE = "mean",
  default = "median")` works.
* Records what it did in an attribute, which `setupDfRefRow()` reads back.

### `setupDfRefRow()`

`setupDfRefRow()` builds the reference row `getForestDF*()` expects, mapping
computed baseline
values onto the `dfCovs` geometry. `singleRef = TRUE` gives one row;
`singleRef = FALSE` gives a full overlay preserving background missingness.

* **`contRef = "model"` / `catRef = "model"` read the reference out of the
  control stream**, through the new `model` argument (a `.mod` path, or the list
  `createParamFunction()` returns). This closes a mismatch that was easy to
  miss: a parameter function falls back to the model's own reference when a
  covariate is inactive - the normalisation constant in `(WT/75)`, or the level
  in the branch marked `; Most common` - while `setupDfRefRow()` took the median
  and the modal level from the data. Where those differed, every row in which
  that covariate was inactive was displaced from 1: on the bundled `run7`, by
  `(75/85.4)^0.75 = 0.907` for CL, `(75/85.4) = 0.878` for V and
  `1/(1-0.145) = 1.17` for Frel.
* **It inherits `catRef` and `sep` from the `dfCovs` it is given**, so the two
  calls cannot silently disagree about the encoding. Passing them explicitly
  still overrides.
* `dfCovs` and `dfRefRow` need not cover the same covariates - a covariate
  present in `dfCovs` but absent from the reference row simply takes `missVal`.
* Defaults are unchanged, so existing plots reproduce exactly.
* Two caveats on `"model"`: as a bare scalar it applies to every covariate,
  so one covariate absent from `$PK` fails the call - use the named-list
  form; and a one-hot coded categorical needs a literal level, because `$PK`
  names the dummies (`GENO1`, `GENO3`, `GENO4`) and never the covariate they
  encode.

### `createParamFunction()`

Writing a `paramFunction` by hand means retyping the model's algebra - the
easiest place in the whole workflow to introduce an error no test would catch.
`createParamFunction()` reads the control stream and writes it for you, as
**source text you read, check and evaluate yourself**. Nothing is evaluated on
your behalf.

It reads the block the model defines its parameters in: `$PK` for a PREDPP
model, `$PRED` for one that codes the prediction by hand. What it accepts is
assignments, `IF`/`ELSE IF`/`ELSE`/`END IF`, the `.EQ.`/`.AND.` operator family
and the usual intrinsics.

* **Typical values.** Every `ETA(n)` is set to 0, and in `$PRED` every
  `EPS(n)`/`ERR(n)` with it, so the function returns typical values - what a
  Forest plot needs.

* **Covariate references, with their provenance.** All missing-covariate
  handling is hoisted into one annotated preamble, each reference taken from
  the control stream by four rules in decreasing order of confidence: explicit
  `IF(WT.EQ.-99)` handling; the branch PsN's `scm` marks `; Most common`; the
  branch assigning the identity value; the normalisation constant in `(WT/75)`
  or `(AGE-50)`. A weaker fifth rule proposes the level no `IF()` tests, and
  says so. A covariate matching no rule is an error you fix with `covRef`.

  A wrong reference displaces every row of a plot and nothing downstream can
  detect it, so each preamble line records which rule produced it and from
  which control-stream line.

* **`parameters` prunes the source, not just the return list.** Only the
  statements the requested parameters can reach are emitted, and the covariate
  set narrows with them - so no reference has to be derived, or justified, for
  a covariate that cannot affect the answer.

  One consequence worth knowing: the reference rules read the statements that
  survive pruning, so a covariate both `CL` and `V` use can get a different
  derived reference from `parameters = "V"` than from the whole block. Supply
  it through `covRef` if you need one value across separate calls.

* **`$PRED` models are supported, and require `parameters`.** A `$PK` block
  assigns parameters and stops; a `$PRED` block computes the prediction and the
  residual error in the same lines, so there is no structural way to tell which
  of its assignments you want. Name them. `Y` is an assignment like any other -
  ask for it and you get the typical prediction.

* **What separates, and what does not.** Both exponential-IIV idioms are
  recognised - `P = <expr> * EXP(ETA(n))` and the MU-referenced
  `P = EXP(MU_n + ETA(n))` that IMP and SAEM models use - because
  `EXP(a + eta)` is `EXP(a) * EXP(eta)`. An eta scaled inside the exponent, a
  Box-Cox transformed one, and inter-occasion variability in either of its
  codings do **not** separate, and get no `etaMap` entry rather than a claimed
  one. `tvMap` records each parameter's typical value as the control stream
  itself names it, so `verifyParamFunction()` never has to guess at a
  `TV`-prefixed column - `run7` has `MAT = TVMAT*EXP(ETA(5))`, where `TVMAT`
  really is `MAT`'s typical value, and then `D1 = MAT*(1-TVD1)`, where `TVD1`
  is a dimensionless fraction.

* **Arithmetic is translated faithfully rather than obviously.** `MOD(a, b)`
  emits `a - b * trunc(a / b)`: `%%` binds tighter than `*` and `/` in R while
  `MOD()` is a call, and Fortran truncates towards zero where R's `%%` floors,
  so the obvious translation is wrong twice over for a negative argument.

* **`OMEGA(i,j)` and `SIGMA(i,j)` resolve from the `.ext` file.** A model that
  builds a correlation by hand reads its own variance matrix -
  `SD_FREL = SQRT(OMEGA(2,2))`. These are constants once the model is
  estimated, so they fold to literals when `extFile` is supplied, and the call
  says what it needs when it is not. The `.ext` holds the lower triangle and
  the matrix is symmetric, so `OMEGA(2,3)` and `OMEGA(3,2)` both resolve.

  They fold to the *final estimates*, so the generated function holds them
  fixed. Across bootstrap or SIR samples, where the THETAs vary, these do not -
  each sample's own variance estimate is not used. For a parameter that reads
  `OMEGA` only to build a correlation this is usually what you want; where the
  variance genuinely drives the parameter, its uncertainty is not represented.

* **It refuses what it cannot translate faithfully** - `$DES`, compartment
  amounts `A(n)`, `DO`, `CALL`, and verbatim FORTRAN unless
  `ignoreVerbatim = TRUE` - naming the file and line rather than guessing. It
  also refuses a symbol read before anything assigns it: NONMEM neither
  initialises these variables nor clears them between data records, so such a
  model reads whatever the previous subject left behind, which a function
  evaluated one row at a time cannot reproduce.

  A **FREM** model is refused too, and that one is a trap rather than a
  limitation: its `$PK` translates perfectly, but its covariate effects live in
  `$OMEGA`, so the result would describe none of the covariates the model was
  built for. Use `PMXFrem::createFREMParamFunction()`.

* **`covRef` also pins a symbol NONMEM supplies rather than `$INPUT`.** A
  `$MIX` model reads `MIXNUM`, the subpopulation index, which the block never
  assigns and `$INPUT` never declares; `covRef = list(MIXNUM = 1)` says which
  subpopulation to generate for, and giving `MIXNUM` its own column in `dfCovs`
  plots each in turn. Pinning is always deliberate - the generator refuses
  rather than choose for you - and some symbols cannot be pinned meaningfully
  at all: `NEWIND` changes within an individual. Everything that needs pinning
  is reported in **one** error with a call you can paste straight in.

* **`secondary`** attaches quantities the block does not contain - AUC, `Cmax`,
  an event probability. Each entry is a line of R code
  (`secondary = list(AUC = "500 / CL")`) or the path to an `.R` file of
  arbitrary code, including a `deSolve` or `mrgsolve` simulation, whose text is
  **inlined** so the result stays self-contained. Entries are spliced inside
  `local({ ... })`, see `thetas`, `df` and every structural parameter by name,
  and are emitted in order so a later one may use an earlier one. An entry may
  be `list(source = <string>, dose = 100, tau = 12)`, emitting the named
  constants immediately before the code, so a reusable secondary file can be
  parametrised at the call site. Secondary names are appended to
  `functionListName` and recorded in `secondaryNames`, with `primaryNames`
  holding the structural parameters alone.

* **`verifyParamFunction()`** turns "do I trust this translation?" into a pass
  or a fail: it evaluates the generated function over a NONMEM `$TABLE` and
  compares it with the values NONMEM wrote. Table output is written after
  `IGNORE`/`ACCEPT`, so taking both parameters and covariates from the table
  removes any need to reproduce the `$DATA` filtering. It returns a single
  `TRUE`/`FALSE`, usable directly in an `if`, with the per-parameter table in
  `attr(., "checks")`.

  What it cannot check is a reference value: it takes covariate values from the
  table, so the reference - which applies only when a covariate is *absent* -
  is never on the path it tests. Read the preamble.

[How the NONMEM model is parsed](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-nonmem-parsing.html) is the full account: the grammar, both
`ETA` idioms and the ones that do not separate, the reference rules in order
with a worked example of one reading a structure backwards, every refusal and
why it is the honest answer, and the gaps stated rather than hidden. The
[Secondary parameters](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-secondary-parameters.html) vignette walks through the secondary-parameter
workflow end to end.

### `setupCovExpressionsList()`

`setupCovExpressionsList()` is the empirical counterpart of `setupDfCovs()`,
producing the
`covExpressionsList` and matching `cdfCovsNames` that `getForestDFemp()`
consumes.

* Continuous covariates split at the `probs` tails or the median
  (`contSplit = "median"`); multi-level categoricals emit one row per level.
* **Every generated expression is checked against the data** and the call stops
  if one selects fewer than `minSubjects` subjects (default 10), so an empty or
  near-empty Forest plot row cannot reach the figure unnoticed.
* It applies the same level and quantile rules as `getCovStats()` but implements
  them separately, because the two return different shapes - a change to one
  does not follow into the other.

### `getSamples()`

`getSamples()` is the uncertainty source. Not new, but substantially fixed and
now documented.

* **SIR files were being read through the wrong code path** - see **Changes to
  existing behaviour** below. This is the single change in 1.3.0 most likely to
  move a number in a plot you have already made.
* Gains `quiet` (default `FALSE`), which prints a message when a SIR file is
  detected, so the input type it chose is visible rather than assumed.

## Other notable additions

### Faster result assembly

`getForestDFSCM()` and `getForestDFemp()` built their per-parameter-vector
result with a per-cell `data.frame()` and a growing `bind_rows()` loop. They now
fill typed column vectors and construct the data frame once. Output is
unchanged; the assembly step alone is roughly 50-80x faster on a representative
shape - 10.7 s down to about 0.15 s. Effect on total runtime depends on how
expensive `functionList` is.

### `filterByModel()`

`filterByModel()` applies the `IGNORE` and `ACCEPT` statements in a control
stream's `$DATA`
record, returning the analysis data - the rows the model actually used.

* Summarising the raw file instead describes a population the model never saw:
  on the bundled `run7`, the file holds 964 subjects and the model reads 754,
  and CRCL's 5th percentile moves from 74.9 to 77.7.
* **Columns are matched by position, not by name.** NONMEM skips the header and
  takes its names from `$INPUT`, so the two can disagree - in the bundled data
  file they do, `$INPUT`'s `DV` being the file's `LNDV` - and matching by name
  would silently read the wrong column. Columns beyond `$INPUT` are ignored,
  `DROP` columns still occupy a position, and either half of a `SYNONYM=REAL`
  pair may be used.
* `.EQ.`/`.NE.` compare as text and `.EQN.`/`.NEN.` numerically, matching
  NONMEM. Where a text comparison selects nothing but the numeric reading
  would, the call warns rather than silently returning an empty result - that
  is almost always a value written in a different form from the data, such as
  `IGNORE=(OCC.EQ.1)` against a column holding `1.0`.

  For the numeric family, a `"."` placeholder is read as `0`, which is what
  NM-TRAN does with it: a NONMEM data file writes `.` wherever a field does not
  apply, and a data set full of them is ordinary rather than broken. The
  coercion is reported unless `quiet = TRUE`, because it is a reading of the
  data rather than something the file says literally, and it is applied only to
  the numeric comparisons - a column used by both families is still compared as
  text where the model compares it as text. Any other non-numeric value is
  refused: `.` is a convention, `"unknown"` is a data problem, and coercing it
  would drop rows the model kept.

  That reading is confirmed against NONMEM rather than assumed. On a model
  filtering 400,000 records down to 232,911 through twenty `IGNORE` conditions,
  one of them a numeric comparison on a column full of placeholders, the
  record, subject and observation counts reproduce the `.lst` exactly.
* `useInputNames = TRUE` returns the data renamed as NONMEM sees it.
* The single-character form (`IGNORE=C`) is a rule about the raw record rather
  than the data and cannot be applied; it is skipped with a warning, silently
  for the conventional `@` and `#`.

#### `verifyFilterByModel()`

`verifyFilterByModel()` checks a filtered data set against NONMEM's own
account of the run. Every `.lst` prints how many records, subjects and
observations were read once `$DATA` had been applied; those three numbers
come from NONMEM rather than from a reading of the control stream, so they
catch a filter that is wrong in the same way your own reading of `$DATA` is
wrong. It needs only artefacts already on disk and never re-runs anything:

```r
v <- verifyFilterByModel("run7.mod")
#> PASS - verifyFilterByModel: 3/3 check(s) (35131 record(s) removed)
#>         CHECK FILTERED NONMEM DIFF PASS
#>       RECORDS    33885  33885    0 TRUE
#>      SUBJECTS      754    754    0 TRUE
#>  OBSERVATIONS     6284   6284    0 TRUE
```

It returns a single `TRUE`/`FALSE`, usable directly in an `if`, with the
table in `attr(., "checks")`. The `.lst` and the `$DATA` file are found
beside the model unless given. Two things it will not let pass quietly: a
`$DATA` record that removed nothing agrees with the `.lst` whatever the
filter did, so the result is marked `informative = FALSE`; and where the
`.lst`'s echoed `$INPUT` or filter clauses no longer match the control
stream on disk, the counts describe a different model and it says so. That
comparison is by **position**, not by name: a `=DROP` item still occupies a
data column, so a `$INPUT` that gained or lost one since the run reads
identically by name while everything after it has shifted.

Finding the echoed control stream no longer depends on a `NM-TRAN MESSAGES`
line, which not every `.lst` carries. Several markers for the start of
NONMEM's own output are tried, and failing those the last `$`-record line
bounds it - a fixed window would truncate a long control stream mid-record
and compare half of one. A `.lst` whose echo cannot be read at all now says
so, rather than passing silently as if it had matched.


### `oneHotEncode()`

`oneHotEncode()` adds `<covariate><sep><level>` dummy columns following the
convention shared
with NONMEM FREM data sets and `PMXFrem::createFREMModel()`. It is idempotent on
data that already carries them.

* `getForestDFSCM()` and `getForestDFemp()` gain `oneHot` / `oneHotSep`, so raw
  multi-level columns are encoded before the parameter function runs and one
  function written against the dummies serves both the parametric and the
  empirical workflow.
* `oneHot` defaults to `NULL`; output is unchanged without it.

### `addStamp()`

`addStamp()` appends a caption recording where and when a figure was produced -
the working
directory name, and inside a `knitr` chunk the input file and chunk label, plus
the creation time - so a plot pasted into a report can be traced back to the
code that made it.

* Segments that are unknown are omitted rather than left as empty path elements.
* It is a native reimplementation of the stamp `PhRame::add_stamp()` applies,
  written against `ggplot2` alone so the public packages can use it; `PhRame` is
  internal, so depending on it would break for external users. Only the "return
  the annotated plot" behaviour is reproduced - saving and printing are left to
  the caller.

### The parser is reusable

`nmParsePK()` returns the parsed tree - `$PK` or `$PRED` - with `ETA()` and
`EPS()` intact, together with the covariates and their references, the THETA
count, the `etaMap` and the `tvMap`, without emitting source. `nmDeparse()`,
`nmFormatNum()` and `nmResolveSecondary()` are exported alongside it, so other
packages can build their own emitters on the same tested parser; `nmDeparse()`
takes `etaValue` and `epsValue` so an emitter decides what the random effects
become. PMXFrem uses it for FREM parameter functions.

## Changes to existing behaviour

* **A folded random effect now collapses the term it multiplies.** Setting
  `ETA(n)` to 0 leaves a literal `0` behind, and the constant folding now
  applies `0 * x -> 0`. In practice this is what makes inter-occasion
  variability disappear cleanly: `EXP(ETA(1) + ETA(2)*(1-OCC) + ETA(3)*OCC)`
  emits `CL <- TVCL` rather than `CL <- TVCL * exp(0 * (1 - OCC) + 0 * OCC)`,
  the occasion variable stops being a covariate that needs a reference value,
  and `tvMap` can record `TVCL` as the typical value so
  `verifyParamFunction()` can check the parameter after all.

  The one behavioural edge: `0 * NA` folds to `0` instead of propagating `NA`.
  For a typical-value function that is the right answer - the term contributes
  nothing whatever the occasion is.

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
  carry neither column, so the check was always false for them.

* **The statistics table shows more digits on the relative scale.** With neither
  `sigdigits` nor `decimals` set - the default - the table now uses 2 decimal
  places on the relative scale, and 2 significant digits on the absolute scale
  as before. Relative covariate effects cluster around 1, where 2 significant
  digits discarded the second decimal, so a row that read `1.0` now reads
  `1.03`. Only the printed numbers change; points, intervals and axis labels are
  unaffected. Calls passing `sigdigits` explicitly are unchanged. `forestPlot()`
  gains `decimals` alongside `sigdigits`; supply one or the other.

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

  **Check any `cdfCovsNames` you pass positionally.** Row labels are matched to
  rows by position, so for a covariate whose data order was not already sorted,
  the labels for its two levels now attach to the opposite levels - silently,
  with no error. On the bundled data exactly one covariate is affected
  (`ETHNIC`, which appears as `1, 0`); `NCIL`, `RACEL2` and `SEX` were already
  in sorted order and are unchanged. If you generate labels with
  `setupDfCovs()` or `setupCovExpressionsList()` rather than writing them out,
  nothing changes.

* **`getCovStats()` refuses a covariate with no non-missing value** instead of
  returning an empty result. `createInputForestData()` emits no rows for such a
  covariate, so it silently vanished from the plot. `setupDfCovs()` inherits the
  error. `refValues()` already behaved this way.

* **The statistics panel label gains its separating space.** `forestPlot()`'s
  `statisticsLabel` defaulted to `"Statistics:"` and is prepended verbatim, so
  every panel was titled `Statistics:CL (L/h)`. The default is now
  `"Statistics: "`. A label you pass yourself is still used exactly as given, so
  include the space you want.

* `setupForestPlotData()` is no longer exported; it is an implementation detail
  of `forestPlot()`, which is unaffected.

* The `table1` dependency is gone. Its only use was `signif_pad()`, now a small
  base-R helper verified to produce identical output over a wide fuzz range.
  `rlang` is declared in `Imports` and `withr` in `Suggests`; both were already
  used and neither was listed.

## Bug fixes

* **Reversed relative confidence intervals.** The `Q*_REL_REFFUNC` and
  `Q*_REL_REFFINAL` columns had their limits swapped when the parameter function
  returned a negative reference value. The relative quantiles are now computed
  from the ratio directly rather than by dividing the absolute quantiles by a
  possibly negative reference.

* **`getCovStats()` leaked `NA` into level counting and quantiles.**
  `x != missVal` is `NA` where `x` is `NA`, and indexing rows by a logical `NA`
  keeps an all-`NA` row rather than dropping it. A binary covariate with a
  genuine `NA` therefore silently changed shape - the `NA` counted as a third
  level, returning a nested one-hot list where the documentation promises a
  sorted vector - and a continuous one reached `quantile()`, whose `na.rm` is
  `FALSE`, and failed with an error naming neither the covariate nor the cause.
  `setupDfCovs()` inherited both.

* **`getForestDFSCM()` with a single covariate** failed, because `dfCovs[i, ]`
  dropped to a vector. The remaining accesses now use `drop = FALSE`. (Original
  ticket dates to 1.0.5/1.0.6.)

* **`tibble` inputs** to `getForestDFSCM()` / `getForestDFemp()` failed with an
  unclear `vctrs` error, because the internal code relies on base-R `[` dropping
  a single-column selection to a vector. `dfCovs`, `dfRefRow` and `dfData` are
  now coerced with `as.data.frame()` on entry.

* **Windows / PSOCK parallelisation** in `getForestDFSCM()` / `getForestDFemp()`:
  with `ncores > 1` on a PSOCK platform the `foreach` loop could fail with
  "object not found", because its static global detection does not reliably
  follow the internal closure's free variables. The local environment is now
  exported explicitly, mirroring the fix in `PMXFrem::getExplainedVar()`. The
  cluster is also torn down through `on.exit()`, so it is released if the
  function errors mid-run. Fork-based parallelism and single-core runs are
  unaffected; results are identical.

* **`1:length(x)` in loop bounds.** `1:0` is `c(1, 0)`, so a loop over an empty
  vector ran twice with out-of-range indices instead of not running. Replaced
  with `seq_along()` / `seq_len()` throughout.

* **Programmatic evaluation.** `getCovStats()`'s `idVar` is now evaluated with
  `rlang::sym()` rather than `rlang::ensym()`, so the function can be wrapped
  and called programmatically without scoping errors.

* **`setupDfCovs()` and `setupDfRefRow()` were not reachable.** Both were added
  without `NAMESPACE` entries or help pages, so `PMXForest::setupDfCovs()`
  failed. Both are now exported and documented.

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

## Quality

* **`filterByModel()` is checked against NONMEM's own account of the data.**
  The `.lst` prints how many records, subjects and observations NONMEM read
  after applying `$DATA`. The suite now reads those three numbers out of the
  bundled `.lst` files and asserts `filterByModel()` reproduces them, on `run7`
  (a `$DATA` with three `IGNORE` conditions) and on the time-to-event model
  (a different structure with no `MDV` column, so the observation count goes
  through the `EVID` fallback). It is the only check here that does not come
  from us, so it is the only one that can catch the filter agreeing with a
  hand-written subset while both are wrong. The numbers are read from the
  `.lst` rather than written into the test, and a companion test asserts that
  `run7`'s filter removes something, so the check cannot quietly become vacuous
  if the bundled data is ever replaced.

* Test coverage is 96.4% (measured for this release), with a `make coverage`
  target that fails below a 95% floor. The remaining gaps are the parallel
  (`ncores > 1`) branches and a few unreachable defensive guards.

* Added a `.lintr` configuration; the package is now clean under `lintr` and
  `styler`.

* The README uses the current API and `system.file()` paths, so its example runs
  from an installed package rather than only from the source tree.

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

