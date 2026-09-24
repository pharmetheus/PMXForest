# PMXForest (development version)

## Other notable additions

- `filterByModel()` gains `dropNoObs`, which also removes subjects left with
  no observation record - subjects with only dose or other event records.
  NONMEM reads them but they do not inform the estimates, so they need not
  enter the covariate summaries. Observations are identified the way NONMEM
  identifies them: from `MDV`, `EVID`, or the dose items `AMT`, `RATE` and
  `SS`, whichever the model declares.

## Changes to existing behaviour

- `verifyFilterByModel()` counts observations under the `$INPUT` names rather
  than the data file's, uses the same rule as `dropNoObs`, and reads an empty
  `MDV` field as 0, as NM-TRAN does.

# PMXForest 1.3.0

Every bug and issue reported against PMXForest has been addressed in this
release. Beyond that, 1.3.0 is about removing the hand-written steps from the
Forest plot workflow:

- the parameter function can now be generated from the control stream
- the covariate table and the reference row can be built from the data

The other big user facing change is that the whole package is documented with
runnable examples and a three-tier set of vignettes.

Two fixes change numbers in plots you have already made - SIR-based uncertainty
and the statistics table. Both are described under **Changes to existing
behaviour**; read that section before regenerating a figure you have shipped.

## Full documentation

The package now has a completely reworked set of vignettes:

* **A three-tier vignette set.**

  - [Quick start](https://rpkgs-docs.pmx.one/PMXForest/articles/Part1-quick-start.html)
    gets a plot out of a NONMEM model in five steps
  - [Walkthrough](https://rpkgs-docs.pmx.one/PMXForest/articles/Part2-walkthrough.html)
    builds a publication figure end to end
  - Deep Dives go into more details
    - [Preparing Forest plot inputs](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-forest-plot-inputs.html)
    - [How the NONMEM model is parsed](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-nonmem-parsing.html)
    - [Secondary parameters](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-secondary-parameters.html)
    - [R-coded models](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-r-coded-models.html)
    - [Time-to-event models](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-tte-models.html)

* **Runnable examples everywhere.** Every user-facing exported function has an
  `@examples` section that runs against the bundled `SimVal` output, so the help
  page is something you can execute rather than only read.

* **Cross-references and links work.** Help pages render their code references
  as working links to the functions they name, instead of showing the markup.

## Convenience functions

Four of the five inputs `getForestDFSCM()` needs used to be assembled by hand.
Each was a place to make a silent mistake: a covariate table that did not match
the parameter function, a reference row taken from the data while the function
fell back to the model, algebra retyped out of `$PK`. There are now convenience
functions that build those inputs from the model and the data instead.

The convenience functions are starting points, not straitjackets. If the output
from these functions does not match what a particular plot needs, the intended
workflow is to take what the function returns and modify it programmatically -
shown in
[Preparing Forest plot inputs](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-forest-plot-inputs.html).

### `setupDfCovs()`

`setupDfCovs()` is a wrapper for `getCovStats()` and `createInputForestData()`
and outputs one row of `dfCovs` per Forest plot row. The output is also used as
input to `setupDfRefRow()`.

### `setupDfRefRow()`

`setupDfRefRow()` builds the reference row `getForestDF*()` expects, mapping
computed baseline values onto the `dfCovs` geometry. `singleRef = TRUE` gives
one row; `singleRef = FALSE` gives a full overlay preserving background
missingness.

### `createParamFunction()`

Writing a `paramFunction` by hand means retyping the model's algebra - the
easiest place in the whole workflow to introduce an error no test would catch.
`createParamFunction()` reads the control stream and writes a `paramFunction`
for you, as **source text you read, check and evaluate yourself**.

`createParamFunction()` takes an argument `secondary` that indicates how
secondary parameters should be derived from the primary. Each entry is a line of
R code (`secondary = list(AUC = "500 / CL")`) or the path to an `.R` file of
arbitrary code, including a `deSolve` or `mrgsolve` simulation, whose text is
**inlined** so the result stays self-contained.

**`verifyParamFunction()`** can be used to verify the correctness of the
generated `paramFunction`, by evaluating it over a NONMEM `$TABLE` and comparing
the result with the values NONMEM wrote.

[How the NONMEM model is parsed](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-nonmem-parsing.html) is the
full account of how `createParamFunction()` parses NONMEM models. The
[Secondary parameters](https://rpkgs-docs.pmx.one/PMXForest/articles/Part3-deep-dive-secondary-parameters.html) vignette
walks through the secondary-parameter workflow end to end.

### `setupCovExpressionsList()`

`setupCovExpressionsList()` is the empirical counterpart of `setupDfCovs()`,
producing the `covExpressionsList` and matching `cdfCovsNames` that
`getForestDFemp()` consumes.

### `getSamples()`

`getSamples()` is the uncertainty source. Not new, but substantially fixed and
now documented.

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
stream's `$DATA` record, returning the analysis data - the rows the model
actually used.

#### `verifyFilterByModel()`

`verifyFilterByModel()` checks a filtered data set against NONMEM's own
account of the run. Every `.lst` prints how many records, subjects and
observations were read once `$DATA` had been applied; those three numbers
come from NONMEM rather than from a reading of the control stream, so they
catch a filter that is wrong.

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
  change can be larger.

* **The statistics table shows more digits on the relative scale.** With neither
  `sigdigits` nor `decimals` set - the default - the table now uses 2 decimal
  places on the relative scale, and 2 significant digits on the absolute scale
  as before.

* **`refLevels` is deprecated in favour of `catRef`.** `getCovStats()`,
  `setupDfCovs()`, `setupDfRefRow()` and `setupCovExpressionsList()` all take
  `catRef`, symmetric with `contRef`. `refLevels` still works and forwards to
  `catRef` with a deprecation warning; supplying both is an error.

* **`getCovStats()` sorts binary covariate levels.** They were returned in order
  of appearance in the data, so row order in a Forest plot depended on which
  subject came first. Now always `c(0, 1)`.

* **`getCovStats()` refuses a covariate with no non-missing value** instead of
  returning an empty result.

* **The statistics panel label gains its separating space.** `forestPlot()`'s
  `statisticsLabel` defaulted to `"Statistics:"` and is prepended verbatim, so
  every panel was titled `Statistics:CL (L/h)`. The default is now
  `"Statistics: "`.

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

