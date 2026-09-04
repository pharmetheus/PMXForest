# PMXForest 1.2.15.9007 (development version)

## Breaking Changes
* **`refLevels` is deprecated, replaced by `catRef`.** `getCovStats()`, `setupDfCovs()`, `setupDfRefRow()` and `setupCovExpressionsList()` all gain `catRef`, symmetric with `contRef`. `refLevels` still works and forwards to `catRef`, with a deprecation warning; supplying both is an error.
* **An explicit `catRef` (or `refLevels`) now sets the reference *value* as well as the one-hot column names.** Previously `refLevels` only chose which level was dropped when encoding, while the reference state was always the most common level in the data - so `refLevels = list(RACE = 1)` produced columns keyed on level 1 but a background sitting at the modal level. The two are now one setting: the level you name is both the level dropped when encoding and the state the reference takes. Calls that do **not** pass `refLevels`/`catRef` are unaffected - the encoding still drops the lowest level and the reference still comes from the mode.
* `getCovStats()` now explicitly sorts binary covariate levels alphanumerically (e.g., always returning `c(0, 1)`) rather than returning them based on their order of appearance in the dataset. This ensures deterministic and reproducible row ordering in downstream forest plots.

## New Features
* **Streamlined Covariate Setup:** Introduced `setupDfCovs()`, a high-level wrapper that natively pipelines `getCovStats()` and `createInputForestData()`. This significantly reduces user friction during standard workflow setup.
* **Flexible Covariate Background States:** `setupDfCovs()` includes `additionalCovs` and `useMissVal` arguments. This allows users to include supplementary covariates (such as those from FREM workflows) and toggle whether inactive cells hold `missVal` or computed baseline references. Reference values are strictly computed on deduplicated data (one record per `idVar`) to prevent longitudinal sampling skew.
* **Automated Reference Row Generation:** Added `setupDfRefRow()` to generate the reference data frame required by `getForestDF` functions. Using a geometry-matching strategy, it maps computed baseline reference values (mode, median, or mean) onto the `dfCovs` structure. It supports both single-row outputs (`singleRef = TRUE`) and full matrix overlays (`singleRef = FALSE`) that perfectly preserve background missingness.
* **Shared one-hot encoder:** Added `oneHotEncode()`, which adds `<covariate><sep><level>` dummy columns (numeric `0`/`1`) to a data frame following the convention shared with NONMEM FREM data sets and `PMXFrem::createFREMModel()`. It is idempotent on data that already carries the dummy columns.
* **Explicit reference levels:** `getCovStats()`, `setupDfCovs()`, and `setupDfRefRow()` gained `refLevels` (a named list, e.g. `list(GENO = 2)`) and `sep`, so the generated one-hot columns can be aligned with a model whose reference category is not the lowest level. Defaults are unchanged (lowest level as reference, `"_"` separator).
* **`oneHot` in the forest-data functions:** `getForestDFSCM()` and `getForestDFemp()` gained optional `oneHot` and `oneHotSep` arguments. When `oneHot` is supplied, raw multi-level categorical columns are one-hot encoded before the parameter functions run - in `dfCovs`/`dfRefRow` for `getForestDFSCM()`, in `dfData` for `getForestDFemp()` - so a single parameter function written against the dummy columns serves both workflows. `oneHotSep` (default `"_"`) sets the name separator; use `""` for columns named like `GENO1`. `oneHot` defaults to `NULL` (no encoding, output unchanged).
* **Statistics-table number format:** `forestPlot()` gained a `decimals` argument alongside `sigdigits`; supply one or the other. With neither set (the default), the statistics table now uses **2 decimal places on the relative scale** and 2 significant digits on the absolute scale. Relative covariate effects cluster around 1, where 2 significant digits (the previous behaviour) discarded the second decimal - the default relative-scale table now shows values such as `1.03` rather than `1.0`. Only the printed numbers change; points, intervals and axis labels are unaffected. Calls that pass `sigdigits` explicitly are unchanged.
* **Filtering a data set the way the model does:** Added `filterByModel()`, which reads the `IGNORE` and `ACCEPT` statements from a control stream's `$DATA` record and applies them, returning the rows the model actually used. Summarising the raw analysis file instead gives quantiles and reference values for a population the model never saw: on the bundled `run7` model the file holds 964 subjects but the model reads 754, and CRCL's 5th percentile moves from 74.9 to 77.7.

  Columns are matched **by position, not by name**. NONMEM skips the header line and takes its names from `$INPUT`, so the two can disagree - in the bundled data file they do, with `$INPUT`'s `DV` being the file's `LNDV` - and filtering by name would silently read the wrong column. Columns beyond `$INPUT` are ignored, `DROP` columns still occupy a position, and a `SYNONYM=REAL` pair can be referred to by either name. `useInputNames = TRUE` returns the data as NONMEM sees it, renamed. The single-character form (`IGNORE=C`) is a rule about the raw record rather than the data and cannot be applied; it is skipped with a warning, silently for the conventional `@` and `#` header markers.
* **Reference values can be taken from the model, or set per covariate:** `contRef` and `catRef` in `setupDfCovs()` and `setupDfRefRow()` now accept either a single setting applied to every covariate, or a named list with an entry per covariate and an optional `default` component. `contRef` takes a number, `"mean"`, `"median"` or `"model"`; `catRef` takes a level, `"mode"`, `"lowest"` or `"model"`. For example `contRef = list(WT = 75, AGE = "mean", default = "median")`. `"model"` reads the reference out of the NONMEM control stream through the new `model` argument, which accepts either a `.mod` path or the list returned by `createParamFunction()`.

  This addresses a mismatch that was easy to miss. A parameter function must return something when a covariate is inactive, and it falls back to the model's own reference - the normalisation constant in `(WT/75)`, or the level in the branch PsN's scm marks `; Most common`. `setupDfRefRow()` meanwhile took the median weight and the modal level from the data. Where those differed, every row in which that covariate was inactive was displaced from 1: on the bundled `run7` model, by `(75/85.4)^0.75 = 0.907` for CL, `(75/85.4)^1 = 0.878` for V, and `1/(1-0.145) = 1.17` for Frel. Setting `contRef = "model"` and `catRef = "model"` takes both from `$PK` so the reference row and the parameter function cannot disagree. The defaults are unchanged, so existing plots are reproduced exactly.
* **Parameter functions generated from the control stream:** Added `createParamFunction()`, which translates the `$PK` block of a NONMEM control stream into R source for a `paramFunction` of the form `function(thetas, df, ...)`. It returns the source as text for you to read, check and edit - nothing is evaluated. Every `ETA(n)` is set to 0, so the function returns typical values. All missing-covariate handling is hoisted into a single annotated preamble at the top of the generated function, with each reference value taken from the control stream itself: explicit `IF(WT.EQ.-99)` handling, the branch PsN's scm marks `; Most common`, the branch assigning the identity value, or the normalisation constant in `(WT/75)` / `(AGE-50)`. A weaker fifth rule proposes the level no `IF()` tests and warns; a covariate matching no rule is an error, and can be supplied through `covRef`. The parser accepts assignments, `IF` statements and closed-form arithmetic, and refuses anything else - `$DES`, compartment amounts `A(n)`, verbatim FORTRAN, `DO` loops, `CALL` - naming the file and line rather than guessing.

  Secondary parameters (AUC, Cmax, event probabilities - anything reached through `$ERROR` or `$DES`) are supplied through the `secondary` argument: a named list where each entry is either a line of R code (`secondary = list(AUC = "df$DOSE / CL")`) or the path to an `.R` file of arbitrary code, including a `deSolve` or `mrgsolve` simulation (`secondary = list(CMAX = "cmax.R")`). A file's text is **inlined** into the generated source, so the result stays a self-contained artifact. Each entry is spliced in inside `local({ ... })`, so it sees `thetas` / `df` / `...` and every structural parameter by name (covariate columns as `df$NAME`) while its own temporaries do not leak; entries are emitted in order so a later one may use an earlier one. The names are appended to `functionListName` (so `getForestDFSCM()` picks them up) and recorded in the new `secondaryNames` element, with `primaryNames` holding the `$PK` parameters alone. A string that looks like a path (`.R` / `.r`) but does not exist is an error rather than a mistyped snippet, and every entry is `parse()`d up front. An entry's value may also be a list `list(source = <string>, dose = 100, tau = 12, ...)` - the `source` is the code/file, and every other named atomic element is emitted as `name <- value` immediately before it inside the same `local({ })`, so a shipped or reusable secondary file can be parametrised at the call site instead of hard-coding constants; a bare string is just `source` with no constants. The resolver is exported as `nmResolveSecondary()` for reuse. Without `secondary`, the emitted source keeps the marked extension point for adding them by hand.
* **The `$PK` parse is now reusable:** the front end of `createParamFunction()` is exported as `nmParsePK()`, which reads a control stream and returns the parsed `$PK` statement tree (with `ETA()` references intact), the covariates and their reference values, the THETA count and the exponential-IIV `etaMap` - without emitting any source. `nmDeparse()` (render one expression node as R source, choosing the substitution for `ETA()`) and `nmFormatNum()` are exported alongside it. This lets other packages build their own emitters on the same tested parser; PMXFrem uses it for FREM parameter functions. `createParamFunction()` now calls `nmParsePK()` internally and is otherwise unchanged.
* **Checking a generated parameter function:** Added `verifyParamFunction()`, which evaluates the generated function over a NONMEM `$TABLE` file and compares the result with the values NONMEM itself wrote. Table output is written after `IGNORE`/`ACCEPT` have been applied, so taking both the parameters and the covariates from the table removes any need to reproduce the `$DATA` filtering. Where a parameter is written `P = <expr> * EXP(ETA(n))`, the tabled individual value is divided by `exp(ETA(n))` to recover the typical value; a `TVP` column is used directly. Anything missing from the table produces a loud warning and a `PASS` of `NA` rather than a silent comparison. It returns a single `TRUE` / `FALSE` (`TRUE` only if every requested parameter was checked and passed, so it can be used directly in an `if`) with the per-parameter table attached as `attr(., "checks")`. By default it checks only the `$PK` parameters; `secondary` parameters are skipped, since a secondary quantity such as AUC is generally not a `$TABLE` column and not reducible to a typical value.
* **Empirical covariate setup:** Added `setupCovExpressionsList()`, the empirical-workflow counterpart of `setupDfCovs()`. It turns a data frame plus a vector of covariate names into the named `covExpressionsList` and a matching `cdfCovsNames` label vector consumed by `getForestDFemp()`, reusing the deduplicated level and quantile logic of `getCovStats()`. Continuous covariates split at the `probs` quantile tails (or the median via `contSplit = "median"`); multi-level categoricals emit one row per level (`includeReference = FALSE` drops the reference). Each entry in `additionalCovs` gets its own rows and, in addition, a fixed condition on it (a level for categoricals; a `prob` or `value` split with a direction for continuous) is combined into every other covariate's expression. Every generated expression is checked against the deduplicated data and the function stops if any selects fewer than `minSubjects` (default 10) subjects.

## Bug Fixes
* **`getForestDFSCM()` with a single covariate:** A `dfCovs` with only one covariate column failed - `dfCovs[i, ]` dropped to a vector, so the covariate column in the result was named `dfCovs[i, ]` (or, more recently, the call errored in `getCovNameString()` with "argument of length 0"). The three remaining `dfCovs[i, ]` accesses now use `drop = FALSE`; single- and multi-covariate `dfCovs` behave the same. (Original ticket dates back to PMXForest 1.0.5/1.0.6.)
* **`tibble` inputs to `getForestDFSCM()` / `getForestDFemp()`:** Passing `dfCovs` (or `dfRefRow`, or `dfData`) as a `tibble` failed with an unclear `vctrs` "Can't subset columns past the end" error, because the internal code relies on base-R `[` dropping a single-column selection to a vector, which a tibble does not do. These arguments are now coerced with `as.data.frame()` on entry, so `tibble` and `data.frame` inputs behave identically. (`PMXFrem::getForestDFFREM()` has the same pattern and needs the same fix.)
* **SIR raw_results handling in `getSamples()`:** The SIR-detection check looked for a `samples_order` column that PsN does not produce (the column is `sample_order`), so SIR files were routed through the bootstrap code path. That path filters on `ofv != 0` and therefore returned the full SIR *proposal* distribution instead of the importance-resampled parameter vectors (`resamples == 1`). `getSamples()` now returns the resampled vectors, with the final estimates prepended as the first row (as for the other inputs). **Forest plots built from SIR files will change**: for a well-converged SIR run the confidence intervals shift by a few percent (and not systematically in one direction); for a poorly initialised run the change can be larger. `getSamples()` also gains a `quiet` argument (default `FALSE`) that prints a message when a SIR file is detected.
* **Reversed relative confidence intervals:** Fixed a bug in `getForestDFSCM()` and `getForestDFemp()` where the `Q*_REL_REFFUNC` and `Q*_REL_REFFINAL` columns had their lower and upper limits swapped when the `functionList` function returned a negative reference value. The relative quantiles are now computed from the ratio directly instead of dividing the absolute quantiles by the (possibly negative) reference, so the interval endpoints stay correctly ordered.
* **Programmatic Evaluation:** Fixed a Non-Standard Evaluation (NSE) bug in `getCovStats()`. The `idVar` argument is now safely evaluated using `rlang::sym()` instead of `rlang::ensym()`, allowing the function to be properly wrapped and called programmatically without scoping errors.
* **Unexported helpers:** `setupDfCovs()` and `setupDfRefRow()` were added without a `NAMESPACE` entry or help page and were therefore not reachable as `PMXForest::setupDfCovs()` / `setupDfRefRow()`. Documentation has been regenerated so both functions are exported and documented.

## Documentation
* Restructured the vignettes into a three-tier set: a Quick-Start, an end-to-end Walkthrough, and four deep dives (forest-plot inputs, R-coded models, secondary parameters, time-to-event models). The vignettes now build with `rmarkdown::html_document` instead of `bookdown`.
* New **"Secondary Parameters"** deep dive: a two-part walkthrough of the `secondary` argument. Part 1 derives `AUC` / an elimination rate / a half-life in closed form and runs at build time; Part 2 computes a steady-state `Cmax` with an `mrgsolve` simulation and is shown but not executed (its output is illustrative), with the mrgsolve version requirements spelled out.
* The "forest-plot inputs" deep dive replaces the earlier "covariate data preparation" vignette. It follows the three inputs `getForestDFSCM()` needs - `dfCovs`, the parameter function, and `dfRefRow` - end to end, and shows `createParamFunction()` and the `contRef`/`catRef` `"model"` option keeping the reference row and the parameter function in agreement.
* Every exported function now has a runnable `@examples` section that works on the bundled `SimVal` model output (`inst/extdata/SimVal`) rather than synthetic data or `\dontrun` snippets. Documented the previously undocumented `forestPlot()` arguments `setSignEff`, `size`, and `xlim`.

## Internal
* Declared `rlang` in `Imports`. It was already used through `rlang::sym()` in the deduplication step of `getCovStats()`, `setupDfCovs()`, `setupDfRefRow()` and `setupCovExpressionsList()` but was not listed, which `R CMD check` flagged.
* Dropped the `table1` dependency. Its only use was `signif_pad()` in the statistics-table formatting, now a small base-R helper (verified to produce identical output).
* `setupForestPlotData()` is no longer exported (`@keywords internal`). It is an implementation detail of `forestPlot()`, which is unaffected.
* Test coverage raised from 94.9% to 98.7% (every source file now at or above 98%). Added a `make coverage` target that fails below a 95% floor. The remaining gaps are the parallel (`ncores > 1`) branches and a few unreachable defensive guards.

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

