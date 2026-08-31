# PMXForest (development version)

## Breaking Changes
* `getCovStats()` now explicitly sorts binary covariate levels alphanumerically (e.g., always returning `c(0, 1)`) rather than returning them based on their order of appearance in the dataset. This ensures deterministic and reproducible row ordering in downstream forest plots.

## New Features
* **Streamlined Covariate Setup:** Introduced `setupDfCovs()`, a high-level wrapper that natively pipelines `getCovStats()` and `createInputForestData()`. This significantly reduces user friction during standard workflow setup.
* **Flexible Covariate Background States:** `setupDfCovs()` includes `additionalCovs` and `useMissVal` arguments. This allows users to include supplementary covariates (such as those from FREM workflows) and toggle whether inactive cells hold `missVal` or computed baseline references. Reference values are strictly computed on deduplicated data (one record per `idVar`) to prevent longitudinal sampling skew.
* **Automated Reference Row Generation:** Added `setupDfRefRow()` to generate the reference data frame required by `getForestDF` functions. Using a geometry-matching strategy, it maps computed baseline reference values (mode, median, or mean) onto the `dfCovs` structure. It supports both single-row outputs (`singleRef = TRUE`) and full matrix overlays (`singleRef = FALSE`) that perfectly preserve background missingness.
* **Shared one-hot encoder:** Added `oneHotEncode()`, which adds `<covariate><sep><level>` dummy columns (numeric `0`/`1`) to a data frame following the convention shared with NONMEM FREM data sets and `PMXFrem::createFREMModel()`. It is idempotent on data that already carries the dummy columns.
* **Explicit reference levels:** `getCovStats()`, `setupDfCovs()`, and `setupDfRefRow()` gained `refLevels` (a named list, e.g. `list(GENO = 2)`) and `sep`, so the generated one-hot columns can be aligned with a model whose reference category is not the lowest level. Defaults are unchanged (lowest level as reference, `"_"` separator).
* **`oneHot` in the forest-data functions:** `getForestDFSCM()` and `getForestDFemp()` gained optional `oneHot` and `oneHotSep` arguments. When `oneHot` is supplied, raw multi-level categorical columns are one-hot encoded before the parameter functions run - in `dfCovs`/`dfRefRow` for `getForestDFSCM()`, in `dfData` for `getForestDFemp()` - so a single parameter function written against the dummy columns serves both workflows. `oneHotSep` (default `"_"`) sets the name separator; use `""` for columns named like `GENO1`. `oneHot` defaults to `NULL` (no encoding, output unchanged).
* **Statistics-table number format:** `forestPlot()` gained a `decimals` argument alongside `sigdigits`; supply one or the other. With neither set (the default), the statistics table now uses **2 decimal places on the relative scale** and 2 significant digits on the absolute scale. Relative covariate effects cluster around 1, where 2 significant digits (the previous behaviour) discarded the second decimal - the default relative-scale table now shows values such as `1.03` rather than `1.0`. Only the printed numbers change; points, intervals and axis labels are unaffected. Calls that pass `sigdigits` explicitly are unchanged.

## Bug Fixes
* **`getForestDFSCM()` with a single covariate:** A `dfCovs` with only one covariate column failed - `dfCovs[i, ]` dropped to a vector, so the covariate column in the result was named `dfCovs[i, ]` (or, more recently, the call errored in `getCovNameString()` with "argument of length 0"). The three remaining `dfCovs[i, ]` accesses now use `drop = FALSE`; single- and multi-covariate `dfCovs` behave the same. (Original ticket dates back to PMXForest 1.0.5/1.0.6.)
* **`tibble` inputs to `getForestDFSCM()` / `getForestDFemp()`:** Passing `dfCovs` (or `dfRefRow`, or `dfData`) as a `tibble` failed with an unclear `vctrs` "Can't subset columns past the end" error, because the internal code relies on base-R `[` dropping a single-column selection to a vector, which a tibble does not do. These arguments are now coerced with `as.data.frame()` on entry, so `tibble` and `data.frame` inputs behave identically. (`PMXFrem::getForestDFFREM()` has the same pattern and needs the same fix.)
* **SIR raw_results handling in `getSamples()`:** The SIR-detection check looked for a `samples_order` column that PsN does not produce (the column is `sample_order`), so SIR files were routed through the bootstrap code path. That path filters on `ofv != 0` and therefore returned the full SIR *proposal* distribution instead of the importance-resampled parameter vectors (`resamples == 1`). `getSamples()` now returns the resampled vectors, with the final estimates prepended as the first row (as for the other inputs). **Forest plots built from SIR files will change**: for a well-converged SIR run the confidence intervals shift by a few percent (and not systematically in one direction); for a poorly initialised run the change can be larger. `getSamples()` also gains a `quiet` argument (default `FALSE`) that prints a message when a SIR file is detected.
* **Reversed relative confidence intervals:** Fixed a bug in `getForestDFSCM()` and `getForestDFemp()` where the `Q*_REL_REFFUNC` and `Q*_REL_REFFINAL` columns had their lower and upper limits swapped when the `functionList` function returned a negative reference value. The relative quantiles are now computed from the ratio directly instead of dividing the absolute quantiles by the (possibly negative) reference, so the interval endpoints stay correctly ordered.
* **Programmatic Evaluation:** Fixed a Non-Standard Evaluation (NSE) bug in `getCovStats()`. The `idVar` argument is now safely evaluated using `rlang::sym()` instead of `rlang::ensym()`, allowing the function to be properly wrapped and called programmatically without scoping errors.
* **Unexported helpers:** `setupDfCovs()` and `setupDfRefRow()` were added without a `NAMESPACE` entry or help page and were therefore not reachable as `PMXForest::setupDfCovs()` / `setupDfRefRow()`. Documentation has been regenerated so both functions are exported and documented.

## Documentation
* Restructured the vignettes into a three-tier set: a Quick-Start, an end-to-end Walkthrough, and three deep dives (covariate data preparation, R-coded models, time-to-event models). The vignettes now build with `rmarkdown::html_document` instead of `bookdown`.
* Every exported function now has a runnable `@examples` section that works on the bundled `SimVal` model output (`inst/extdata/SimVal`) rather than synthetic data or `\dontrun` snippets. Documented the previously undocumented `forestPlot()` arguments `setSignEff`, `size`, and `xlim`.

## Internal
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

