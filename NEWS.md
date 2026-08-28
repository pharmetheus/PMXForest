# PMXForest (development version)

## Breaking Changes
* `getCovStats()` now explicitly sorts binary covariate levels alphanumerically (e.g., always returning `c(0, 1)`) rather than returning them based on their order of appearance in the dataset. This ensures deterministic and reproducible row ordering in downstream forest plots.

## New Features
* **Streamlined Covariate Setup:** Introduced `setupDfCovs()`, a high-level wrapper that natively pipelines `getCovStats()` and `createInputForestData()`. This significantly reduces user friction during standard workflow setup.
* **Flexible Covariate Background States:** `setupDfCovs()` includes `additionalCovs` and `useMissVal` arguments. This allows users to include supplementary covariates (such as those from FREM workflows) and toggle whether inactive cells hold `missVal` or computed baseline references. Reference values are strictly computed on deduplicated data (one record per `idVar`) to prevent longitudinal sampling skew.
* **Automated Reference Row Generation:** Added `setupDfRefRow()` to generate the reference data frame required by `getForestDF` functions. Using a geometry-matching strategy, it maps computed baseline reference values (mode, median, or mean) onto the `dfCovs` structure. It supports both single-row outputs (`singleRef = TRUE`) and full matrix overlays (`singleRef = FALSE`) that perfectly preserve background missingness.

## Bug Fixes
* **Reversed relative confidence intervals:** Fixed a bug in `getForestDFSCM()` and `getForestDFemp()` where the `Q*_REL_REFFUNC` and `Q*_REL_REFFINAL` columns had their lower and upper limits swapped when the `functionList` function returned a negative reference value. The relative quantiles are now computed from the ratio directly instead of dividing the absolute quantiles by the (possibly negative) reference, so the interval endpoints stay correctly ordered.
* **Programmatic Evaluation:** Fixed a Non-Standard Evaluation (NSE) bug in `getCovStats()`. The `idVar` argument is now safely evaluated using `rlang::sym()` instead of `rlang::ensym()`, allowing the function to be properly wrapped and called programmatically without scoping errors.
* **Unexported helpers:** `setupDfCovs()` and `setupDfRefRow()` were added without a `NAMESPACE` entry or help page and were therefore not reachable as `PMXForest::setupDfCovs()` / `setupDfRefRow()`. Documentation has been regenerated so both functions are exported and documented.

## Documentation
* Restructured the vignettes into a three-tier set: a Quick-Start, an end-to-end Walkthrough, and three deep dives (covariate data preparation, R-coded models, time-to-event models). The vignettes now build with `rmarkdown::html_document` instead of `bookdown`.

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

