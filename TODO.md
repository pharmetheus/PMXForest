# PMXForest — design backlog

Items parked for future consideration. Not a substitute for GitHub issues; move an
item there when it becomes active work.

## `refModelValues()` uses a weaker covariate rule than the generator

`createParamFunction()` discovers covariates as "read before assigned, therefore
a data item", which is what NONMEM does — data items are populated before the
block runs. That is what makes `IF(WT.EQ.-99) WT = 75` a covariate carrying its
own reference value rather than a local.

`refModelValues()` (`R/refValues.R`), reached through
`setupDfRefRow(contRef = "model")`, instead uses "used anywhere **and** never
assigned anywhere":

```r
covs <- intersect(syms$used, setdiff(nmInputNames(mod), syms$assigned))
```

So a covariate the model guards — the exact idiom the primary reference rule
exists to read — is excluded before the rules are consulted, and
`setupDfRefRow()` silently falls back to the data median. Confirmed while adding
`$PRED` support: a model with `IF(WT.EQ.-99) WT = 75` and `CL = THETA(1)*(WT/75)`
gives `contRef = "model"` the median (69 on the test data) rather than 75.

The fix is to use `nmFirstUnboundUse()` as the main path does. It is not a
one-liner: `refModelValues()` parses the block itself and does not build the
`seen`/`assigned` sets that walk needs, and changing which covariates it returns
changes `setupDfRefRow()` output for existing users. Worth doing deliberately,
with the before/after on a real model.

## `forestPlot()` passes `size` to line layers — use `linewidth`

ggplot2 3.4.0 (November 2022) replaced `size` with `linewidth` for lines and
deprecated the old spelling. `forestPlot()` still passes `size` to two line
layers: the reference line (`geom_vline(size = ref_line_size)`,
`R/forestPlot.R:346`) and the confidence interval
(`geom_errorbarh(size = ci_line_size)`, `:349`). Every `forestPlot()` call warns
on ggplot2 >= 3.4.0; that is 130 of the 139 warnings in the local test suite.

- Change both to `linewidth =`. The user-facing arguments `ref_line_size` and
  `ci_line_size` keep their names, so nothing changes for callers.
- Add `ggplot2 (>= 3.4.0)` to `Imports`. `DESCRIPTION` sets no minimum today,
  and on an older ggplot2 `linewidth` is dropped with only an "Ignoring unknown
  parameters" warning, so the lines would silently lose their width. Check what
  the internal package manager serves before setting the floor.
- Unverified, check while in there: newer ggplot2 releases also deprecate
  `geom_errorbarh()` in favour of `geom_errorbar(orientation = "y")`.
- Plot tests assert on `ggplot_build()` data; check that the line width still
  reaches the built layer (`linewidth` column) rather than only that the
  warning went away.
