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
