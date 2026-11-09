# PMXForest — design backlog

Items parked for future consideration. Not a substitute for GitHub issues; move an
item there when it becomes active work.

## Release plan: 1.4.1

Agreed 2026-11-09. Development on `epic/1.4.1`, released through #71. In
scope:

1. `forestPlot()` passes `statisticsLabels` by its full name; `Title` without
   the full stop that `R CMD check` notes. Done in this PR.
2. The walkthrough content from the 2026-10-26 cold read (below).
3. Reference conditioning of SCM Forest plots, and `filterByModel()` as a
   check on exclusions made elsewhere (friction item 3).
4. Additive etas in `nmEtaMap()`, so `createParamFunction()` maps
   `P = TVP + ETA(n)` and `verifyParamFunction()` checks it (PMXFrem T28).

Cold read of the walkthrough (2026-10-26), to be done in item 2: guidance
on choosing the uncertainty source (covariance, bootstrap, small bootstrap,
SIR) and parametric vs empirical plots; show `print_template`'s output and
what a mismatched reference row looks like rather than describing them; the
`catRef`/`sep` warning is given twice (Sections 1 and 2) without showing the
failure; and the subject and percentile numbers in Section 0 are asserted, not
printed. Done in item 2, which also corrected what `noVar` does.

## Walkthrough: second cold read (2026-11-09), left for later

Not addressed in 1.4.1:

- Where `"model"` references come from (the `WT/75` constant, the
  `; Most common` comments) is left implicit, and the displacement a
  data-median reference causes is asserted, not shown.
- `covRef` for `GENO1`/`GENO3`: "no `IF()` tests directly" does not match the
  `IF(GENO1.EQ.1)` lines in `$PK`; the warning it refers to is not shown.
- `POINT`: which parameter vector gives it is not said; `noBaseThetas` reads as
  "base thetas" but is all of them.
- Section 6: one record per subject drops time-varying covariates such as
  FOOD; the empirical cut points differ from the parametric ones without
  comment; `probs` and `cstrPackages` are not introduced.
- `thetas` in the `verifyParamFunction()` example is never defined.
- Style: em-dashes and possessives on things throughout the older text;
  "Section n" cross-references; NCI, FDC and the 2D6 levels undefined.

## Old branches on the public mirror

`Test3.5.3` and `miscMinorFixes` held nothing that `main` lacks and were
deleted on 2026-11-09. `noVarIssue` (1 commit, 2025-12-19) and `master`
(203 commits, 2023-06-25) are kept until someone reviews their commits.

## Friction log for 1.3.x - user feedback, October 2026

From the internal friction log. Items 1, 2 and 6 were answered in 1.4.0;
item 5 was praise.

3. **Why `filterByModel()` when pmxdata already excludes records?** The report
   already imports data and excludes records through pmxdata and record types;
   a second route risks two answers. Positioning `filterByModel()` as a check
   on the pmxdata result (where they differ, the modeller fixes it) may answer
   it. *1.4.1, item 3.*
4. **Why `createParamFunction()` when the report already has an mrgsolve
   translation?** A second translation doubles the checking, and the
   `secondary-parameters.R` hook is hard to understand and to debug. Asks for a
   vignette: building the parameter function from an existing mrgsolve model,
   fed `thetas` and `dfCovs`; the reporter offered code. Waits for that code.
7. **Figure captions should be quicker to edit**, e.g. built with `paste0()`
   from the covariate and parameter names. Concerns the report template, not
   the package; pass it on.
