# Stress-testing PMXForest against your own models

PMXForest 1.3.0 adds two pieces of code that read NONMEM control streams:

- `createParamFunction()` — translates a `$PK` block into R source, so you no
  longer retype the model's algebra by hand
- `filterByModel()` — applies the `$DATA` `IGNORE`/`ACCEPT` statements, so a
  Forest plot describes the population the model actually saw

Both have only been exercised against the models we have in house, and those
are nearly all the same shape: one 1-compartment oral model with different
covariate decorations. Your project directories contain idioms nobody here
would think to invent, which is the point of asking.

`verifyPMXForest.R` points the new code at your runs and reports what happened.

## What it does not do

- **It never re-runs a model.** Everything is read from artefacts already on
  disk: the control stream, the `.lst`, the `.ext`, and any `$TABLE` output.
- **It never writes into a run directory.** Output goes to one directory you
  name, plus the R session's `tempdir()`.

## Running it

```r
source("verifyPMXForest.R")
verifyPMXForest("~/projects/xyz/Models")
```

or from a shell:

```sh
Rscript verifyPMXForest.R ~/projects/xyz/Models ~/projects/abc/Models
```

Point it at directories (scanned recursively for `.mod`, `.ctl`, `.con`) or at
individual files. More models is better; a hundred is not too many.

### Getting the right version installed

The script checks that PMXForest **1.2.15.9009** is loaded and stops if not.
To install it:

```r
PMXRenv::activate.unqualified.packages()
PMXRenv::install.unqualified.packages("PMXForest", repoName = "development")
```

`--install` does this for you, but note it **overwrites the PMXForest in your
unqualified library** — a change on disk, not just in the session. The script
prints the version and path it is about to replace so you can put it back, and
you will need a fresh R session afterwards.

## What it checks

| Tier | Check | Needs |
|---|---|---|
| A | `filterByModel()` reproduces the record, subject and observation counts NONMEM printed in the `.lst` | `.mod` + `.lst` + the `$DATA` file |
| B | `createParamFunction()` translates the model, the code parses and runs — or refuses honestly | `.mod` only |
| C | `verifyParamFunction()` reproduces the parameter values NONMEM tabled | a `$TABLE` containing the parameters |
| D | the generated function agrees with a `paramFunction` you wrote by hand | your function (optional) |

Tier A is the strongest check that runs everywhere: NONMEM records what it
actually did, so the comparison is exact and needs no `$TABLE`.

**A refusal is a pass.** The generator is supposed to refuse `$DES`,
compartment amounts, verbatim FORTRAN, `DO` and `CALL` rather than guess. What
we want to hear about is a *crash*, or a refusal of something it should have
handled.

### Tier D, if you have a hand-written paramFunction

This is the most valuable check, and it needs no `$TABLE`. If you already wrote
a `paramFunction` for a past Forest plot, put it next to the model as
`<model>.paramFunction.R`, or list them in a CSV:

```
modFile,funFile,funName
/path/run7.mod,/path/myParamFunction.R,paramFunction
```

```r
verifyPMXForest("~/projects", manifest = "map.csv")
```

## What to send back

Three files are written:

- `...-report.csv` — **this is the one to send back.** Redacted: no paths, no
  covariate names, no subject counts. Model structure (ADVAN, THETA count) and
  pass/fail per tier only.
- `...-local.rds` — everything verbatim, including paths, covariate names and
  raw error messages. **Stays on your machine.**
- `...-console.log` — what was printed.

The script prints the first rows of the CSV before you send it, so you can see
exactly what is in it. If your work is not sensitive and you would rather send
the detail, run with `full = TRUE`.

## If something looks wrong

Send the `-report.csv` and say which `modelId` looked odd. The `modelId` maps
back to a path in your local `.rds`, so we can ask you a precise question
without you having to send the model.

Please do **not** send control streams or data without checking what is in
them first.
