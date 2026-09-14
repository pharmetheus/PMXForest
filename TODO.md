# PMXForest — design backlog

Items parked for future consideration. Not a substitute for GitHub issues; move an
item there when it becomes active work.

## `$PRED` models

`createParamFunction()` refuses a `$PRED` model: `nmParsePK()` looks for a `$PK`
record and stops when there is none. The refusal is honest rather than
principled — most of the parser would work unchanged, since `$PRED` is the same
statement language.

What makes it more than a record-name change is that `$PRED` has no separation
between structure and residual error. A `$PK` block assigns parameters and
stops; a `$PRED` block computes the prediction and the `Y` in the same lines, so
"which assignments are parameters" stops being answerable from the block alone.
`parameters` would have to be mandatory, and the covariate discovery — which
rests on "read before assigned, therefore a data item" — needs re-examining
against a block that reads `DV` and `EPS()`.

Worth doing when someone asks for it with a real model in hand. Until then the
refusal names `$PRED` and says it must be handled by hand, which is the right
answer for a user who would otherwise get plausible-looking wrong source.
