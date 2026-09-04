# Secondary parameter: steady-state Cmax after repeated oral dosing, via mrgsolve.
#
# Inlined by createParamFunction(secondary = list(CMAX = <this>)) into the body
# of the generated parameter function, wrapped in local({ }). It therefore sees,
# by name:
#   * thetas            - the NONMEM THETA vector
#   * df                - the one-row covariate data frame (columns as df$NAME)
#   * ...               - anything getForestDFSCM() forwards
#   * CL, V, MAT, FREL, KA, D1, ...  - the structural $PK parameters, already
#                         computed above in the generated function
# The value of the LAST expression below becomes CMAX.
#
# Regimen constants. Pass them with the config-list form to override, e.g.
#   secondary = list(CMAX = list(source = "<this>", dose = 100, tau = 12, n = 10))
# When a constant is supplied it is bound (dose <- 100, ...) immediately above
# this code, so exists(name, inherits = FALSE) sees it; otherwise the default
# below is used. inherits = FALSE keeps a like-named $PK variable from leaking in.
dose <- if (exists("dose", inherits = FALSE)) dose else 80    # amount per dose
tau  <- if (exists("tau",  inherits = FALSE)) tau  else 24    # dosing interval (h)
n    <- if (exists("n",    inherits = FALSE)) n    else 7      # doses to steady state

## ---- model: compiled once, then re-used from an on-disk cache --------------
## mcode_cache() keeps the compiled shared object between calls, so the C++
## build cost is paid once per session rather than on every bootstrap sample.
.code <- "
$PARAM CL = 1, V = 10, KA = 1
$CMT   ABS CENT
$ODE
dxdt_ABS  = -KA*ABS;
dxdt_CENT =  KA*ABS - (CL/V)*CENT;
$TABLE double CP = CENT / V;
$CAPTURE CP
"
.mod <- mrgsolve::mcode_cache("pmxforest_secondary_1cmt_oral", .code)

## ---- covariate-specific parameters and regimen ---------------------------
.amt <- dose * FREL          # bioavailable amount; FREL comes from $PK
.mod <- mrgsolve::param(.mod, CL = CL, V = V, KA = KA)
.ev  <- mrgsolve::ev(amt = .amt, ii = tau, addl = n - 1, cmt = "ABS")

## mrgsim_df() returns a plain data.frame. Do NOT use as.data.frame() on a
## mrgsim() result here: `mrgsims` is an S4 class and its as.data.frame method
## only dispatches when mrgsolve is *attached* (library(mrgsolve)); this file
## runs with mrgsolve merely loaded via ::, so as.data.frame() would fall
## through to the default method and error.
.sim <- mrgsolve::mrgsim_df(
  .mod, events = .ev, end = n * tau, delta = 0.1
)

## ---- Cmax within the final (steady-state) dosing interval ---------------
.lastInterval <- subset(.sim, time >= (n - 1) * tau)
max(.lastInterval$CP)
