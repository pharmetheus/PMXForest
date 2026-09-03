# Secondary parameter: steady-state Cmax after q24h oral dosing, via mrgsolve.
#
# This file is inlined by createParamFunction(secondary = list(CMAX = "<this>"))
# into the body of the generated parameter function, wrapped in local({ }). It
# therefore sees, by name:
#   * thetas            - the NONMEM THETA vector
#   * df                - the one-row covariate data frame (columns as df$NAME)
#   * ...               - anything getForestDFSCM() forwards
#   * CL, V, MAT, FREL, KA, D1, ...  - the structural $PK parameters, already
#                         computed above in the generated function
# The value of the LAST expression below becomes CMAX.

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
.tau   <- 24                 # dosing interval (h)
.ndose <- 7                  # doses - enough to reach steady state here
.amt   <- 80 * FREL          # bioavailable amount; FREL comes from $PK

.mod <- mrgsolve::param(.mod, CL = CL, V = V, KA = KA)
.ev  <- mrgsolve::ev(amt = .amt, ii = .tau, addl = .ndose - 1, cmt = "ABS")

.sim <- as.data.frame(mrgsolve::mrgsim(
  .mod, events = .ev, end = .ndose * .tau, delta = 0.1
))

## ---- Cmax within the final (steady-state) dosing interval ---------------
.lastInterval <- subset(.sim, time >= (.ndose - 1) * .tau)
max(.lastInterval$CP)
