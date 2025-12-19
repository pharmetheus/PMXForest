#' @importFrom utils globalVariables
NULL

# This file defines global variables to satisfy R CMD check when using
# non-standard evaluation in dplyr, ggplot2, and other packages.

utils::globalVariables(c(
  "COVEFF", "COVNAME", "COVNUM", "GROUPNAME", "GROUPNAME2", "GROUPNAMELABEL",
  "ITER", "ITERATION", "NAME", "PARAMETER", "PARAMETERLABEL", "REF",
  "STATISTIC", "STATISTICSLABEL", "SUBJ", "THETA1", "lowcilabel",
  "upcilabel", "meanlabel", "ofv", "point", "q1", "q2", "resamples",
  "xmax", "xmin", "ymax", "ymin"
))
