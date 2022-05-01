This directory contains the files that are the basis for the PMXForest vignettes and examples. The actual runs were made in PMX-REP-PMX-2/Analysis/Model/SimVal, copied to PMX-REP-PMX-2/Analysis/Model/SimVal-for-Forest and then copied to this location.

run7.mod was developed using scmplus and stage-wise filtering.

run22-3.mod is a frem model based on the model with the mechanistic and structural covariates in run7.mod.

bs7.dir and bs22-3.dir are boostrap runs for the two models with 175 samples.

bs7-n30.dir and bs22-3n30.dir are boostrap runs for the two models with 30 samples, suitable for the small bootstrap approach to uncertainty.

sir7.dir is a sir analysis for run7.mod. Sir didn't work for 22-3, at least not with the setting used.
