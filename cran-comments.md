# APRScenario 0.0.5.0 — submission comments

## Reason for this release

New exported function `wz_posterior()`. It computes the posterior of the model
parameters updated by the scenario conditions (Waggoner and Zha, 1999;
Antolin-Diaz, Petrella and Rubio-Ramirez, 2021, Algorithm 1) from existing
posterior draws, together with effective-sample-size and Kullback-Leibler
diagnostics. The vignette has a new section illustrating it. No existing
function changed behaviour.

## Test environments

* Local: Ubuntu 26.04.1 LTS, R 4.5.2, gcc 15.2.0,
  `R CMD check --as-cran --run-donttest`, with the suggested packages
  `bsvars` 4.0 and `bsvarSIGNs` 3.0.

## R CMD check results

0 ERRORs, 0 WARNINGs, 2 NOTEs.

Both NOTEs are specific to the local check machine:

* "unable to verify current time" — no network time service is available
  there.
* "Compilation used the following non-portable flag(s):
  -mno-omit-leaf-frame-pointer" — this flag is injected by the local R
  `Makeconf`, not by the package's `src/Makevars`.

## Downstream dependencies

None (no packages depend on APRScenario).
