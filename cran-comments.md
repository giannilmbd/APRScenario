# APRScenario 0.0.4.0 — resubmission comments

## Reason for this release

This is a **critical bug-fix release**. The moving-average (impulse-response)
recursion in `mat_forc()`/`big_b_and_M()` was mistranscribed from Appendix A of
Antolin-Diaz, Petrella and Rubio-Ramirez (2021, JME): the code computed
`M_i = sum_j M_j B_j` instead of `M_i = sum_j M_{i-j} B_j`. Impulse responses
and all scenario/conditional-forecast output at horizons >= 3 produced by
versions <= 0.0.3.1 were incorrect. The corrected matrices reproduce
`bsvars::compute_impulse_responses()` at machine precision; regression tests
pinning this are included (`tests/testthat/test_irf_recursion.R`).

The release also adds cross-platform parallelization over posterior draws
(forked processes on Unix/macOS, PSOCK clusters on Windows), contributed by
Tito Quadri, who joins as package author. All parallel defaults are serial;
tests use at most 2 workers.

## Test environments

* Local: Ubuntu 24.04, R 4.4.x (gcc/g++ 13), `R CMD check --as-cran --run-donttest`
* win-builder (release and devel)

## R CMD check results

0 ERRORs, 0 WARNINGs on machines with `qpdf` available.

Local NOTEs, both environment-specific:

* "unable to verify current time" — no network time service available on the
  check machine.
* "Compilation used the following non-portable flag(s): -mno-omit-leaf-frame-pointer"
  — flag injected by the local R Makeconf, not by the package's Makevars.

## Downstream dependencies

None (no packages depend on APRScenario).
