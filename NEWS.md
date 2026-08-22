# APRScenario 0.0.4.2

## Bug fixes

* Vignette: removed the call to `bsvars::forecast()` and the version-dispatch
  helper introduced in 0.0.4.1. The forecast object was never used -- the
  vignette computes its unconditional forecast with `forc_h()` -- so the
  vignette no longer depends on which package owns the `forecast()` generic.
  No exported function of APRScenario calls `forecast()`. The vignette builds
  under both `bsvarSIGNs` 2.0 and 3.0.

* Vignette: removed a `dev.new()` / `grid::grid.newpage()` pair that opened a
  graphics device during vignette building, and a
  `Sys.setlocale("LC_TIME", "English.UTF-8")` call using a Windows-only locale
  name that failed silently on other platforms.

## Acknowledgements

* Thanks to PhoenixEmik (https://github.com/PhoenixEmik) for reporting the
  `bsvarSIGNs` 3.0 breakage and independently diagnosing its cause
  (issue #2, pull request #3).

# APRScenario 0.0.4.1

## Bug fixes

* Updated the vignette to use the `generics::forecast()` generic exported by
  `bsvarSIGNs` 3.0 while retaining compatibility with older `bsvarSIGNs`
  versions, which dispatch forecasts through `bsvars::forecast()`.

# APRScenario 0.0.4.0

## Bug fixes

* Fixed a critical transcription error in the moving-average (IRF) recursion of
  `mat_forc()` / `big_b_and_M()`: the code computed `M_i = sum_j M_j B_j`
  instead of `M_i = sum_j M_{i-j} B_j` (Antolin-Diaz, Petrella and
  Rubio-Ramirez, 2021, Appendix A). Impact and one-step-ahead responses were
  unaffected, but all IRF matrices two or more steps after the shock were
  incorrect. Consequently, results from `scenarios()`, `forc_h()` and `KL()`
  at forecast horizons of 3 or more periods change with this release and
  should be recomputed. The corrected matrices now reproduce
  `bsvars::compute_impulse_responses()` at machine precision (see the new
  regression tests).

* Fixed the shock-to-horizon pairing in the simulation loop of `forc_h()`
  (the shock dated `t+tt` was multiplied by the IRF of `tt-1` steps instead of
  `h-tt` steps). Per-horizon forecast quantiles were unaffected in
  distribution; the joint distribution of simulated paths across horizons is
  now coherent.

## New features

* Cross-platform parallelization, contributed by Tito Quadri. `big_b_and_M()`
  and `scenarios()` gain `n_cores` and `parallel` arguments; work is
  distributed over posterior draws using forked processes on Unix/macOS or a
  PSOCK cluster on Windows (`parallel = "auto"` picks the backend). The
  computation is also restructured to build all horizons in a single pass per
  draw, which makes it substantially faster even in the default serial mode.

* `forc_h()` restructured like `scenarios()`: forecast matrices are built in a
  single pass per draw and the shock simulation uses batched per-draw matrix
  products, both parallelized over posterior draws on all platforms
  (roughly 15x faster serially, more with several workers). Given the same
  seed, results are identical for any `max_cores` value and backend.
* Windows compatibility of the remaining fork-based steps: `gen_mats()` and
  `KL()` now fall back to serial execution on Windows instead of erroring
  when more than one core is requested (these steps are computationally
  light; the heavy draw-level computation in
  `big_b_and_M()`/`scenarios()`/`forc_h()` is parallelized on all platforms
  via PSOCK clusters).
* `DESCRIPTION` gains `URL` and `BugReports` fields.

## Authors

* Tito Quadri joins as package author.
