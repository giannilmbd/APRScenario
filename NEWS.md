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
