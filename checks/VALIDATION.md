# Validation ledger

Each row is a check: a script plus its log in `checks/`. Tolerances were
written into the script before the first full run. `checks/` is in
`.Rbuildignore`.

| Date | Estimator | Oracle | Tolerance | Result | Script / log |
|---|---|---|---|---|---|
| 2026-11-20 | Forecast moments b, Sigma_y = Phi (I x Sigma) Phi' (conjugate test library) | 200,000 simulated VAR paths from one parameter draw | max z(mean) < 4; max relative covariance error < 0.02 | PASS: 1.22; 0.0032 | `wz/wz_conjugate_validation.R`, `.log` |
| 2026-11-20 | Importance reweighting, hard observable conditions, w ∝ N(f; Cb, C Sy C') | Exact Waggoner-Zha / APR Algorithm 1 Gibbs sampler, conjugate NIW VAR(1), n=2, h=4, T=80, 40,000 iterations | max z < 4 over 34 statistics (parameters, per-draw conditional means and variances) | PASS: 2.29; ESS 71% | same |
| 2026-11-20 | Importance reweighting, noisy-signal density conditions, w ∝ N(f; Cb, C Sy C' + Omega) | Exact Gibbs, noisy-signal model | max z < 4 | PASS: 2.10; ESS 78% | same |
| 2026-11-20 | APR Algorithm 1 with step 1 by resampling 20,000 draws from p(theta \| Y) (SIRG), hard conditions | Exact Gibbs | max z < 4 | PASS: 2.95 | same |
| 2026-11-20 | SIRG, noisy-signal conditions | Exact Gibbs | max z < 4 | PASS: 3.92 (closest to the bound) | same |
| 2026-11-20 | SIRG, APR density conditions, Omega diagonal | Exact Gibbs | max z < 4 | PASS: 2.26 | same |
| 2026-11-20 | SIRG, APR density conditions, Omega = D D' (package default) | Exact Gibbs | max z < 4 | PASS: 2.16 | same |

Non-vacuity, all cases: the exact-Gibbs parameter means differ from the
p(theta | Y) means by z = 65 to 88. The test therefore discriminates between
the two posteriors.

Scope. These checks cover the algebra and sampling of the weights and of
SIRG in a reduced-form conjugate VAR. They do not cover the bsvarSIGNs
adapter (`gen_mats`, `big_b_and_M`), structural (shock) restrictions, or
weight degeneracy in larger systems. Rows for the package functions follow.

## Package functions

| Date | Estimator | Oracle | Tolerance | Result | Script / log |
|---|---|---|---|---|---|
| 2026-11-20 | `gen_mats` -> `big_b_and_M`: forecast means and M'M, 20,000 conjugate draws passed as a bsvarSIGNs-shaped posterior | `wz_conjugate_lib.R` | max relative diff < 1e-8 | PASS: 8.3e-16 | `wz/wz_package_validation.R`, `.log` |
| 2026-11-20 | `scenarios()` per-draw mu_y, Sigma_y, APR density (diagonal Omega and D D') | y-space form of APR eqs. 13-14 in `wz_conjugate_lib.R` | < 1e-8 | PASS: 1.0e-14 | same |
| 2026-11-20 | `wz_posterior` hard and noisy: weights, per-draw moments, kl_param; hard moments vs `scenarios(Sigma_g = 0)`; noisy kl_draw in y vs shock coordinates | oracle library; APR eq. 16 | < 1e-8 | PASS: max 8.7e-15 | same |
| 2026-11-20 | `wz_posterior` hard and noisy, posterior means | exact Gibbs (`wz_conjugate_results.rds`) | max z < 4 | PASS: 1.67, 1.63 | same |
| 2026-11-20 | `wz_posterior` apr (C++ `wz_sirg_core`), diagonal Omega and D D'; standard errors from 25 independent replications of 1,000 iterations | exact Gibbs | max z < 4 | PASS: 2.49, 2.56 | same |
| 2026-11-20 | `wz_posterior` apr, kl_draw under D D' | mu_eps' mu_eps / 2 | < 1e-8 | PASS: 2.9e-16 | same |

Not covered: real bsvarSIGNs output (bsvars and bsvarSIGNs are not installed on
the machine that ran these checks), structural scenarios with shock
restrictions, and systems larger than n = 2, h = 4.
