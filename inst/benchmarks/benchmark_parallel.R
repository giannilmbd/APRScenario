# Benchmark of the parallelized big_b_and_M()/scenarios() (contributed by Tito Quadri).
# Not run by R CMD check; run interactively.
#
# Findings that motivated the design (originally on Windows, 8 cores):
# - big_b_and_M() is the computational bottleneck of scenarios();
# - the per-draw single-pass restructuring is faster than the old per-horizon
#   code even with 1 core;
# - parallelization over draws scales well for large n_draws/p.

library(APRScenario)

data(NKdata)

S <- 1000 # posterior draws
h <- 24   # forecast horizon
p <- 10   # VAR lags

spec <- bsvarSIGNs::specify_bsvarSIGN$new(as.matrix(NKdata[, 2:4]), p = p)
posterior <- bsvars::estimate(spec, S = S)
matrices <- gen_mats(posterior = posterior, specification = spec)

n_var   <- dim(posterior$posterior$B)[1]
n_p     <- (dim(posterior$posterior$A)[2] - 1) / n_var
n_draws <- dim(posterior$posterior$B)[3]

t1 <- system.time(
  serial <- big_b_and_M(h = h, n_draws = n_draws, n_var = n_var, n_p = n_p,
                        matrices = matrices, n_cores = 1)
)

n_workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
t2 <- system.time(
  par_auto <- big_b_and_M(h = h, n_draws = n_draws, n_var = n_var, n_p = n_p,
                          matrices = matrices, n_cores = n_workers)
)

cat(sprintf("big_b_and_M: serial %.2fs | %d workers (auto backend) %.2fs (x%.1f)\n",
            t1["elapsed"], n_workers, t2["elapsed"], t1["elapsed"] / t2["elapsed"]))

# equality check between serial and parallel results
stopifnot(max(abs(unlist(serial) - unlist(par_auto))) < 1e-12)

# whole scenarios() call
t3 <- system.time(
  scen <- scenarios(h = h, path = rep(1, h), obs = 1, free_shocks = NA,
                    posterior = posterior, matrices = matrices, n_cores = n_workers)
)
cat(sprintf("scenarios() with %d workers: %.2fs\n", n_workers, t3["elapsed"]))
