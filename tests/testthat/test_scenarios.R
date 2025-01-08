# test_scenarios.R

test_that("get_mat works correctly", {
  h <- 4
  gen_mats()
  expect_true(is.array(M))
})

test_that("forc_h works correctly", {
  gen_mats()
  h <- 4
  y_h_all <- forc_h(h, n_sim = 200)
  expect_true(is.array(y_h_all[[1]]))
})

test_that("scenarios returns correct array", {
  gen_mats()
  h <- 4
  n_sim <- 200
  obs <- c(2)
  TT <- nrow(X0)
  path <- X0[(TT - h + 1):TT, obs]
  bvarSign_path <- X0[(TT - h + 1):TT, ]
  bvarSign_path[,-obs] <- NA
  shocks <- NA
  tmp <- scenarios(h, path, obs, shocks)
  mu_eps <- tmp[[1]]
  expect_true(is.array(mu_eps))
})

test_that("SimScen generates correct array", {
  gen_mats()
  h <- 4
  n_sim <- 200
  obs <- c(2)
  TT <- nrow(X0)
  path <- X0[(TT - h + 1):TT, obs]
  bvarSign_path <- X0[(TT - h + 1):TT, ]
  bvarSign_path[,-obs] <- NA
  shocks <- NA
  tmp <- scenarios(h, path, obs, shocks)
  mu_eps <- tmp[[1]]
  Sigma_eps <- tmp[[2]]
  mu_y <- tmp[[3]]
  Sigma_y <- tmp[[4]]
  big_b <- tmp[[5]]
  big_M <- tmp[[6]]
  y_h <- SimScen(mu_eps, Sigma_eps, mu_y, Sigma_y, big_b, big_M, n_sim,h)
  expect_true(is.array(y_h))
})

test_that("KL returns correct matrix", {
  h <- 4
  gen_mats()
  y_h_all <- forc_h(h, n_sim = 200)
  n_sim <- 200
  obs <- c(2)
  TT <- nrow(X0)
  path <- X0[(TT - h + 1):TT, obs]
  bvarSign_path <- X0[(TT - h + 1):TT, ]
  bvarSign_path[,-obs] <- NA
  shocks <- NA
  tmp <- scenarios(h, path, obs, shocks)
  mu_eps <- tmp[[1]]
  Sigma_eps <- tmp[[2]]
  q <- KL(Sigma_eps, mu_eps,h)
  expect_true(is.matrix(q))
})


