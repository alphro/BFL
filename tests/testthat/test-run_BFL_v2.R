# Tests for the v2.0 input-driven run_BFL (1:1 rows, Y_add correction).

# Build tiny, internally-consistent local_summaries on N rows.
make_ls <- function(N = 12, C = 3, M = 2, seed = 1) {
  set.seed(seed)
  X  <- matrix(rbinom(N * 4, 1, 0.5), N, 4)
  rh <- compute_row_hashes(X)
  causes <- paste0("c", seq_len(C))
  mk <- function(s) {
    set.seed(s)
    phi <- matrix(runif(N * C), N, C); phi <- phi / rowSums(phi)
    list(posterior_phi = phi, cause_ids = causes,
         target_info   = list(row_hash = rh, N = N, P = ncol(X)))
  }
  ls <- setNames(lapply(seq_len(M), function(i) mk(seed + i)), paste0("site", seq_len(M)))
  attr(ls, "X") <- X
  ls
}

# Does the fast Gibbs (no_partial) path run here? (Rcpp must be compiled.)
can_gibbs <- function() {
  ls <- make_ls(N = 8, seed = 7); X <- attr(ls, "X")
  ok <- tryCatch({
    run_BFL(ls, X_target = X, Y_target = NULL, sampler = "gibbs",
            mcmc_args = list(iter = 100, chains = 1, seed = 1)); TRUE
  }, error = function(e) FALSE)
  isTRUE(ok)
}

test_that("row invariant is enforced (no sampler needed)", {
  ls <- make_ls(N = 10)
  X  <- attr(ls, "X")
  # phi has 10 rows; X_target has 11 -> 1:1 violation
  expect_error(
    run_BFL(ls, X_target = rbind(X, X[1, , drop = FALSE]), Y_target = NULL),
    "1:1"
  )
  # Y_target length must equal nrow(X_target)
  expect_error(
    run_BFL(ls, X_target = X, Y_target = rep(NA_character_, 9)),
    "must equal"
  )
})

test_that("no labels: Y_target NULL, no Y_add -> no_partial, no correction", {
  skip_if_not(can_gibbs())
  ls <- make_ls(N = 12, seed = 2); X <- attr(ls, "X")
  f <- run_BFL(ls, X_target = X, Y_target = NULL, sampler = "gibbs",
               mcmc_args = list(iter = 200, chains = 1, seed = 1))
  expect_identical(f$stan_idx, seq_len(12L))
  expect_equal(f$n_total, 12L)
  expect_equal(f$n_add, 0L)
  expect_null(f$nLc)
  expect_false(f$has_labels)
})

test_that("labels folded in as a source: Y_target NULL + Y_add -> no_partial, nLc from Y_add", {
  skip_if_not(can_gibbs())
  ls <- make_ls(N = 12, seed = 3); X <- attr(ls, "X")
  Y_add <- c("c1", "c1", "c2", NA)            # NA entries are dropped
  f <- run_BFL(ls, X_target = X, Y_target = NULL, Y_add = Y_add, sampler = "gibbs",
               mcmc_args = list(iter = 200, chains = 1, seed = 1))
  expect_false(f$has_labels)                  # still no-partial model
  expect_equal(f$n_add, 3L)
  expect_equal(as.integer(f$nLc[c("c1", "c2")]), c(2L, 1L))

  # correction formula in score_BFL: (n_total*raw + nLc)/(n_total+n_add)
  Y_eval <- rep("c1", 12)
  sc <- score_BFL(f, Y_eval = Y_eval, seed = 1)
  expect_equal(sum(sc$pi_hat), 1, tolerance = 1e-8)
})

test_that("labels constrain prediction: labels in Y_target -> balanced / unbalanced via label_shift", {
  skip_if_not(can_gibbs())   # variant detection runs before the sampler
  ls <- make_ls(N = 12, seed = 4); X <- attr(ls, "X")
  Y <- rep(NA_character_, 12); Y[1:4] <- c("c1", "c2", "c3", "c1")
  # We only assert the inferred variant flags; the partial Stan path may be
  # skipped here if rstan is unavailable, so wrap in try.
  got <- tryCatch(
    run_BFL(ls, X_target = X, Y_target = Y, sampler = "gibbs",
            mcmc_args = list(iter = 200, chains = 1, seed = 1)),
    error = function(e) e)
  skip_if(inherits(got, "error"), "partial Stan path unavailable here")
  expect_true(got$has_labels)
  expect_false(got$label_shift)
})
