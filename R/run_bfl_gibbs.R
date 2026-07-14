#' Run a BFL model via Rcpp Gibbs sampler
#'
#' Dispatcher for all three BFL variants. Runs one independent Gibbs chain per
#' element of \code{1:mcmc_args$chains}, seeds each at \code{mcmc_args$seed +
#' chain}, and concatenates post-warmup draws into the same format that
#' \code{rstan::extract()} would produce.
#'
#' @section Variant behaviour:
#' \describe{
#'   \item{\code{"no_partial"}}{All causes unknown. Calls
#'     \code{gibbs_bfl_cpp} with all-zero Y_known.}
#'   \item{\code{"balanced"}}{With partial labels, a single prevalence pi. Calls
#'     \code{gibbs_bfl_cpp(single_pi = TRUE)}.}
#'   \item{\code{"unbalanced"}}{With partial labels, separate pi (unlabeled) and
#'     pi_O (labeled). Calls \code{gibbs_bfl_cpp(single_pi = FALSE)}.
#'     Returns \code{pi_O} in addition to \code{pi} and \code{lambda}.}
#' }
#'
#' @section Prior for lambda_k:
#' Controlled by \code{gibbs_args$logistic_normal}:
#' \describe{
#'   \item{\code{FALSE} (default)}{Conjugate \code{Dir(1,...,1)} — all steps
#'     are closed-form.}
#'   \item{\code{TRUE}}{Logistic-normal matching Stan's prior. Active beta
#'     components updated via random-walk MH with step size
#'     \code{gibbs_args$mh_scale} (default 0.5). Acceptance rates printed per
#'     chain.}
#' }
#'
#' @param stan_data  Named list matching the Stan data block.
#' @param mcmc_args  List: iter (default 2000), chains (default 4),
#'   seed (default 12345).
#' @param gibbs_args List: logistic_normal (default FALSE), mh_scale
#'   (default 0.5).
#' @param variant    One of \code{"no_partial"}, \code{"balanced"},
#'   \code{"unbalanced"}.
#'
#' @return List: pi [S x K], lambda [S x K x M], stan_fit (NULL), and pi_O
#'   [S x K] for the unbalanced variant.
#' @keywords internal
run_bfl_gibbs <- function(stan_data,
                          mcmc_args  = list(),
                          gibbs_args = list(),
                          variant    = c("no_partial", "balanced", "unbalanced")) {

  variant         <- match.arg(variant)
  n_iter          <- as.integer(mcmc_args$iter   %||% 2000L)
  n_chains        <- as.integer(mcmc_args$chains %||% 4L)
  seed            <- mcmc_args$seed              %||% 12345L
  logistic_normal <- gibbs_args$logistic_normal  %||% FALSE
  mh_scale        <- gibbs_args$mh_scale         %||% 0.5

  n_warmup <- n_iter %/% 2L
  n_keep   <- n_iter - n_warmup
  K        <- stan_data$num_causes
  M        <- stan_data$M
  phi_vec  <- as.vector(stan_data$phi)
  N        <- stan_data$N

  # Build Y arguments: no_partial uses sentinel all-zero Y_known
  if (variant == "no_partial") {
    Y_known <- integer(N)
    Y_idx_r <- rep(1L, N)
    single_pi <- TRUE
  } else {
    Y_known   <- stan_data$Y_known
    Y_idx_r   <- stan_data$Y_idx
    single_pi <- (variant == "balanced")
  }

  # ---- Pre-Gibbs dimension + sanity checks -----------------------------------
  # These stop() before the C++ call so the user sees a clear error message
  # rather than a segfault.
  expected_phi_len <- N * stan_data$C_max * M
  if (length(phi_vec) != expected_phi_len)
    stop(sprintf(
      "run_bfl_gibbs: phi_vec length %d does not match N=%d * C_max=%d * M=%d = %d.\n  Check that stan_data$phi has dim c(N, C_max, M).",
      length(phi_vec), N, stan_data$C_max, M, expected_phi_len))

  if (!identical(dim(stan_data$model_presence), c(K, M)))
    stop(sprintf(
      "run_bfl_gibbs: model_presence has dim [%s], expected [%d x %d] (num_causes x M).",
      paste(dim(stan_data$model_presence), collapse = " x "), K, M))

  if (length(stan_data$causes) != K)
    stop(sprintf(
      "run_bfl_gibbs: length(causes) = %d but num_causes = %d.",
      length(stan_data$causes), K))

  if (max(stan_data$causes) > stan_data$C_max)
    stop(sprintf(
      "run_bfl_gibbs: max(causes) = %d > C_max = %d.  Cause indices exceed the phi array column range.",
      max(stan_data$causes), stan_data$C_max))

  n_bad_phi <- sum(is.na(phi_vec) | is.infinite(phi_vec))
  if (n_bad_phi > 0)
    stop(sprintf("run_bfl_gibbs: phi contains %d NA/Inf value(s).", n_bad_phi))

  if (any(phi_vec < 0, na.rm = TRUE))
    stop(sprintf("run_bfl_gibbs: phi contains %d negative value(s).",
                 sum(phi_vec < 0, na.rm = TRUE)))

  # ---- Distinguish structural absence from numerical underflow ----------------
  #
  # Structural absence: model_presence[k, m] = 0 means source m genuinely does
  #   not have cause k (e.g. an absent cause expanded to zero by
  #   expand_phi_to_causes). These columns ARE expected to be all-zero and must
  #   NOT be reported as underflow.
  #
  # Numerical underflow: model_presence[k, m] = 1 but all phi[i, k, m] = 0.
  #   This means the source claims to have the cause but every observation's
  #   likelihood underflowed to zero.  This is the pathological case that
  #   degrades inference and should be warned.
  # ---------------------------------------------------------------------------

  phi_arr <- array(phi_vec, dim = c(N, stan_data$C_max, M))

  # Structural absence count (informational only — expected and handled correctly).
  n_structural_absent <- sum(stan_data$model_presence == 0L)
  if (n_structural_absent > 0)
    message(sprintf(
      "run_bfl_gibbs: %d structural absent (source,cause) pairs — model_presence=0, phi=0 by design. These are excluded from the likelihood (correct).",
      n_structural_absent))

  # Dead causes: ALL sources structurally absent for cause k.
  # These will never be sampled; posterior is pi-prior only.
  dead_k <- which(rowSums(stan_data$model_presence) == 0L)
  if (length(dead_k) > 0)
    warning(sprintf(
      "run_bfl_gibbs: %d cause(s) have no active source at all (model_presence all-zero row): cause indices [%s].\n  These causes cannot be sampled; their posterior is prior-only.\n  This is expected if all %d sources lack those causes entirely.",
      length(dead_k), paste(dead_k, collapse = ", "), M))

  # True underflow: model_presence=1 but all phi rows are zero for that (cause, source).
  underflow_count <- 0L
  underflow_pairs <- character(0)
  for (m in seq_len(M)) {
    for (k in seq_len(K)) {
      if (stan_data$model_presence[k, m] == 1L) {
        c_idx <- stan_data$causes[k]   # 1-indexed column in phi array
        col_sum <- sum(phi_arr[, c_idx, m])
        if (col_sum == 0) {
          underflow_count <- underflow_count + 1L
          underflow_pairs <- c(underflow_pairs,
                               sprintf("(cause=%d,src=%d)", k, m))
        }
      }
    }
  }
  if (underflow_count > 0)
    warning(sprintf(
      "run_bfl_gibbs: %d (cause,source) pair(s) have model_presence=1 but all-zero phi (true numerical underflow): %s.\n  Inference for those cause/source combinations will use only other active sources.",
      underflow_count,
      paste(head(underflow_pairs, 10), collapse = "; ")))
  # ---------------------------------------------------------------------------

  chain_results <- vector("list", n_chains)

  for (ch in seq_len(n_chains)) {
    set.seed(seed + ch)
    chain_results[[ch]] <- gibbs_bfl_cpp(
      phi             = phi_vec,
      model_presence  = stan_data$model_presence,
      causes_r        = stan_data$causes,
      Y_known         = Y_known,
      Y_idx_r         = Y_idx_r,
      N               = N,
      C_max           = stan_data$C_max,
      M               = M,
      num_causes      = K,
      n_iter          = n_iter,
      n_warmup        = n_warmup,
      single_pi       = single_pi,
      logistic_normal = logistic_normal,
      mh_scale        = mh_scale
    )
    .report_mh(chain_results[[ch]], ch, logistic_normal)
  }

  pi_all     <- .combine_pi(chain_results, n_keep, n_chains)
  lambda_all <- .combine_lambda(chain_results, n_keep, n_chains, K, M)
  out <- list(pi = pi_all, lambda = lambda_all, stan_fit = NULL)

  if (!single_pi) {
    S_total  <- n_keep * n_chains
    pi_O_all <- matrix(0.0, S_total, K)
    for (ch in seq_len(n_chains)) {
      rows <- seq_len(n_keep) + (ch - 1L) * n_keep
      pi_O_all[rows, ] <- chain_results[[ch]]$pi_O
    }
    out$pi_O <- pi_O_all
  }

  out
}

# ---- Private helpers -------------------------------------------------------

.combine_pi <- function(chain_results, n_keep, n_chains) {
  do.call(rbind, lapply(chain_results, `[[`, "pi"))
}

.combine_lambda <- function(chain_results, n_keep, n_chains, K, M) {
  S_total    <- n_keep * n_chains
  lambda_all <- array(0.0, dim = c(S_total, K, M))
  for (ch in seq_len(n_chains)) {
    rows <- seq_len(n_keep) + (ch - 1L) * n_keep
    lambda_all[rows, , ] <- array(chain_results[[ch]]$lambda,
                                  dim = c(n_keep, K, M))
  }
  lambda_all
}

.report_mh <- function(result, ch, logistic_normal) {
  if (!logistic_normal) return(invisible(NULL))
  rates <- result$mh_accept_rate
  message(sprintf(
    "Chain %d MH acceptance: mean = %.3f, range = [%.3f, %.3f]",
    ch, mean(rates, na.rm = TRUE), min(rates, na.rm = TRUE), max(rates, na.rm = TRUE)
  ))
}
