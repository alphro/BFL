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
#'   \item{\code{"balanced"}}{Partial labels, single prevalence pi. Calls
#'     \code{gibbs_bfl_cpp(single_pi = TRUE)}.}
#'   \item{\code{"unbalanced"}}{Partial labels, separate pi (unlabeled) and
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
