#' Run the no-partial-labels BFL model via Rcpp Gibbs sampler
#'
#' Drop-in replacement for the Stan path inside \code{run_bfl_stan()} when
#' \code{variant = "no_partial"}.  Runs one independent chain per element of
#' \code{1:mcmc_args$chains}, each seeded at \code{mcmc_args$seed + chain},
#' then concatenates the post-warmup draws into the same output format that
#' \code{rstan::extract()} would produce.
#'
#' @section Prior for lambda_k:
#' Controlled by \code{gibbs_args$logistic_normal}:
#' \describe{
#'   \item{\code{FALSE} (default)}{Conjugate \code{Dir(1,...,1)} prior.
#'     Full conditional is closed-form; no tuning required.}
#'   \item{\code{TRUE}}{Logistic-normal prior matching Stan:
#'     \code{beta_raw[k,m] ~ N(0,1)}, \code{lambda_k = softmax(beta_eff)}.
#'     Active beta components are updated via random-walk MH with step size
#'     \code{gibbs_args$mh_scale} (default 0.5). Target acceptance rate is
#'     roughly 0.3-0.5; decrease \code{mh_scale} if acceptance is too low.
#'     Acceptance rates are printed per chain.}
#' }
#'
#' @param stan_data  Named list matching the Stan data block.
#' @param mcmc_args  List of shared MCMC controls: iter (default 2000),
#'   chains (default 4), seed (default 12345).
#' @param gibbs_args List of Gibbs-specific options: logistic_normal (default
#'   FALSE), mh_scale (default 0.5).
#'
#' @return List with:
#'   \describe{
#'     \item{pi}{Matrix \[S x num_causes\] of posterior cause-fraction draws.}
#'     \item{lambda}{Array \[S x num_causes x M\] of posterior model-weight draws.}
#'     \item{stan_fit}{NULL.}
#'   }
#' @keywords internal
run_bfl_gibbs <- function(stan_data, mcmc_args = list(), gibbs_args = list()) {

  n_iter          <- mcmc_args$iter   %||% 2000L
  n_chains        <- mcmc_args$chains %||% 4L
  seed            <- mcmc_args$seed   %||% 12345L
  logistic_normal <- gibbs_args$logistic_normal %||% FALSE
  mh_scale        <- gibbs_args$mh_scale        %||% 0.5

  n_iter   <- as.integer(n_iter)
  n_chains <- as.integer(n_chains)
  n_warmup <- n_iter %/% 2L
  n_keep   <- n_iter - n_warmup

  K       <- stan_data$num_causes
  M       <- stan_data$M
  phi_vec <- as.vector(stan_data$phi)  # column-major flatten

  chain_results <- vector("list", n_chains)
  for (ch in seq_len(n_chains)) {
    set.seed(seed + ch)
    chain_results[[ch]] <- gibbs_no_partial_cpp(
      phi             = phi_vec,
      model_presence  = stan_data$model_presence,
      causes_r        = stan_data$causes,
      N               = stan_data$N,
      C_max           = stan_data$C_max,
      M               = M,
      num_causes      = K,
      n_iter          = n_iter,
      n_warmup        = n_warmup,
      logistic_normal = logistic_normal,
      mh_scale        = mh_scale
    )
    if (logistic_normal) {
      rates <- round(chain_results[[ch]]$mh_accept_rate, 3)
      message("Chain ", ch, " MH acceptance rates (per cause): ",
              paste(rates, collapse = ", "))
    }
  }

  # Combine chains: pi -> [S_total x K], lambda -> [S_total x K x M]
  pi_all     <- do.call(rbind, lapply(chain_results, `[[`, "pi"))
  S_total    <- n_keep * n_chains
  lambda_all <- array(0.0, dim = c(S_total, K, M))
  for (ch in seq_len(n_chains)) {
    rows <- seq_len(n_keep) + (ch - 1L) * n_keep
    lambda_all[rows, , ] <- array(chain_results[[ch]]$lambda,
                                  dim = c(n_keep, K, M))
  }

  list(pi = pi_all, lambda = lambda_all, stan_fit = NULL)
}
