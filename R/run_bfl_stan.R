#' Run a BFL Stan model and extract posterior draws
#'
#' Internal helper used by \code{run_BFL()}. Selects the appropriate Stan program
#' for the requested variant, runs MCMC via \pkg{rstan}, and returns posterior draws.
#'
#' @param stan_data Named list matching the Stan data block.
#' @param stan_args List with Stan sampling args (iter, chains, seed, ...).
#' @param variant One of "no_partial", "balanced", "unbalanced".
#'
#' @return A list with posterior draws (pi, lambda, and optionally pi_O) and the Stan fit.
#'
#' @keywords internal
run_bfl_stan <- function(
    stan_data,
    stan_args,
    variant = c("no_partial", "balanced", "unbalanced")
) {

  variant <- match.arg(variant)

  # Select Stan file by variant
  stan_file <- switch(
    variant,
    no_partial = "no_partial_labels.stan",
    balanced   = "partial_labels_shift_false.stan",
    unbalanced = "partial_labels_shift_true.stan"
  )

  stan_path <- system.file("stan", stan_file, package = "BFL")
  if (stan_path == "") stop("Stan file not found: ", stan_file)

  fit <- rstan::sampling(
    rstan::stan_model(stan_path),
    data = stan_data,
    iter = stan_args$iter,
    chains = stan_args$chains,
    seed = stan_args$seed,
    refresh = 0
  )

  post <- rstan::extract(fit)

  out <- list(
    pi = post$pi,
    lambda = post$lambda,
    stan_fit = fit
  )

  if (variant == "unbalanced") {
    out$pi_O <- post$pi_O
  }

  out
}
