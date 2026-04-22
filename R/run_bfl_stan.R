#' Run a BFL model and extract posterior draws
#'
#' Internal helper used by \code{run_BFL()}. Dispatches to the Gibbs sampler
#' or Stan depending on \code{sampler} and \code{variant}.
#'
#' @param stan_data Named list matching the Stan data block.
#' @param mcmc_args List of shared MCMC controls: iter, chains, seed, and
#'   optionally init (Stan only, default "random").
#' @param gibbs_args List of Gibbs-specific options: logistic_normal, mh_scale.
#'   Ignored when sampler = "stan".
#' @param variant One of "no_partial", "balanced", "unbalanced".
#' @param sampler One of "gibbs" or "stan". Gibbs is only available for
#'   "no_partial"; partial-label variants always fall through to Stan.
#'
#' @return A list with posterior draws (pi, lambda, and optionally pi_O) and
#'   the Stan fit object (NULL for Gibbs).
#'
#' @keywords internal
run_bfl_stan <- function(
    stan_data,
    mcmc_args  = list(),
    gibbs_args = list(),
    variant    = c("no_partial", "balanced", "unbalanced"),
    sampler    = c("gibbs", "stan")
) {
  variant <- match.arg(variant)
  sampler <- match.arg(sampler)

  # Gibbs sampler path: all three variants supported
  if (sampler == "gibbs") {
    return(run_bfl_gibbs(stan_data, mcmc_args, gibbs_args, variant = variant))
  }

  # Stan path
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
    data    = stan_data,
    iter    = mcmc_args$iter,
    chains  = mcmc_args$chains,
    seed    = mcmc_args$seed,
    init    = if (is.null(mcmc_args$init)) "random" else mcmc_args$init,
    refresh = 100
  )

  post <- rstan::extract(fit)

  out <- list(pi = post$pi, lambda = post$lambda, stan_fit = fit)
  if (variant == "unbalanced") out$pi_O <- post$pi_O
  out
}
