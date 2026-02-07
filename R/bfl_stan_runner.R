# PURPOSE:
# Run a BFL Stan model (no-partial, balanced, or unbalanced)
# and extract posterior summaries.

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
    balanced   = "partial_labels_false.stan",
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
    pi = colMeans(post$pi),
    lambda = post$lambda,
    stan_fit = fit
  )

  if (variant == "unbalanced") {
    out$pi_O <- colMeans(post$pi_O)
  }

  out
}
