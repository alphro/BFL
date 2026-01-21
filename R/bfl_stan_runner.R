# PURPOSE:
# Runs the BFL Stan model and extracts posterior summaries.

run_bfl_stan <- function(stan_data, stan_args) {

  stan_path <- system.file("stan", "no_partial_labels.stan", package = "BFL")
  if (stan_path == "") stop("Stan file not found.")

  model <- rstan::stan_model(stan_path)

  fit <- rstan::sampling(
    model,
    data = stan_data,
    iter = stan_args$iter,
    chains = stan_args$chains,
    seed = stan_args$seed,
    refresh = 0
  )

  post <- rstan::extract(fit)

  list(
    pi = colMeans(post$pi),
    lambda = post$lambda,
    stan_fit = fit
  )
}
