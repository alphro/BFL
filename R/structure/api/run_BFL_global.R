#' Run global BFL aggregation
#'
#' Aggregates local LCVA summaries into a global BFL model using Stan.
#' This function never sees individual-level data.
#'
#' @param local_summaries Named list of local summary objects
#' @param stan_args List of arguments passed to rstan::sampling
#'
#' @export
run_BFL_global <- function(
    local_summaries,
    stan_args = list(iter = 2000, chains = 4, seed = 12345)
) {

  validate_local_summaries(local_summaries)

  aligned <- align_local_summaries(local_summaries)

  stan_data <- build_bfl_stan_data(aligned)

  fit <- run_bfl_stan(stan_data, stan_args)

  phi_list <- aligned$aligned_phi  # list of P x C matrices, one per model/site
  w <- colMeans(fit$lambda)
  w <- w / sum(w)

  phi_global <- Reduce(`+`, Map(function(phi_m, wm) wm * phi_m, phi_list, w))
  phi_global <- pmin(pmax(phi_global, 1e-12), 1 - 1e-12)

  list(
    pi = fit$pi,
    lambda = fit$lambda,
    phi = phi_global,
    causes = aligned$global_causes,
    stan_fit = fit$stan_fit
  )
}
