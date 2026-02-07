#' Run global BFL aggregation
#'
#' Aggregates local LCVA summaries into a global BFL model using Stan.
#' This function never sees individual-level data.
#'
#' @param local_summaries Named list of local summary objects
#' @param stan_args List of arguments passed to rstan::sampling
#'
#' @export
run_BFL <- function(
    local_summaries,
    stan_args = list(iter = 2000, chains = 4, seed = 12345)
) {

  validate_local_summaries(local_summaries)

  aligned <- align_local_summaries(local_summaries)

  stan_data <- build_bfl_stan_data(aligned)

  fit <- run_bfl_stan(stan_data, stan_args)

  # Build a global phi (N x C) by weighting each site's phi
  phi_list_na <- aligned$aligned_phi
  phi_list <- lapply(phi_list_na, function(mat) { mat[is.na(mat)] <- 0; mat })

  N <- nrow(phi_list[[1]])
  C <- ncol(phi_list[[1]])
  M <- length(phi_list)

  # posterior-mean lambda: dims [iter, C, M] -> C x M
  w <- apply(fit$lambda, c(2, 3), mean)
  w <- w / rowSums(w)

  phi_global <- matrix(0, nrow = N, ncol = C)
  for (m in seq_len(M)) {
    phi_global <- phi_global + phi_list[[m]] * matrix(rep(w[, m], each = N), nrow = N, ncol = C)
  }

  list(
    pi = fit$pi,
    lambda = fit$lambda,
    phi = phi_global,
    causes = aligned$global_causes,
    stan_fit = fit$stan_fit
  )
}
