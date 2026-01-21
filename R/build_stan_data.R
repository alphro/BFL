build_bfl_stan_data <- function(aligned) {

  phi_list <- aligned$aligned_phi
  M <- length(phi_list)
  N <- nrow(phi_list[[1]])
  C <- ncol(phi_list[[1]])  # safer than length(global_causes)

  # phi: N x C x M
  phi_array <- array(0, dim = c(N, C, M))
  for (m in seq_len(M)) {
    stopifnot(nrow(phi_list[[m]]) == N, ncol(phi_list[[m]]) == C)
    phi_array[, , m] <- phi_list[[m]]
  }

  # model_presence: C x M (cause present in model if column not all-zero)
  model_presence <- matrix(0L, nrow = C, ncol = M)
  for (m in seq_len(M)) {
    present <- colSums(phi_list[[m]] != 0) > 0
    model_presence[, m] <- as.integer(present)
  }

  list(
    N = as.integer(N),
    C_max = as.integer(C),
    M = as.integer(M),
    num_causes = as.integer(C),

    # IMPORTANT: integers for Stan, not character labels
    causes = as.integer(seq_len(C)),

    model_presence = model_presence,
    count = as.integer(rowSums(model_presence)),
    phi = phi_array
  )
}
