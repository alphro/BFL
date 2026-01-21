# PURPOSE:
# Build the exact data list expected by the BFL Stan model.

build_bfl_stan_data <- function(aligned) {

  phi_list <- aligned$aligned_phi
  M <- length(phi_list)
  N <- nrow(phi_list[[1]])
  C <- length(aligned$global_causes)

  phi_array <- array(0, dim = c(N, C, M))
  model_presence <- matrix(0, nrow = C, ncol = M)

  for (m in seq_len(M)) {
    nonzero <- which(colSums(phi_list[[m]]) > 0)
    phi_array[, nonzero, m] <- phi_list[[m]][, nonzero]
    model_presence[nonzero, m] <- 1
  }

  list(
    N = N,
    C_max = C,
    M = M,
    num_causes = C,
    causes = aligned$global_causes,
    model_presence = model_presence,
    count = rowSums(model_presence),
    phi = phi_array
  )
}
