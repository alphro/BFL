#' Summarize BFL site contributions and similarity
#'
#' Provides a diagnostic summary of a global BFL fit:
#' - lambda weights per cause × site
#' - estimated target prevalence per cause (pi)
#' - correlation of per-site posterior likelihoods (phi)
#'
#' @param global_fit Output from \code{run_BFL()}.
#'
#' @return A list with components:
#' \describe{
#'   \item{lambda_table}{Data.frame of posterior mean lambda: cause × site.}
#'   \item{pi_summary}{Numeric vector: mean target prevalence per cause.}
#'   \item{phi_correlation}{Matrix: pairwise correlations between sites' phi flattened across observations and causes.}
#' }
#' @export
eval_BFL <- function(global_fit) {
  
  stopifnot(!is.null(global_fit$lambda),
            !is.null(global_fit$phi),
            !is.null(global_fit$causes),
            !is.null(global_fit$model_names))
  
  # --------------------------
  # Lambda table (cause × site)
  # --------------------------
  S <- dim(global_fit$lambda)[1]  # posterior samples
  C <- dim(global_fit$lambda)[2]  # causes
  M <- dim(global_fit$lambda)[3]  # sites/models
  
  lambda_mean <- apply(global_fit$lambda, c(2,3), mean)  # C × M
  rownames(lambda_mean) <- global_fit$causes
  colnames(lambda_mean) <- global_fit$model_names
  
  lambda_table <- as.data.frame(lambda_mean)
  lambda_table$cause <- rownames(lambda_table)
  lambda_table <- lambda_table[, c("cause", global_fit$model_names)]
  
  # --------------------------
  # pi summary
  # --------------------------
  pi_summary <- colMeans(global_fit$pi)
  names(pi_summary) <- global_fit$causes
  
  # --------------------------
  # Phi correlation across sites
  # --------------------------
  phi_mat <- global_fit$phi   # N × C × M
  phi_flat <- sapply(seq_len(M), function(m) {
    as.vector(phi_mat[,,m])
  })
  colnames(phi_flat) <- global_fit$model_names
  phi_correlation <- cor(phi_flat)
  
  list(
    lambda_table = lambda_table,
    pi_summary = pi_summary,
    phi_correlation = phi_correlation
  )
}