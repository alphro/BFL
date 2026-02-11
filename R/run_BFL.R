#' Run global BFL aggregation
#'
#' Aggregates per-model likelihood summaries into a global Bayesian
#' ensemble using Stan. Each model contributes a matrix of
#' per-observation, per-class likelihoods, and BFL learns
#' posterior prevalence (pi) and model weights (lambda).
#'
#' @param local_summaries Named list of model summary objects.
#'   Each element must contain per-observation likelihoods aligned
#'   to a common target dataset.
#' @param X_target Target feature matrix/data.frame. Used only to
#'   define a canonical row order via hashing.
#' @param Y_target Optional target labels (use NA for unknown).
#'   If NULL, runs the no-partial-label model.
#' @param label_shift Logical; if TRUE and Y_target is provided,
#'   uses the unbalanced/label-shift variant.
#' @param stan_args List of arguments passed to rstan::sampling.
#'
#' @return A list with components:
#' \describe{
#'   \item{pi}{Posterior draws of class prevalence (S x C).}
#'   \item{lambda}{Posterior draws of per-class model weights (S x C x M).}
#'   \item{phi}{Per-model likelihoods aligned to target rows and global classes (N x C x M).}
#'   \item{causes}{Character vector of global class labels (length C).}
#'   \item{model_names}{Character vector of model names (length M).}
#'   \item{row_hash}{Row hashes for X_target (length N).}
#'   \item{stan_fit}{Underlying rstan fit object.}
#' }
#'
#' @details
#' This function is model-agnostic. It does not assume the local
#' summaries were produced by LCVA. Any method that supplies
#' per-observation class likelihoods in the required structure
#' can be used as input.
#'
#' @export
run_BFL <- function(
    local_summaries,
    X_target,
    Y_target = NULL,
    label_shift = FALSE,
    stan_args = list(iter = 2000, chains = 4, seed = 12345)
) {

  validate_local_summaries(local_summaries)

  # ------------------------------------------------------------
  # Canonical target row hashes from X_target
  # ------------------------------------------------------------
  X_target <- as.matrix(X_target)
  storage.mode(X_target) <- "numeric"
  ref_row_hash <- hash_rows(X_target)

  # ------------------------------------------------------------
  # Align rows (to X_target) + align causes (global cause set)
  # ------------------------------------------------------------
  aligned <- align_local_summaries(local_summaries, ref_row_hash = ref_row_hash)

  # ------------------------------------------------------------
  # Build Stan data depending on label availability
  # ------------------------------------------------------------
  if (is.null(Y_target)) {
    stan_data <- build_bfl_stan_data(aligned)
  } else {
    y <- as.character(Y_target)
    stopifnot(length(y) == nrow(X_target))

    Y_known <- as.integer(!is.na(y))

    cause_to_idx <- setNames(seq_along(aligned$global_causes), aligned$global_causes)
    Y_idx <- unname(cause_to_idx[y])

    Y_idx[which(Y_known == 0)] <- 1L

    if (any(Y_known == 1 & is.na(Y_idx))) {
      bad <- unique(y[Y_known == 1 & is.na(Y_idx)])
      stop("Y_target contains labels not in aligned$global_causes: ",
           paste(bad, collapse = ", "))
    }

    if (label_shift) {
      stan_data <- build_bfl_stan_data_unbalanced(aligned, Y_known = Y_known, Y_idx = as.integer(Y_idx))
    } else {
      stan_data <- build_bfl_stan_data_balanced(aligned, Y_known = Y_known, Y_idx = as.integer(Y_idx))
    }
  }

  # ------------------------------------------------------------
  # Decide which Stan variant to run
  # ------------------------------------------------------------
  if (is.null(Y_target)) {

    # No labels provided
    variant <- "no_partial"

  } else {

    if (label_shift) {
      variant <- "unbalanced"
    } else {
      variant <- "balanced"
    }

  }

  fit <- run_bfl_stan(stan_data, stan_args, variant = variant)

  # ------------------------------------------------------------
  # Return per-model phi (N x C x M), NOT collapsed phi_global
  # ------------------------------------------------------------
  phi_list_na <- aligned$aligned_phi
  model_names <- names(phi_list_na)

  # Replace NA with 0 (matches paper-style "missing cause in model => 0 contribution")
  phi_list <- lapply(phi_list_na, function(mat) {
    mat <- as.matrix(mat)
    mat[is.na(mat)] <- 0
    mat
  })

  N <- nrow(phi_list[[1]])
  C <- ncol(phi_list[[1]])
  M <- length(phi_list)

  phi_arr <- array(0, dim = c(N, C, M))
  for (m in seq_len(M)) {
    phi_arr[, , m] <- phi_list[[m]]
  }

  list(
    pi = fit$pi,                 # expect S x C (posterior draws)
    lambda = fit$lambda,         # expect S x C x M (posterior draws)
    phi = phi_arr,               # N x C x M  (phi per model)
    causes = aligned$global_causes,
    model_names = model_names,
    row_hash = ref_row_hash,
    stan_fit = fit$stan_fit
  )
}
