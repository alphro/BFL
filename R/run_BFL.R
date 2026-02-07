#' Run global BFL aggregation
#'
#' Aggregates local LCVA summaries into a global BFL model using Stan.
#' This function never sees individual-level data from *other* sites, but
#' it can take target-site X_target (for row hashing / alignment) and optional
#' target labels for partial-label variants.
#'
#' @param local_summaries Named list of local summary objects
#' @param X_target Target symptom matrix/data.frame used only to define canonical row order via hashing
#' @param Y_target Optional target labels (use NA for unknown). If NULL, runs no-partial-label model.
#' @param label_shift Logical; if TRUE and Y_target is provided, uses unbalanced/label-shift variant.
#' @param stan_args List of arguments passed to rstan::sampling
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
    # no-partial-label (current behavior)
    stan_data <- build_bfl_stan_data(aligned)

  } else {
    # partial labels: build Y_known + Y_idx (index into aligned$global_causes)
    y <- as.character(Y_target)
    stopifnot(length(y) == nrow(X_target))

    Y_known <- as.integer(!is.na(y))

    cause_to_idx <- setNames(seq_along(aligned$global_causes), aligned$global_causes)
    Y_idx <- unname(cause_to_idx[y])

    # Stan requires integers everywhere; fill unknowns with 1 (ignored when Y_known==0)
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
    row_hash = ref_row_hash,
    stan_fit = fit$stan_fit
  )
}
