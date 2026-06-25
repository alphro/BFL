# ============================================================
# build_bfl_stan_data.R
# ============================================================
# PURPOSE:
#   Build Stan data lists for global BFL models.
#
# DESIGN:
#   - build_bfl_stan_data(): shared core (used by ALL BFL variants)
#   - build_bfl_stan_data_balanced(): adds partial-label fields
#   - build_bfl_stan_data_unbalanced(): adds partial-label fields
#
# IMPORTANT:
#   The difference between balanced vs unbalanced lives in the
#   STAN MODEL, not the data structure.
# ============================================================


# ------------------------------------------------------------
# Core builder: used by ALL Stan models
# ------------------------------------------------------------
# Contributes to:
#   - no-partial-label Stan model
#   - balanced partial-label Stan model
#   - unbalanced / label-shift Stan model
#
# This function MUST remain schema-stable.
# ------------------------------------------------------------

#' Build Stan data (core) for BFL
#'
#' Internal helper that converts an aligned local-summaries object into the
#' schema used by all BFL Stan variants.
#'
#' @param aligned Output of `align_local_summaries()`. Must contain:
#'   - `aligned_phi`: named list of N x C matrices (one per local model).
#'   - `model_presence`: integer matrix (C x M) of structural cause membership
#'     produced by `align_local_summaries()`.
#'
#' @return A Stan data list shared across all BFL variants.
#'
#' @keywords internal
build_bfl_stan_data <- function(aligned) {

  phi_list <- aligned$aligned_phi
  M <- length(phi_list)
  N <- nrow(phi_list[[1]])
  C <- ncol(phi_list[[1]])

  # ----------------------------------------------------------
  # phi array: N x C x M
  # ----------------------------------------------------------
  phi_array <- array(0, dim = c(N, C, M))
  for (m in seq_len(M)) {
    stopifnot(
      nrow(phi_list[[m]]) == N,
      ncol(phi_list[[m]]) == C
    )
    phi_array[, , m] <- phi_list[[m]]
  }

  # ----------------------------------------------------------
  # model_presence: C x M  (structural cause membership)
  #
  # Taken directly from aligned$model_presence, which is built
  # in align_local_summaries() from cause_ids — NOT from phi
  # values.  A source that structurally contains a cause but
  # whose phi has numerically underflowed to all-zeros still
  # gets model_presence=1, so the underflow is visible to the
  # run_bfl_gibbs() diagnostic rather than being silently
  # reclassified as structural absence.
  # ----------------------------------------------------------
  model_presence <- aligned$model_presence
  if (is.null(model_presence)) {
    stop(
      "aligned$model_presence is missing. ",
      "Re-run align_local_summaries() to regenerate it."
    )
  }
  stopifnot(identical(dim(model_presence), c(C, M)))

  list(
    # Dimensions
    N = as.integer(N),
    C_max = as.integer(C),
    M = as.integer(M),
    num_causes = as.integer(C),

    # Cause indexing (canonicalized to 1..C)
    causes = as.integer(seq_len(C)),

    # Model structure (structural, from cause_ids — not phi)
    model_presence = model_presence,
    count = as.integer(rowSums(model_presence)),

    # Likelihood inputs
    phi = phi_array
  )
}


# ------------------------------------------------------------
# Balanced partial-label builder
# ------------------------------------------------------------
# Contributes to:
#   - balanced partial-label Stan model
#
# Adds:
#   - Y_known[N]
#   - Y[N]  (cause index in 1..num_causes)
# ------------------------------------------------------------

#' Build Stan data (balanced partial labels) for BFL
#'
#' Internal helper that augments the core BFL Stan data with partial-label fields
#' for the *balanced* model variant.
#'
#' @param aligned Output of `align_local_summaries()`.
#' @param Y_known Integer vector (length N): 1 if label observed, 0 if missing.
#' @param Y_idx Integer vector (length N): cause index in 1..C (dummy value allowed
#'   where Y_known == 0).
#'
#' @return Stan data list for the balanced partial-label Stan model.
#'
#' @keywords internal
build_bfl_stan_data_balanced <- function(aligned, Y_known, Y_idx) {

  stan_data <- build_bfl_stan_data(aligned)

  stopifnot(length(Y_known) == stan_data$N)
  stopifnot(length(Y_idx)   == stan_data$N)

  stan_data$Y_known <- as.integer(Y_known)
  stan_data$Y_idx   <- as.integer(Y_idx)

  stan_data
}


# ------------------------------------------------------------
# Unbalanced / label-shift partial-label builder
# ------------------------------------------------------------
# Contributes to:
#   - unbalanced (label-shift) Stan model
#
# NOTE:
#   Data structure is IDENTICAL to balanced case.
#   The difference is entirely in the Stan model:
#     - balanced: one prevalence pi
#     - unbalanced: pi + pi_O
# ------------------------------------------------------------

#' Build Stan data (unbalanced partial labels) for BFL
#'
#' Internal helper that augments the core BFL Stan data with partial-label fields
#' for the *unbalanced/label-shift* model variant.
#'
#' @param aligned Output of `align_local_summaries()`.
#' @param Y_known Integer vector (length N): 1 if label observed, 0 if missing.
#' @param Y_idx Integer vector (length N): cause index in 1..C (dummy value allowed
#'   where Y_known == 0).
#'
#' @return Stan data list for the unbalanced partial-label Stan model.
#'
#' @keywords internal
build_bfl_stan_data_unbalanced <- function(aligned, Y_known, Y_idx) {

  stan_data <- build_bfl_stan_data(aligned)

  stopifnot(length(Y_known) == stan_data$N)
  stopifnot(length(Y_idx)   == stan_data$N)

  stan_data$Y_known <- as.integer(Y_known)
  stan_data$Y_idx   <- as.integer(Y_idx)

  stan_data
}
