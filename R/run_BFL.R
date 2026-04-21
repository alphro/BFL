#' Run global BFL aggregation
#'
#' Aggregates per-model per-record class scores into a global Bayesian ensemble
#' using Stan. Automatically infers the correct Stan model and CSMF correction
#' from the structure of \code{local_summaries} and \code{Y_target} — no model
#' type parameter required.
#'
#' @section How BFL variants map to inputs:
#' \describe{
#'   \item{Base}{Build \code{local_summaries} predicting on all N rows of
#'     \code{X_target}. Pass \code{X_target} (full N rows) to
#'     \code{run_BFL()}. Leave \code{Y_target = NULL}.}
#'   \item{Domain}{Build \code{local_summaries} predicting on
#'     \code{X_unlabeled} only (the rows whose causes are unknown).
#'     Also add a self-site trained on the labeled rows, predicting on
#'     \code{X_unlabeled}. Pass \code{X_target = X_full} (all N rows) and
#'     \code{Y_target} with the known labels (NA for unlabeled rows). Stan
#'     uses the no-partial model because all labeled rows fall outside
#'     \code{stan_idx}; the labels are used only for the automatic CSMF
#'     correction \eqn{(n_{stan} \cdot \hat\pi + n_{Lc}) / N}.}
#'   \item{Partial}{Build \code{local_summaries} predicting on all N rows.
#'     Pass \code{X_target} (full N rows) and \code{Y_target} with \code{NA}
#'     for unlabeled records. Stan uses the partial-label model.}
#'   \item{Mix}{Same as Partial, but also include an extra self-site in
#'     \code{local_summaries} trained on the labeled target rows, predicting
#'     on all N rows. Stan uses the partial-label model.}
#' }
#'
#' @param local_summaries Named list of local model summary objects. Each
#'   element must contain:
#'   \itemize{
#'     \item \code{posterior_phi}: numeric matrix (N x C).
#'     \item \code{cause_ids}: character vector of length C.
#'     \item \code{target_info$row_hash}: character vector of length N
#'       (produced by \code{compute_row_hashes()}).
#'   }
#' @param X_target Full target feature matrix/data.frame (all N rows). Used
#'   to compute reference row hashes and infer CSMF correction metadata.
#' @param Y_target Optional character/factor vector of length N giving target
#'   labels. Use \code{NA} for unlabeled records. If \code{NULL}, the
#'   no-partial-label Stan model is used.
#' @param label_shift Logical; if \code{TRUE} and \code{Y_target} is provided,
#'   uses the unbalanced/label-shift Stan variant. Default \code{FALSE}.
#' @param sampler One of \code{"gibbs"} (default) or \code{"stan"}.
#'   \code{"gibbs"} uses the conjugate Rcpp Gibbs sampler, which is much faster
#'   and applies when no partial labels enter the sampler (\code{no_partial}
#'   variant).  \code{"stan"} forces the original Stan/NUTS path.  For
#'   partial-label variants (\code{balanced}, \code{unbalanced}) Stan is always
#'   used regardless of this argument.
#' @param mcmc_args List of MCMC controls shared by both samplers: \code{iter}
#'   (default 2000), \code{chains} (default 4), \code{seed} (default 12345),
#'   and \code{init} (Stan only, default \code{"random"}).
#' @param gibbs_args List of Gibbs-specific tuning options, ignored when
#'   \code{sampler = "stan"}: \code{logistic_normal} (logical, default
#'   \code{FALSE} — use conjugate Dirichlet prior; set \code{TRUE} to match
#'   Stan's logistic-normal prior via an MH step) and \code{mh_scale}
#'   (random-walk step size, default \code{0.5}).
#'
#' @return An object of class \code{"BFL"} — a list with:
#' \describe{
#'   \item{pi}{Posterior draws of class prevalence (S x C).}
#'   \item{lambda}{Posterior draws of per-class model weights (S x C x M).}
#'   \item{phi}{Per-model scores aligned to global classes (N x C x M).}
#'   \item{causes}{Character vector of global class labels (length C).}
#'   \item{model_names}{Character vector of model/site names (length M).}
#'   \item{row_hash}{Reference row hashes for X_target (length N).}
#'   \item{stan_idx}{Integer indices of rows that entered Stan.}
#'   \item{n_total}{Total number of target records (N).}
#'   \item{nLc}{Named integer vector of labeled-row counts outside Stan
#'     (used for CSMF correction; NULL if Stan saw all rows).}
#'   \item{has_labels}{Logical: whether partial labels were passed to Stan.}
#'   \item{label_shift}{Logical: whether the label-shift variant was used.}
#' }
#'
#' @seealso \code{\link{compute_row_hashes}}, \code{\link{run_lcva}},
#'   \code{\link{score_BFL}}
#'
#' @export
run_BFL <- function(
    local_summaries,
    X_target,
    Y_target    = NULL,
    label_shift = FALSE,
    sampler     = c("gibbs", "stan"),
    mcmc_args   = list(iter = 2000, chains = 4, seed = 12345),
    gibbs_args  = list()
) {
  sampler <- match.arg(sampler)

  # ------------------------------------------------------------------
  # 1. Validate local_summaries
  # ------------------------------------------------------------------
  validate_local_summaries(local_summaries)

  # ------------------------------------------------------------------
  # 2. Reference row hashes from the full X_target (all N rows)
  # ------------------------------------------------------------------
  X_target     <- as.matrix(X_target)
  storage.mode(X_target) <- "numeric"
  n_total      <- nrow(X_target)
  ref_hashes   <- compute_row_hashes(X_target)

  # ------------------------------------------------------------------
  # 3. Determine which rows entered Stan (stan_idx)
  #
  #    Two cases:
  #
  #    (a) local_summaries cover all N rows of X_target (Base / Partial /
  #        Mix): stan_hashes == ref_hashes, so stan_idx = 1..N via the
  #        identical fast-path in .find_stan_idx.
  #
  #    (b) local_summaries cover only the unlabeled subset of X_target
  #        (Domain): Y_target has NA exactly at those rows. Using hash
  #        inference here is fragile — with binary symptom data the same
  #        168-feature profile can appear at both labeled and unlabeled
  #        positions, so the consume-first matcher may return labeled-row
  #        positions instead of the true unlabeled ones. We therefore
  #        derive stan_idx directly from Y_target in this case.
  # ------------------------------------------------------------------
  stan_hashes    <- local_summaries[[1]]$target_info$row_hash
  n_stan_expect  <- length(stan_hashes)

  if (n_stan_expect < n_total && !is.null(Y_target)) {
    # Domain path: Stan rows == NA rows in Y_target
    stan_idx <- sort(which(is.na(Y_target)))
    if (length(stan_idx) != n_stan_expect)
      stop("Number of NA rows in Y_target (", length(stan_idx), ") does not match ",
           "the number of rows in local_summaries (", n_stan_expect, "). ",
           "For Domain, local_summaries must be built on exactly the unlabeled rows.")
  } else {
    # Base / Partial / Mix: hash matching (hits the identical fast-path)
    stan_idx <- .find_stan_idx(ref_hashes, stan_hashes)
  }
  n_stan <- length(stan_idx)

  # ------------------------------------------------------------------
  # 4. Infer CSMF correction metadata
  #    If Stan saw fewer rows than N, labeled rows outside Stan give nLc
  # ------------------------------------------------------------------
  needs_correction <- n_stan < n_total
  nLc              <- NULL

  if (needs_correction && !is.null(Y_target)) {
    outside_idx  <- setdiff(seq_len(n_total), stan_idx)
    outside_Y    <- Y_target[outside_idx]
    labeled_outside <- outside_Y[!is.na(outside_Y)]
    if (length(labeled_outside) > 0) {
      nLc <- table(factor(as.character(labeled_outside)))
    }
  }

  # ------------------------------------------------------------------
  # 5. Determine Stan variant
  #    Check for labels only within the rows Stan actually sees (stan_idx).
  #    This means Domain can pass Y_target for CSMF correction purposes
  #    while still using the no-partial Stan model — because the labeled
  #    rows sit outside stan_idx and Stan never sees their labels.
  #    - no labels in Stan rows       → no_partial
  #    - labels in Stan rows + label_shift=FALSE  → balanced
  #    - labels in Stan rows + label_shift=TRUE   → unbalanced
  # ------------------------------------------------------------------
  has_labels <- if (is.null(Y_target)) FALSE
                else any(!is.na(Y_target[stan_idx]))
  variant    <- if (!has_labels) {
    "no_partial"
  } else if (isTRUE(label_shift)) {
    "unbalanced"
  } else {
    "balanced"
  }

  # ------------------------------------------------------------------
  # 6. Align local summaries to global cause set
  # ------------------------------------------------------------------
  aligned <- align_local_summaries(local_summaries, ref_row_hash = stan_hashes)

  # ------------------------------------------------------------------
  # 7. Build Stan data
  # ------------------------------------------------------------------
  stan_data <- .build_stan_data(aligned, Y_target, stan_idx, variant)

  # ------------------------------------------------------------------
  # 8. Run sampler
  # ------------------------------------------------------------------
  fit <- run_bfl_stan(stan_data, mcmc_args, gibbs_args,
                      variant = variant, sampler = sampler)

  # ------------------------------------------------------------------
  # 9. Build phi array (N x C x M)
  # ------------------------------------------------------------------
  phi_list    <- lapply(aligned$aligned_phi, as.matrix)
  model_names <- names(phi_list)
  N  <- nrow(phi_list[[1]])
  C  <- ncol(phi_list[[1]])
  M  <- length(phi_list)

  phi_arr <- array(0, dim = c(N, C, M))
  for (m in seq_len(M)) phi_arr[, , m] <- phi_list[[m]]

  # ------------------------------------------------------------------
  # 10. Return classed BFL object
  # ------------------------------------------------------------------
  structure(
    list(
      pi          = fit$pi,
      lambda      = fit$lambda,
      phi         = phi_arr,
      causes      = aligned$global_causes,
      model_names = model_names,
      row_hash    = ref_hashes,
      stan_idx    = stan_idx,
      n_total     = n_total,
      nLc         = nLc,
      has_labels  = has_labels,
      label_shift = isTRUE(label_shift)
    ),
    class = "BFL"
  )
}


# ------------------------------------------------------------------
# Internal: consume-first row-hash matching
#
# Finds exactly length(stan_hashes) positions in ref_hashes by
# consuming each match at most once. This correctly handles the case
# where the same hash appears at both labeled and unlabeled positions
# in X_target (common with binary symptom data).
#
# The identical-vectors fast path handles Base / Partial / Mix where
# local_summaries were built on all N rows of X_target.
# ------------------------------------------------------------------
.find_stan_idx <- function(ref_hashes, stan_hashes) {

  # Fast path: local_summaries cover all N rows in the same order
  if (identical(ref_hashes, stan_hashes)) return(seq_along(ref_hashes))

  # Build a lookup: hash → ordered list of positions in ref_hashes
  hash_positions <- tapply(
    seq_along(ref_hashes), ref_hashes,
    function(x) x, simplify = FALSE
  )

  consumed <- list()   # hash → how many positions already consumed

  idx <- vapply(stan_hashes, function(h) {
    positions <- hash_positions[[h]]
    if (is.null(positions))
      stop("A row in local_summaries has no matching row in X_target. ",
           "Make sure X_target contains all rows that local_summaries ",
           "were built on.")
    n_used <- if (is.null(consumed[[h]])) 0L else consumed[[h]]
    if (n_used >= length(positions))
      stop("More rows in local_summaries match hash '", h, "' than ",
           "exist in X_target. Check that X_target and local_summaries ",
           "were built from the same data.")
    consumed[[h]] <<- n_used + 1L
    positions[[n_used + 1L]]
  }, integer(1))

  sort(idx)
}


# ------------------------------------------------------------------
# Internal: build Stan data from aligned summaries + Y_target
# ------------------------------------------------------------------
.build_stan_data <- function(aligned, Y_target, stan_idx, variant) {

  if (variant == "no_partial") {
    return(build_bfl_stan_data(aligned))
  }

  # Partial-label variants: subset Y_target to Stan rows
  y        <- as.character(Y_target[stan_idx])
  Y_known  <- as.integer(!is.na(y))

  cause_to_idx <- setNames(seq_along(aligned$global_causes), aligned$global_causes)
  Y_idx        <- unname(cause_to_idx[y])
  Y_idx[Y_known == 0] <- 1L

  if (any(Y_known == 1 & is.na(Y_idx))) {
    bad <- unique(y[Y_known == 1 & is.na(Y_idx)])
    stop("Y_target contains labels not in aligned$global_causes: ",
         paste(bad, collapse = ", "))
  }

  if (variant == "unbalanced") {
    build_bfl_stan_data_unbalanced(aligned,
                                   Y_known = Y_known,
                                   Y_idx   = as.integer(Y_idx))
  } else {
    build_bfl_stan_data_balanced(aligned,
                                 Y_known = Y_known,
                                 Y_idx   = as.integer(Y_idx))
  }
}
