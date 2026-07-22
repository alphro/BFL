#' Run global BFL aggregation
#'
#' Aggregates per-model per-record class scores into a global Bayesian ensemble.
#' The Stan model is decided purely by \strong{what you pass in} — there
#' is no model-type parameter and no row inference. The rows of
#' \code{local_summaries}, \code{X_target} and \code{Y_target} must line up 1:1;
#' the caller assembles any labeled rows (and the self-model column) beforehand.
#'
#' @section What you pass in decides the model:
#' There is no variant flag. The Stan model and CSMF correction are chosen from
#' three inputs only:
#' \describe{
#'   \item{\code{Y_target}}{\code{NULL} (or all \code{NA}) selects the
#'     no-partial-label model. A vector with labels (\code{NA} on unlabeled
#'     rows, labels on the rest) selects the partial-label model.}
#'   \item{\code{Y_add}}{when supplied, its held-out labels drive the CSMF
#'     correction \eqn{(n_{total}\hat\pi + n_{Lc})/(n_{total}+n_{add})}; when
#'     \code{NULL}, no correction is applied.}
#'   \item{\code{label_shift}}{with partial labels, \code{FALSE} uses the
#'     balanced Stan variant, \code{TRUE} the unbalanced (shift) variant.}
#' }
#'
#' Any self-model or additional local model is just another entry in
#' \code{local_summaries}. The caller assembles every row 1:1 with
#' \code{X_target} (appending labeled rows where relevant) before calling; the
#' package reads only the three inputs above and infers nothing else.
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
#' @param Y_target Optional character/factor vector of length \code{nrow(X_target)}
#'   giving target labels. Use \code{NA} for unlabeled records. If \code{NULL}
#'   (or all \code{NA}), the no-partial-label Stan model is used.
#' @param Y_add Optional character/factor vector of held-out labels for rows that
#'   are \emph{not} in \code{X_target}, used only for the CSMF correction
#'   \eqn{(n_{total}\hat\pi + n_{Lc})/(n_{total}+n_{add})}. Only meaningful when
#'   \code{label_shift = FALSE}. \code{NULL} = no CSMF correction.
#' @param label_shift Logical; if \code{TRUE} and \code{Y_target} is provided,
#'   uses the unbalanced/label-shift Stan variant. Default \code{FALSE}.
#' @param sampler One of \code{"gibbs"} (default) or \code{"stan"}. Both samplers
#'   support all three variants (\code{no_partial}, \code{balanced},
#'   \code{unbalanced}). \code{"gibbs"} uses the conjugate Rcpp Gibbs sampler,
#'   which is much faster; \code{"stan"} runs the matching Stan/NUTS model. The
#'   variant itself is chosen from the inputs (see above), not from this
#'   argument, so partial labels work under either sampler.
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
#'   \item{stan_idx}{Integer indices of rows that entered the model
#'     (always \code{1:nrow(X_target)}).}
#'   \item{n_total}{Number of model rows (\code{nrow(X_target)}).}
#'   \item{n_add}{Number of held rows behind \code{nLc} (0 if none).}
#'   \item{nLc}{Named integer vector of cause counts from \code{Y_add}
#'     (used for CSMF correction; NULL if \code{Y_add} not supplied).}
#'   \item{has_labels}{Logical: whether partial labels were passed to Stan.}
#'   \item{label_shift}{Logical: whether the label-shift variant was used.}
#' }
#'
#' @seealso \code{\link{compute_row_hashes}}, \code{\link{score_BFL}}
#'
#' @export
run_BFL <- function(
    local_summaries,
    X_target,
    Y_target    = NULL,
    Y_add       = NULL,
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
  # 2. Reference row hashes from X_target (used for the row_hash field
  #    and downstream scoring, not to locate Stan rows — rows are 1:1).
  # ------------------------------------------------------------------
  X_target     <- as.matrix(X_target)
  storage.mode(X_target) <- "numeric"
  n_total      <- nrow(X_target)
  ref_hashes   <- compute_row_hashes(X_target)
  stan_hashes  <- local_summaries[[1]]$target_info$row_hash

  # ------------------------------------------------------------------
  # 3. Invariant: local_summaries rows == X_target rows == length(Y_target).
  #    The caller assembles every row (appending any labeled rows, and adding
  #    a self-model column to local_summaries where relevant) BEFORE calling.
  #    Because rows are 1:1 there is no inference — every row enters the model.
  # ------------------------------------------------------------------
  n_phi <- nrow(local_summaries[[1]]$posterior_phi)
  if (n_phi != n_total)
    stop("local_summaries has ", n_phi, " rows but X_target has ", n_total,
         ". These must match 1:1 — assemble any labeled rows into ",
         "X_target (and the local summaries) before calling run_BFL().")
  if (!is.null(Y_target) && length(Y_target) != n_total)
    stop("length(Y_target) (", length(Y_target), ") must equal nrow(X_target) (",
         n_total, ").")

  stan_idx <- seq_len(n_total)

  # ------------------------------------------------------------------
  # 4. CSMF-correction metadata comes ONLY from Y_add now (the held-out
  #    labels for rows that are NOT in X_target). nLc = their cause counts.
  # ------------------------------------------------------------------
  nLc <- if (!is.null(Y_add) && length(Y_add) > 0) {
    keep <- !is.na(Y_add)
    if (any(keep)) table(factor(as.character(Y_add[keep]))) else NULL
  } else NULL

  # ------------------------------------------------------------------
  # 5. Determine Stan variant purely from the inputs:
  #    - Y_target NULL (or all NA)        → no_partial
  #    - Y_target has labels, no shift    → balanced
  #    - Y_target has labels, label_shift → unbalanced
  #    A self-model column, if present, is just another local summary;
  #    run_BFL does not treat it specially.
  # ------------------------------------------------------------------
  has_labels <- !is.null(Y_target) && any(!is.na(Y_target))
  variant    <- if (!has_labels) {
    "no_partial"
  } else if (isTRUE(label_shift)) {
    "unbalanced"
  } else {
    "balanced"
  }

  # ------------------------------------------------------------------
  # 6. Align local summaries to global cause set (rows are positional)
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
      n_total     = n_total,                          # model rows (N, N+L, or N+P)
      n_add       = if (is.null(nLc)) 0L else sum(nLc), # held rows behind nLc
      nLc         = nLc,
      has_labels  = has_labels,
      label_shift = isTRUE(label_shift)
    ),
    class = "BFL"
  )
}


# ------------------------------------------------------------------
# Internal: build Stan data from aligned summaries + Y_target
# ------------------------------------------------------------------
.build_stan_data <- function(aligned, Y_target, stan_idx, variant) {

  if (variant == "no_partial") {
    return(build_bfl_stan_data(aligned))
  }

  # When partial labels are present: subset Y_target to Stan rows
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
