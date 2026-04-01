#' LCVA Interface Layer
#'
#' Clean interface around \pkg{LCVA} for fitting and predicting
#' on single or multiple datasets.
#'
#' @name LCVA_interface_layer
#' @keywords internal
NULL


# ============================================================
# Internal helper
# ============================================================

.get_pred_nitr <- function(lcva_args) {
  if (!is.null(lcva_args$pred_Nitr)) lcva_args$pred_Nitr else 4000
}


# ============================================================
# 1) Fit LCVA
# ============================================================

#' Fit LCVA model on a single training dataset
#'
#' @param X_train Numeric matrix/data.frame (N x P).
#' @param Y_train Vector of training causes.
#' @param lcva_args Optional LCVA hyperparameters (K, Nitr, thin, seed).
#'
#' @return A list with \code{fit} and \code{cause_ids}.
#'
#' @examples
#' \dontrun{
#' model <- fit_lcva(X_train, Y_train, lcva_args = list(Nitr = 200))
#' }
#'
#' @export
fit_lcva <- function(X_train, Y_train, lcva_args = list()) {

  X_train <- as.matrix(X_train)
  storage.mode(X_train) <- "numeric"

  Y_fac     <- factor(Y_train)
  cause_ids <- levels(Y_fac)
  Y_int     <- as.integer(Y_fac)

  fit <- LCVA::LCVA.train(
    X       = X_train,
    Y       = Y_int,
    Domain  = rep(1, length(Y_int)),
    K       = if (!is.null(lcva_args$K))    lcva_args$K    else 5,
    model   = "S",
    Nitr    = if (!is.null(lcva_args$Nitr)) lcva_args$Nitr else 2000,
    thin    = if (!is.null(lcva_args$thin)) lcva_args$thin else 2,
    seed    = if (!is.null(lcva_args$seed)) lcva_args$seed else 12345,
    verbose = FALSE
  )

  list(fit = fit, cause_ids = cause_ids)
}


# ============================================================
# 2) Predict LCVA
# ============================================================

#' Predict LCVA on a target dataset
#'
#' @param lcva_model Output from \code{fit_lcva()}.
#' @param X_target Numeric matrix/data.frame (N x P).
#' @param pred_Nitr MCMC iterations for prediction (default \code{4000}).
#'
#' @return A list with \code{posterior_phi}, \code{cause_ids},
#'   \code{pi_pred}, \code{Y_pred}, and \code{target_info}.
#'   \code{target_info} contains \code{row_hash} (per-row hashes),
#'   \code{N}, \code{P}, and \code{dataset_hash} (audit fingerprint).
#'
#' @examples
#' \dontrun{
#' preds <- predict_lcva(model, X_test, pred_Nitr = 400)
#' }
#'
#' @export
predict_lcva <- function(lcva_model, X_target, pred_Nitr = 4000) {

  stopifnot(
    !is.null(lcva_model$fit),
    !is.null(lcva_model$cause_ids)
  )

  X_target <- as.matrix(X_target)
  storage.mode(X_target) <- "numeric"

  out <- LCVA::LCVA.pred(
    fit               = lcva_model$fit,
    X_test            = X_target,
    model             = "C",
    Nitr              = pred_Nitr,
    return_likelihood = TRUE,
    verbose           = FALSE
  )

  posterior_phi <- apply(out$x_given_y_prob, c(2, 3), mean)
  posterior_phi[!is.finite(posterior_phi)] <- 0
  posterior_phi[posterior_phi == 0]        <- .Machine$double.xmin

  pi_pred <- colMeans(out$pi_test)

  Y_pred_factor <- factor(
    get_assignment(out$Y_test),
    levels = seq_along(lcva_model$cause_ids),
    labels = lcva_model$cause_ids
  )

  rm(out); gc()

  # Per-row hashes for alignment in run_BFL()
  row_hash <- compute_row_hashes(X_target)

  # Single whole-matrix hash for audit purposes
  X_audit <- X_target
  rownames(X_audit) <- NULL
  colnames(X_audit) <- NULL
  dataset_hash <- rlang::hash(X_audit)

  list(
    posterior_phi = posterior_phi,
    cause_ids     = lcva_model$cause_ids,
    pi_pred       = pi_pred,
    Y_pred        = Y_pred_factor,
    target_info   = list(
      N            = nrow(X_target),
      P            = ncol(X_target),
      row_hash     = row_hash,     # per-row: used for alignment + CSMF correction
      dataset_hash = dataset_hash  # whole-matrix: audit trail only
    )
  )
}


# ============================================================
# 3) Run LCVA (unified entry point)
# ============================================================

#' Fit LCVA and predict on a target dataset
#'
#' Unified entry point for one or many training sites predicting on a single
#' shared target matrix. Handles both the single-site and multi-site cases
#' with the same call.
#'
#' @section Single-site usage:
#' Pass \code{X_train} as a matrix and \code{Y_train} as a vector. Returns
#' the \code{predict_lcva()} output directly (not wrapped in a list).
#'
#' @section Multi-site usage:
#' Pass \code{X_train} as a \emph{named} list of training matrices and
#' \code{Y_train} as a named list of cause vectors (same names, same order).
#' Returns a named list of \code{predict_lcva()} outputs, one per site —
#' ready to pass directly to \code{run_BFL()} as \code{local_summaries}.
#'
#' @section Building local_summaries for each BFL variant:
#' \describe{
#'   \item{Base}{Pass all N rows of \code{X_target}. Call
#'     \code{run_BFL(local_summaries, X_target)}.}
#'   \item{Domain}{Pass \strong{only the unlabeled rows}
#'     (\code{X_target[unlabeled_idx, ]}) as \code{X_target} here. Also
#'     add a self-site trained on the labeled rows, also predicting on
#'     \code{X_unlabeled}. Then call
#'     \code{run_BFL(local_summaries, X_full, Y_target = Y_partial)}
#'     where \code{X_full} is all N rows and \code{Y_partial} has the
#'     known labels (NA for unlabeled). Stan uses the no-partial model
#'     because labeled rows are outside \code{stan_idx}; the labels feed
#'     the automatic CSMF correction only.}
#'   \item{Partial}{Pass all N rows of \code{X_target}. Call
#'     \code{run_BFL(local_summaries, X_target, Y_target = Y_partial)}.}
#'   \item{Mix}{Pass all N rows for source sites; also add a self-site
#'     trained on labeled rows predicting on all N rows. Call
#'     \code{run_BFL(local_summaries, X_target, Y_target = Y_partial)}.}
#' }
#'
#' @param X_train Training symptom matrix (N_train \eqn{\times} P) \emph{or}
#'   a named list of training matrices (one per site).
#' @param Y_train Training cause vector (length N_train) \emph{or} a named
#'   list of cause vectors matching \code{X_train}.
#' @param X_target Target symptom matrix (N \eqn{\times} P). Always pass all
#'   N rows; \code{run_BFL()} infers the correct subset via row hashing.
#' @param lcva_args Optional named list of LCVA hyperparameters passed to
#'   \code{fit_lcva()} and \code{predict_lcva()} (K, Nitr, thin, seed,
#'   pred_Nitr).
#'
#' @return
#' \itemize{
#'   \item \strong{Single site}: the \code{predict_lcva()} output directly
#'     (\code{posterior_phi}, \code{cause_ids}, \code{pi_pred}, \code{Y_pred},
#'     \code{target_info}).
#'   \item \strong{Multi-site}: a named list of \code{predict_lcva()} outputs,
#'     one per site.
#' }
#'
#' @examples
#' \dontrun{
#' # --- Single site ---
#' s1 <- run_lcva(X1, Y1, X_target)
#'
#' # --- Multiple sites in one call ---
#' local_summaries <- run_lcva(
#'   X_train  = list(site1 = X1, site2 = X2, site3 = X3),
#'   Y_train  = list(site1 = Y1, site2 = Y2, site3 = Y3),
#'   X_target = X_target
#' )
#'
#' # --- Domain: add a site trained on labeled target data ---
#' local_summaries_domain <- c(
#'   local_summaries,
#'   list(target_lbl = run_lcva(X_lbl, Y_lbl, X_target))
#' )
#'
#' # --- Pass directly to run_BFL ---
#' fit_base <- run_BFL(local_summaries, X_target)
#' fit_dom  <- run_BFL(local_summaries_domain, X_target)
#' fit_par  <- run_BFL(local_summaries, X_target, Y_target = Y_target)
#' fit_mix  <- run_BFL(local_summaries_domain, X_target, Y_target = Y_target)
#' }
#'
#' @export
run_lcva <- function(X_train, Y_train, X_target, lcva_args = list()) {

  multi_site <- is.list(X_train) && !is.data.frame(X_train)

  if (multi_site) {
    # ------------------------------------------------------------------
    # Multi-site: X_train and Y_train are both named lists
    # ------------------------------------------------------------------
    if (!is.list(Y_train))
      stop("When X_train is a list, Y_train must also be a named list.")

    nms <- names(X_train)
    if (is.null(nms) || any(nms == ""))
      stop("X_train list must be fully named (one name per training site).")
    if (!identical(sort(names(Y_train)), sort(nms)))
      stop("names(X_train) and names(Y_train) must match.")

    pred_Nitr <- .get_pred_nitr(lcva_args)

    lapply(nms, function(nm) {
      model <- fit_lcva(X_train[[nm]], Y_train[[nm]], lcva_args)
      predict_lcva(lcva_model = model, X_target = X_target, pred_Nitr = pred_Nitr)
    }) |> setNames(nms)

  } else {
    # ------------------------------------------------------------------
    # Single site: X_train is a matrix, Y_train is a vector
    # ------------------------------------------------------------------
    model     <- fit_lcva(X_train, Y_train, lcva_args)
    pred_Nitr <- .get_pred_nitr(lcva_args)
    predict_lcva(lcva_model = model, X_target = X_target, pred_Nitr = pred_Nitr)
  }
}
