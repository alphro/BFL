#' LCVA Interface Layer
#'
#' This file provides a clean interface around the \pkg{LCVA}
#' package for fitting and predicting on single or multiple
#' datasets. These wrappers standardize cause handling,
#' posterior summaries, and prediction outputs for downstream use.
#'
#' @name LCVA_interface_layer
#' @keywords internal
NULL

# ============================================================
# 1) Fit LCVA
# ============================================================

#' Fit LCVA model on a single training dataset
#'
#' Trains an LCVA model on labeled symptom data and returns
#' a reusable model object.
#'
#' @param X_train Numeric matrix/data.frame (N x P).
#' @param Y_train Vector of training causes.
#' @param lcva_args Optional LCVA hyperparameters (K, Nitr, thin, seed).
#'
#' @return A list with:
#' \describe{
#'   \item{fit}{LCVA fit object from \code{LCVA::LCVA.train()}.}
#'   \item{cause_ids}{Character vector of cause labels in internal order.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rbinom(40 * 8, 1, 0.3), 40, 8)
#' Y <- sample(c("A","B","C"), 40, replace = TRUE)
#' model <- fit_lcva(X, Y, lcva_args = list(Nitr = 200))
#' }
#'
#' @export
fit_lcva <- function(X_train, Y_train, lcva_args = list()) {

  X_train <- as.matrix(X_train)
  storage.mode(X_train) <- "numeric"

  Y_fac     <- factor(Y_train)
  cause_ids <- levels(Y_fac)
  Y_int     <- as.integer(Y_fac)

  K    <- if (!is.null(lcva_args$K))    lcva_args$K    else 5
  Nitr <- if (!is.null(lcva_args$Nitr)) lcva_args$Nitr else 2000
  thin <- if (!is.null(lcva_args$thin)) lcva_args$thin else 2
  seed <- if (!is.null(lcva_args$seed)) lcva_args$seed else 12345

  fit <- LCVA::LCVA.train(
    X = X_train,
    Y = Y_int,
    Domain = rep(1, length(Y_int)),
    K = K,
    model = "S",
    Nitr = Nitr,
    thin = thin,
    seed = seed,
    verbose = FALSE
  )

  list(
    fit = fit,
    cause_ids = cause_ids
  )
}


# ============================================================
# 2) Predict LCVA
# ============================================================

#' Predict LCVA on a target dataset
#'
#' Runs \code{LCVA::LCVA.pred()} on \code{X_target} using a fitted model and returns
#' target-level posterior summaries and top-1 predicted causes.
#'
#' @param lcva_model Output from \code{fit_lcva()}.
#' @param X_target Numeric matrix/data.frame (N x P).
#' @param pred_Nitr Number of MCMC iterations for prediction (default \code{4000}).
#'
#' @return A list with:
#' \describe{
#'   \item{posterior_phi}{Numeric matrix (N x C) of per-record posterior mean cause
#'     probabilities/scores for the target set (rows align to \code{X_target}).}
#'   \item{cause_ids}{Character vector of cause labels (length C).}
#'   \item{pi_pred}{Numeric vector (length C) of estimated target cause fractions.}
#'   \item{Y_pred}{Factor (length N) of top-1 predicted causes.}
#'   \item{target_info}{List with \code{N}, \code{P}, and \code{row_hash}.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X_train <- matrix(rbinom(40 * 8, 1, 0.3), 40, 8)
#' Y_train <- sample(c("A","B","C"), 40, replace = TRUE)
#' X_test  <- matrix(rbinom(20 * 8, 1, 0.3), 20, 8)
#'
#' model <- fit_lcva(X_train, Y_train, lcva_args = list(Nitr = 200))
#' preds <- predict_lcva(model, X_test, pred_Nitr = 400)
#' dim(preds$posterior_phi)  # N x C
#' }
#'
#' @export
predict_lcva <- function(lcva_model, X_target, pred_Nitr = 4000) {

  stopifnot(!is.null(lcva_model$fit),
            !is.null(lcva_model$cause_ids))

  X_target <- as.matrix(X_target)
  storage.mode(X_target) <- "numeric"

  # Hash copy (avoid mutating original)
  X_hash <- X_target
  rownames(X_hash) <- NULL
  row_hash <- hash_rows(X_hash)

  out <- LCVA::LCVA.pred(
    fit = lcva_model$fit,
    X_test = X_target,
    model = "C",
    Nitr = pred_Nitr,
    return_likelihood = TRUE,
    verbose = FALSE
  )

  posterior_phi <- apply(out$x_given_y_prob, c(2,3), mean)
  posterior_phi[!is.finite(posterior_phi)] <- 0
  posterior_phi[posterior_phi == 0] <- .Machine$double.xmin

  pi_pred <- colMeans(out$pi_test)

  Y_pred <- get_assignment(out$Y_test)
  Y_pred_factor <- factor(
    Y_pred,
    levels = seq_along(lcva_model$cause_ids),
    labels = lcva_model$cause_ids
  )

  rm(out); gc()

  list(
    posterior_phi = posterior_phi,
    cause_ids     = lcva_model$cause_ids,
    pi_pred       = pi_pred,
    Y_pred        = Y_pred_factor,
    target_info = list(
      N = nrow(X_target),
      P = ncol(X_target),
      row_hash = row_hash
    )
  )
}


# ============================================================
# 3) Single LCVA (fit + predict)
# ============================================================

#' Fit and predict LCVA on a single target dataset
#'
#' Convenience wrapper equivalent to calling
#' \code{fit_lcva()} followed by \code{predict_lcva()}.
#'
#' @param X_train Training symptom matrix.
#' @param Y_train Training causes.
#' @param X_target Target symptom matrix.
#' @param lcva_args Optional LCVA hyperparameters.
#'
#' @return Same output structure as \code{predict_lcva()}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X_train <- matrix(rbinom(50 * 8, 1, 0.3), 50, 8)
#' Y_train <- sample(c("A","B","C"), 50, replace = TRUE)
#' X_target <- matrix(rbinom(20 * 8, 1, 0.35), 20, 8)
#'
#' result <- single_lcva(
#'   X_train,
#'   Y_train,
#'   X_target,
#'   lcva_args = list(Nitr = 200, pred_Nitr = 400)
#' )
#' str(result)
#' }
#'
#' @export
single_lcva <- function(X_train, Y_train, X_target, lcva_args = list()) {

  model <- fit_lcva(X_train, Y_train, lcva_args)

  pred_Nitr <- if (!is.null(lcva_args$pred_Nitr))
    lcva_args$pred_Nitr else 4000

  predict_lcva(model, X_target, pred_Nitr = pred_Nitr)
}


# ============================================================
# 4) Multi LCVA (fit once, predict many)
# ============================================================

#' Fit LCVA once and predict on multiple target datasets
#'
#' Internally equivalent to:
#' \preformatted{
#'   model <- fit_lcva(X_train, Y_train, lcva_args)
#'   lapply(targets, function(X) predict_lcva(model, X))
#' }
#'
#' @param X_train Training symptom matrix.
#' @param Y_train Training causes.
#' @param targets Named list of target matrices.
#' @param lcva_args Optional LCVA hyperparameters.
#'
#' @return A list with:
#' \describe{
#'   \item{cause_ids}{Cause labels in training model.}
#'   \item{targets}{Named list of prediction results per site.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#'
#' X_train <- matrix(rbinom(60 * 8, 1, 0.3), 60, 8)
#' Y_train <- sample(c("A","B","C"), 60, replace = TRUE)
#'
#' targets <- list(
#'   site1 = matrix(rbinom(20 * 8, 1, 0.30), 20, 8),
#'   site2 = matrix(rbinom(25 * 8, 1, 0.35), 25, 8),
#'   site3 = matrix(rbinom(30 * 8, 1, 0.25), 30, 8),
#'   site4 = matrix(rbinom(18 * 8, 1, 0.40), 18, 8)
#' )
#'
#' results <- multi_lcva(
#'   X_train,
#'   Y_train,
#'   targets,
#'   lcva_args = list(Nitr = 200, pred_Nitr = 400)
#' )
#'
#' names(results$targets)
#' str(results$targets$site1)
#' }
#'
#' @export
multi_lcva <- function(X_train, Y_train, targets, lcva_args = list()) {

  if (!is.list(targets) || is.null(names(targets)))
    stop("targets must be a named list of matrices.")

  model <- fit_lcva(X_train, Y_train, lcva_args)

  pred_Nitr <- if (!is.null(lcva_args$pred_Nitr))
    lcva_args$pred_Nitr else 4000

  results <- lapply(targets, function(X_target) {
    predict_lcva(model, X_target, pred_Nitr)
  })
  names(results) <- names(targets)

  list(
    cause_ids = model$cause_ids,
    targets   = results
  )
}
