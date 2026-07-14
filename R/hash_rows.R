#' Compute row hashes for a target dataset
#'
#' Produces a character vector of per-row hashes using xxHash via
#' \code{rlang::hash()}. This is the published contract between BFL and
#' any local model: attach the output to \code{target_info$row_hash} so
#' \code{run_BFL()} can align rows and compute CSMF corrections automatically.
#'
#' @param X Numeric matrix or data.frame (N x P).
#'
#' @return Character vector of length N — one hash per row.
#'
#' @details
#' All metadata (rownames, colnames) is stripped before hashing so that
#' identical values always produce identical hashes regardless of how the
#' matrix was named. Storage mode is coerced to \code{"numeric"} so that
#' integer and double matrices with the same values hash identically.
#'
#' For non-LCVA users building local summaries manually, call this on your
#' \code{X_target} and attach the result:
#' \preformatted{
#'   my_summary$target_info$row_hash <- compute_row_hashes(X_target)
#' }
#'
#' @examples
#' \dontrun{
#' X <- matrix(rbinom(100, 1, 0.3), nrow = 20)
#' hashes <- compute_row_hashes(X)
#' length(hashes)  # 20
#'
#' # Identical rows always produce the same hash
#' X2 <- X[sample(nrow(X)), ]
#' all(sort(compute_row_hashes(X)) == sort(compute_row_hashes(X2)))  # TRUE
#' }
#'
#' @export
compute_row_hashes <- function(X) {

  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  rownames(X)     <- NULL
  colnames(X)     <- NULL

  apply(X, 1, rlang::hash)
}
