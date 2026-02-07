#' Hash rows of a matrix or data.frame
#'
#' Internal helper for row-wise alignment of target datasets.
#'
#' @param X Matrix or data.frame
#' @param algo Hash algorithm
#'
#' @keywords internal
hash_rows <- function(X, algo = "xxhash64") {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for hash_rows().")
  }
  
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  rownames(X) <- NULL
  
  row_strings <- apply(
    X,
    1,
    function(r) paste(format(r, scientific = FALSE, trim = TRUE), collapse = "|")
  )
  
  vapply(row_strings, digest::digest, FUN.VALUE = character(1), algo = algo)
}