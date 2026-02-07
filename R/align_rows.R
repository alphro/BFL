# PURPOSE:
# Reorder each local summary's posterior_phi rows to match a reference row_hash.
# For now, requires full coverage (a pure permutation). Later, this is where
# partial coverage / masking logic will live.

align_rows_local_summaries <- function(local_summaries, ref_row_hash) {
  stopifnot(is.character(ref_row_hash))
  stopifnot(length(ref_row_hash) >= 1)
  
  lapply(local_summaries, function(x) {
    stopifnot(!is.null(x$target_info))
    stopifnot(!is.null(x$target_info$row_hash))
    
    idx <- match(ref_row_hash, x$target_info$row_hash)
    
    # For now require full permutation (no missing rows)
    stopifnot(!anyNA(idx))
    
    x$posterior_phi <- x$posterior_phi[idx, , drop = FALSE]
    x$target_info$row_hash <- ref_row_hash
    x$target_info$N <- length(ref_row_hash)
    
    x
  })
}