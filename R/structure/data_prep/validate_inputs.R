validate_local_summaries <- function(local_summaries) {

  stopifnot(is.list(local_summaries))
  stopifnot(length(local_summaries) >= 1)

  for (x in local_summaries) {
    stopifnot(is.matrix(x$posterior_phi))
    stopifnot(is.numeric(x$cause_ids))
    stopifnot(ncol(x$posterior_phi) == length(x$cause_ids))
  }
}

validate_global_fit <- function(global_fit) {
  stopifnot(is.list(global_fit))
  stopifnot(!is.null(global_fit$pi))
  stopifnot(!is.null(global_fit$phi))
}
