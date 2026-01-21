#' Prepare local LCVA summaries for BFL
#'
#' PURPOSE:
#' Runs LCVA on *local individual-level data* and produces a
#' shareable summary that can be sent to the global BFL aggregator.
#'
#' This function is run ON-SITE.
#' Raw data never needs to leave the local machine.
#'
#' INPUT:
#' - X: feature matrix (N x P)
#' - Y: optional cause labels (length N), if supervised LCVA
#'
#' OUTPUT:
#' A local summary object containing ONLY model summaries,
#' not individual-level data.
#'
#' @export
prep_BFL_local_LCVA <- function(
    X,
    Y = NULL,
    K = 5,
    lcva_args = list()
) {
  stop("prep_BFL_local_LCVA not implemented yet")
}
