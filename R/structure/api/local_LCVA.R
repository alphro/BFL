local_LCVA <- function(
    X_train,
    Y_train,
    X_target,
    lcva_args = list()
) {
  stopifnot(!is.null(X_train), !is.null(Y_train), !is.null(X_target))

  X_train  <- as.matrix(X_train);  storage.mode(X_train)  <- "numeric"
  X_target <- as.matrix(X_target); storage.mode(X_target) <- "numeric"
  stopifnot(ncol(X_train) == ncol(X_target))

  # Canonicalize target for hashing (ignore rownames)
  X_hash <- X_target
  rownames(X_hash) <- NULL
  target_hash <- digest::digest(X_hash)

  # IMPORTANT: keep cause identities as NAMES (not 1..C)
  Y_fac <- factor(Y_train)
  cause_ids <- levels(Y_fac)
  Y_int <- as.integer(Y_fac)

  K       <- if (!is.null(lcva_args$K)) lcva_args$K else 5
  Nitr    <- if (!is.null(lcva_args$Nitr)) lcva_args$Nitr else 2000
  thin    <- if (!is.null(lcva_args$thin)) lcva_args$thin else 2
  seed    <- if (!is.null(lcva_args$seed)) lcva_args$seed else 12345
  burn_in <- if (!is.null(lcva_args$burn_in)) lcva_args$burn_in else floor(Nitr / 2)

  fit <- LCVA::LCVA.train(
    X = X_train, Y = Y_int,
    Domain = rep(1, length(Y_int)),
    K = K, model = "S",
    Nitr = Nitr, thin = thin, seed = seed,
    verbose = FALSE
  )

  out <- LCVA::LCVA.pred(
    fit = fit,
    X_test = X_target,
    model = "C",
    Burn_in = burn_in,
    Nitr = max(2 * burn_in, 4000),
    return_likelihood = TRUE,
    verbose = FALSE
  )

  keep <- seq_len(dim(out$x_given_y_prob)[1])
  keep <- keep[keep > burn_in]
  posterior_phi <- apply(out$x_given_y_prob[keep, , , drop = FALSE], c(2,3), mean)

  # only guard real 0s
  posterior_phi[!is.finite(posterior_phi)] <- 0
  posterior_phi[posterior_phi == 0] <- .Machine$double.xmin

  list(
    posterior_phi = posterior_phi,  # N_target x C_site
    cause_ids = cause_ids,
    target_info = list(
      N = nrow(X_target),
      P = ncol(X_target),
      hash = target_hash
    )
  )
}
