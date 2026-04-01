# run_BFL_yz.R
#
# BFL with P(X|Y,Z) inputs — K x M base models.
#
# Instead of marginalising Z before passing to BFL:
#   phi_{ic}^(m) = sum_k theta_{ck} * P(X_i | Z=k)   [standard, N x C]
#
# This passes:
#   phi_{ick}^(m) = theta_{ck} * P(X_i | Z=k)         [K x M entries, N x C each]
#
# BFL's Stan model is unchanged — it just sees K*M "sites".
# Limitation: latent classes are not label-aligned across sites.

# ---- 1. Extract P(X|Y,Z) from one LCVA fit ----
#
# Returns phi_yz (N x C x K): posterior-mean theta-weighted likelihoods.
#
# LCVA phi layouts handled:
#   (a) per-draw 2D:  phi[[s]]     is K x P  (or P x K)
#   (b) per-chain 3D: phi[[chain]] is S_chain x K x P
#   (c) per-draw 3D:  phi[[s]]     is C x K x P  (cause-specific)
#
# Heuristic: if length(phi_list) == length(loglam_list), both are per-draw.
# A 3D phi in that case must be layout (c); otherwise layout (b).
extract_pxyz <- function(lcva_fit, X_target, cause_ids, burn_in = 200) {

  X_target <- as.matrix(X_target)
  storage.mode(X_target) <- "numeric"
  N <- nrow(X_target); P <- ncol(X_target); C <- length(cause_ids)

  # Unwrap one-element wrapper lists
  while (is.list(lcva_fit) && length(lcva_fit) == 1L &&
         (is.null(lcva_fit$phi) || is.null(lcva_fit$loglambda)))
    lcva_fit <- lcva_fit[[1]]

  phi_list    <- lcva_fit$phi
  loglam_list <- lcva_fit$loglambda
  if (is.null(phi_list) || is.null(loglam_list))
    stop("Could not find phi / loglambda in LCVA fit.\n",
         "Available fields: ", paste(names(lcva_fit), collapse = ", "))

  first_phi <- phi_list[[1]]
  phi_ndim  <- length(dim(first_phi))
  first_lam <- loglam_list[[1]]
  lam_ndim  <- length(dim(first_lam))
  same_len  <- (length(phi_list) == length(loglam_list))

  # ---- Stack phi ----
  if (phi_ndim == 2L) {
    # Layout (a): per-draw K x P (or P x K)
    S  <- length(phi_list)
    d1 <- dim(first_phi)[1]; d2 <- dim(first_phi)[2]
    if (d2 == P) {
      K <- d1
      phi_kps <- array(NA_real_, c(K, P, S))
      for (s in seq_len(S)) phi_kps[,,s] <- phi_list[[s]]
    } else if (d1 == P) {
      K <- d2
      phi_kps <- array(NA_real_, c(K, P, S))
      for (s in seq_len(S)) phi_kps[,,s] <- t(phi_list[[s]])
    } else {
      stop("phi[[1]] dim (", paste(dim(first_phi), collapse=","),
           ") does not match P=", P)
    }
    phi_layout <- "kp"

  } else if (phi_ndim == 3L && same_len) {
    # Layout (c): per-draw C x K x P (cause-specific phi)
    dp <- dim(first_phi)
    if (dp[1] == C && dp[3] == P) {
      K <- dp[2]
    } else if (dp[1] == C && dp[2] == P) {
      K <- dp[3]
      phi_list <- lapply(phi_list, function(ph) aperm(ph, c(1L, 3L, 2L)))  # -> C x K x P
    } else {
      stop("phi[[1]] dim (", paste(dp, collapse=","),
           ") does not match C x K x P with C=", C, " P=", P)
    }
    S <- length(phi_list)
    # Stack: C x K x P x S
    phi_ckps <- array(NA_real_, c(C, K, P, S))
    for (s in seq_len(S)) phi_ckps[,,,s] <- phi_list[[s]]
    phi_layout <- "ckp"

  } else if (phi_ndim == 3L && !same_len) {
    # Layout (b): per-chain S_chain x K x P
    dp <- dim(first_phi)
    if (dp[3] == P) {
      K <- dp[2]
    } else if (dp[2] == P) {
      K <- dp[3]
      phi_list <- lapply(phi_list, function(ch) aperm(ch, c(1L, 3L, 2L)))
    } else {
      stop("phi[[chain]] dim (", paste(dp, collapse=","),
           ") does not match P=", P)
    }
    S <- sum(sapply(phi_list, function(ch) dim(ch)[1]))
    phi_kps <- array(NA_real_, c(K, P, S))
    idx <- 1L
    for (ch in phi_list) {
      n_s <- dim(ch)[1]
      for (s in seq_len(n_s)) { phi_kps[,, idx] <- ch[s,,]; idx <- idx + 1L }
    }
    phi_layout <- "kp"

  } else {
    stop("Unexpected phi structure: phi[[1]] has ", phi_ndim, " dimensions.")
  }

  # ---- Stack loglambda -> C x K x S_total (softmax -> theta) ----
  if (lam_ndim == 2L || (lam_ndim == 3L && same_len)) {
    # Per-draw layout
    lam_S <- length(loglam_list)
    theta_chain <- array(NA_real_, c(C, K, lam_S))
    for (s in seq_len(lam_S)) {
      loglam_s <- loglam_list[[s]]
      if (length(dim(loglam_s)) == 3L) {
        if (dim(loglam_s)[3] == 1L) dim(loglam_s) <- dim(loglam_s)[1:2]
        else stop("loglambda[[", s, "]] is 3D with unexpected dims: ",
                  paste(dim(loglam_s), collapse=","))
      }
      if (!all(dim(loglam_s) == c(C, K)))
        stop("loglambda[[", s, "]] dim (", paste(dim(loglam_s), collapse=","),
             ") expected (", C, ",", K, ")")
      ex <- exp(loglam_s - apply(loglam_s, 1, max))
      theta_chain[,,s] <- ex / rowSums(ex)
    }
    S_lam <- lam_S

  } else if (lam_ndim == 3L && !same_len) {
    # Per-chain S_chain x C x K
    lam_S <- sum(sapply(loglam_list, function(ch) dim(ch)[1]))
    theta_chain <- array(NA_real_, c(C, K, lam_S))
    idx <- 1L
    for (ch in loglam_list) {
      n_s <- dim(ch)[1]
      for (s in seq_len(n_s)) {
        loglam_s <- ch[s,,]   # C x K
        ex <- exp(loglam_s - apply(loglam_s, 1, max))
        theta_chain[,, idx] <- ex / rowSums(ex)
        idx <- idx + 1L
      }
    }
    S_lam <- lam_S

  } else {
    stop("Unexpected loglambda structure: loglambda[[1]] has ", lam_ndim, " dimensions.")
  }

  # Reconcile draw counts
  S_phi <- if (phi_layout == "kp") dim(phi_kps)[3] else dim(phi_ckps)[4]
  if (S_phi != S_lam)
    stop("phi draws (", S_phi, ") != loglambda draws (", S_lam, ")")
  S <- S_phi

  # Burn-in
  keep   <- seq.int(burn_in + 1L, S)
  if (length(keep) < 10L)
    stop("burn_in=", burn_in, " leaves <10 draws (S=", S, ")")
  S_keep <- length(keep)

  theta_k <- theta_chain[,, keep, drop = FALSE]   # C x K x S_keep

  # ---- Compute phi_yz (N x C x K), log-normalised across causes ----
  #
  # Raw values  theta[c,k] * P(X_i|Z=k)  underflow to zero (products of 168
  # binary probs are ~10^-50).  We instead compute log-unnormalised scores and
  # normalise across causes for each (i, k), yielding values on [0,1] that
  # match the scale BFL's Stan model expects (analogous to BFL2's posterior_phi).
  #
  # phi_yz[i,c,k] = exp( log_theta[c,k] + log_px[i,c,k]
  #                      - logsumexp_c( log_theta[c',k] + log_px[i,c',k] ) )
  #
  logsumexp_rows <- function(mat) {   # mat: N x C  ->  N-vector
    m <- apply(mat, 1, max)
    m + log(rowSums(exp(mat - m)))
  }

  phi_yz <- array(0, c(N, C, K))

  if (phi_layout == "kp") {
    phi_k <- phi_kps[,, keep, drop = FALSE]   # K x P x S_keep
    for (t in seq_len(S_keep)) {
      phi_t <- pmin(pmax(phi_k[,,t], .Machine$double.xmin), 1 - .Machine$double.xmin)
      log_px_z <- X_target %*% t(log(phi_t)) +
                  (1 - X_target) %*% t(log(1 - phi_t))   # N x K
      for (k in seq_len(K)) {
        log_theta_k  <- log(pmax(theta_k[, k, t], .Machine$double.xmin))  # C
        log_unnorm   <- outer(log_px_z[, k], log_theta_k, "+")             # N x C
        log_z        <- logsumexp_rows(log_unnorm)
        phi_yz[,, k] <- phi_yz[,, k] + exp(log_unnorm - log_z)
      }
    }

  } else {
    phi_k <- phi_ckps[,,, keep, drop = FALSE]   # C x K x P x S_keep
    for (t in seq_len(S_keep)) {
      for (k in seq_len(K)) {
        log_theta_k <- log(pmax(theta_k[, k, t], .Machine$double.xmin))   # C
        log_unnorm  <- matrix(0, N, C)
        for (c_idx in seq_len(C)) {
          phi_ck <- pmin(pmax(phi_k[c_idx, k, , t], .Machine$double.xmin),
                         1 - .Machine$double.xmin)   # P
          log_unnorm[, c_idx] <-
            as.vector(X_target %*% log(phi_ck) +
                      (1 - X_target) %*% log(1 - phi_ck)) + log_theta_k[c_idx]
        }
        log_z        <- logsumexp_rows(log_unnorm)
        phi_yz[,, k] <- phi_yz[,, k] + exp(log_unnorm - log_z)
      }
    }
  }

  phi_yz <- phi_yz / S_keep

  # Floor at eps then renormalise each (i,k) slice.
  # Raw log-likelihoods span >700 log-units across 168 binary features, so
  # many exp(log_norm) values underflow to 0.  Clipping to .Machine$double.xmin
  # leaves phi ~ 2e-308, which makes 1/total_prob_i ~ 5e307 in Stan's autodiff
  # and produces NaN gradients.  A floor of 1e-6 keeps values discriminative
  # while staying far from double-precision limits.
  eps <- 1e-6
  for (k in seq_len(K)) {
    mat <- phi_yz[,, k]                         # N x C
    mat[!is.finite(mat) | mat < eps] <- eps
    phi_yz[,, k] <- mat / rowSums(mat)          # renormalise -> sums to 1 across C
  }

  list(phi_yz = phi_yz, K = K, C = C,
       cause_ids = cause_ids, row_hash = compute_row_hashes(X_target))
}


# ---- 2. Build K x M local_summaries ----
#
# Each (site, latent class k) pair becomes one entry named "<site>__k<k>".
build_local_summaries_yz <- function(pxyz_list, P_target) {
  ls_out <- list()
  for (site in names(pxyz_list)) {
    pxyz <- pxyz_list[[site]]
    for (k in seq_len(pxyz$K)) {
      phi_k         <- pxyz$phi_yz[,,k]          # N x C
      ls_out[[sprintf("%s__k%d", site, k)]] <- list(
        posterior_phi = phi_k,
        cause_ids     = pxyz$cause_ids,
        target_info   = list(row_hash = pxyz$row_hash,
                             N = nrow(phi_k), P = P_target)
      )
    }
  }
  ls_out
}


# ---- 3. Unified entry point ----
#
# Extracts P(X|Y,Z) per site, builds K*M local_summaries, calls run_BFL().
# Returns a standard BFL fit — use score_BFL() as normal.
run_BFL_yz <- function(lcva_fits,
                        cause_ids_list,
                        X_target,
                        Y_target    = NULL,
                        burn_in     = 200,
                        label_shift = FALSE,
                        stan_args   = list(iter = 2000, chains = 4, seed = 12345)) {

  stopifnot(is.list(lcva_fits), !is.null(names(lcva_fits)))
  stopifnot(identical(sort(names(lcva_fits)), sort(names(cause_ids_list))))

  X_target <- as.matrix(X_target)
  P        <- ncol(X_target)

  cat("Extracting P(X|Y,Z) from LCVA fits...\n")
  pxyz_list <- setNames(lapply(names(lcva_fits), function(s) {
    cat(sprintf("  [%s] ", s))
    out <- extract_pxyz(lcva_fits[[s]], X_target, cause_ids_list[[s]], burn_in)
    cat(sprintf("K=%d C=%d N=%d\n", out$K, out$C, nrow(out$phi_yz)))
    out
  }), names(lcva_fits))

  ls_yz <- build_local_summaries_yz(pxyz_list, P_target = P)
  cat(sprintf("Total base models: %d\n\n", length(ls_yz)))

  run_BFL(local_summaries = ls_yz, X_target = X_target,
          Y_target = Y_target, label_shift = label_shift, stan_args = stan_args)
}
