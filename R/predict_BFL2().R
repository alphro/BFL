# ---- helpers (fixed versions of theirs) ----

mode_int <- function(x) {
  tab <- table(x)
  as.integer(names(tab)[which.max(tab)])
}

top1_from_samples <- function(Y_samp_mat) {
  # Y_samp_mat: S x N
  apply(Y_samp_mat, 2, mode_int)
}

top1_acc <- function(Y_top1, Y_true_int) {
  mean(Y_top1 == Y_true_int)
}

CSMF_acc <- function(pi_hat, pi_true) {
  pi_hat <- as.numeric(pi_hat); pi_true <- as.numeric(pi_true)
  if (length(pi_hat) != length(pi_true)) stop("pi_hat and pi_true must have same length")
  if (sum(pi_hat) > 1.0001 || sum(pi_true) > 1.0001) stop("pi vectors must sum to 1")
  C <- length(pi_true)
  1 - (sum(abs(pi_hat - pi_true)) / (2 * (1 - min(pi_true))))
}

# ---- core: paper-style posterior predictive sampling ----
posterior_predict_BFL <- function(phi, lambda, pi, ndraws = NULL, seed = NULL) {
  # phi: N x C x M
  # lambda: S x C x M
  # pi: S x C

  stopifnot(length(dim(phi)) == 3)
  stopifnot(length(dim(lambda)) == 3)
  stopifnot(length(dim(pi)) == 2)

  N <- dim(phi)[1]; C <- dim(phi)[2]; M <- dim(phi)[3]
  S <- dim(lambda)[1]

  stopifnot(dim(lambda)[2] == C, dim(lambda)[3] == M)
  stopifnot(dim(pi)[2] == C)

  draws <- seq_len(S)
  if (!is.null(ndraws)) {
    stopifnot(ndraws >= 1, ndraws <= S)
    if (!is.null(seed)) set.seed(seed)
    draws <- sample(draws, size = ndraws, replace = FALSE)
  }

  S_use <- length(draws)

  # store sampled labels: S_use x N
  Y_samp <- matrix(NA_integer_, nrow = S_use, ncol = N)

  # optional: posterior mean probs over draws (N x C)
  prob_mean <- matrix(0, nrow = N, ncol = C)

  for (ii in seq_len(S_use)) {
    s <- draws[ii]
    lam_s <- lambda[s, , , drop = FALSE]  # 1 x C x M
    lam_s <- matrix(lam_s, nrow = C, ncol = M)
    pi_s <- as.numeric(pi[s, ])

    for (i in seq_len(N)) {
      # phi_i: C x M
      phi_i <- phi[i, , , drop = FALSE]
      phi_i <- matrix(phi_i, nrow = C, ncol = M)

      # score[c] = pi_s[c] * sum_m phi[i,c,m] * lam_s[c,m]
      score <- pi_s * rowSums(phi_i * lam_s)

      # numeric safety
      score[!is.finite(score)] <- 0
      if (all(score <= 0)) score <- rep(1, C)

      p <- score / sum(score)

      prob_mean[i, ] <- prob_mean[i, ] + p
      Y_samp[ii, i] <- sample.int(C, size = 1, prob = p)
    }
  }

  prob_mean <- prob_mean / S_use
  list(Y_samp = Y_samp, prob_mean = prob_mean, draws_used = draws)
}

align_prevalence_to_master <- function(pi_vec, causes_vec, master_causes) {
  out <- setNames(rep(0, length(master_causes)), master_causes)
  out[as.character(causes_vec)] <- as.numeric(pi_vec)
  out
}

#' Emulates the original BFL evaluation procedure (Zoey-aligned)
#'
#' - Top-1: posterior predictive sampling, mode across draws
#' - CSMF: Zoey-style prevalence estimator (with labeled/unlabeled blend when provided)
#' - Prevalence vectors are aligned to the global cause support before scoring
#'
#' @export
predict_BFL2 <- function(global_fit,
                         ndraws = NULL,
                         seed = 12345,
                         Y_true = NULL,
                         Y_partial = NULL) {

  phi <- global_fit$phi          # N x C x M
  lambda <- global_fit$lambda    # S x C x M
  pi <- global_fit$pi            # S x C
  causes <- global_fit$causes    # global cause labels (master support)
  C <- length(causes)

  # ---- posterior predictive sampling for Top-1 ----
  pp <- posterior_predict_BFL(phi = phi, lambda = lambda, pi = pi,
                              ndraws = ndraws, seed = seed)

  Y_top1_idx <- top1_from_samples(pp$Y_samp)
  Y_top1 <- causes[Y_top1_idx]

  out <- list(
    Y_top1 = Y_top1,
    Y_top1_idx = Y_top1_idx,
    prob_mean = pp$prob_mean,
    causes = causes,
    draws_used = pp$draws_used
  )

  # ---- If no truth provided, we can't score ----
  if (is.null(Y_true)) return(out)

  # Map Y_true into global indices (Zoey uses integers; you can pass chars or ints)
  if (is.character(Y_true)) {
    cause_to_idx <- setNames(seq_len(C), causes)
    Y_true_idx <- unname(cause_to_idx[as.character(Y_true)])
    if (anyNA(Y_true_idx)) stop("Y_true contains labels not in global_fit$causes")
  } else {
    Y_true_idx <- as.integer(Y_true)
  }

  # Top-1 accuracy (same as Zoey)
  out$acc_top1 <- mean(Y_top1_idx == Y_true_idx)

  # ---- TRUE prevalence over the *full population* you pass in Y_true ----
  tab_true <- table(factor(Y_true_idx, levels = seq_len(C)))
  pi_true <- as.numeric(tab_true) / sum(tab_true)
  pi_true_named <- setNames(pi_true, causes)

  # ---- Posterior mean prevalence from Stan pi draws ----
  pi_hat <- colMeans(pi)
  pi_hat_named <- setNames(pi_hat, causes)

  # Align (mostly a no-op since both are already on global causes,
  # but keeps us faithful and future-proof)
  pi_hat_aligned <- align_prevalence_to_master(pi_hat_named, names(pi_hat_named), causes)
  pi_true_aligned <- align_prevalence_to_master(pi_true_named, names(pi_true_named), causes)

  out$pi_true <- pi_true_aligned
  out$pi_hat <- pi_hat_aligned
  out$csmf_acc <- CSMF_acc(pi_hat_aligned, pi_true_aligned)

  # ---- Zoey-style blended prevalence when partial labels exist ----
  # Y_partial: length N, with NA for unlabeled
  if (!is.null(Y_partial)) {
    if (length(Y_partial) != length(Y_true)) {
      stop("Y_partial must have same length as Y_true (full population).")
    }

    labeled_mask <- !is.na(Y_partial)
    n0 <- length(Y_partial)
    nU <- sum(!labeled_mask)

    # labeled counts nLc in global cause space
    if (is.character(Y_partial)) {
      ylab_idx <- unname(cause_to_idx[as.character(Y_partial[labeled_mask])])
      if (anyNA(ylab_idx)) stop("Y_partial contains labels not in global_fit$causes")
    } else {
      ylab_idx <- as.integer(Y_partial[labeled_mask])
    }

    nLc <- as.numeric(table(factor(ylab_idx, levels = seq_len(C))))

    pi_blend <- (nU * pi_hat + nLc) / n0
    pi_blend_named <- setNames(pi_blend, causes)

    # Align (again, mostly a no-op, but consistent)
    pi_blend_aligned <- align_prevalence_to_master(pi_blend_named, names(pi_blend_named), causes)

    out$pi_blend <- pi_blend_aligned
    out$csmf_acc_blend <- CSMF_acc(pi_blend_aligned, pi_true_aligned)
  }

  out
}
