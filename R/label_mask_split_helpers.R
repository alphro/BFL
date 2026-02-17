#' Balanced unlabeled-index generator (stratified, validity-checked)
#'
#' Selects indices to be treated as unlabeled (missing Y) at a target proportion.
#' Small classes (<= 10 samples) are handled carefully so that at least one
#' labeled example remains whenever possible. Remaining indices are sampled
#' from larger classes until a user-provided validity check passes.
#'
#' @param Y_true Vector of true class labels (length N). Can be factor/character/integer.
#' @param miss_prop Proportion of observations to mark as unlabeled (in (0,1)).
#' @param check_valid_fn Function with signature \code{function(Y_true, missing_idx)}
#'   returning a single logical indicating whether the sampled split is acceptable.
#' @param small_n Threshold for "small" classes (default 10).
#'
#' @return Integer vector of indices to be treated as unlabeled (missing Y).
#'
#' @keywords internal
split_idx_balanced <- function(Y_true, miss_prop, check_valid_fn, small_n = 10) {
  stopifnot(miss_prop > 0, miss_prop < 1)
  stopifnot(is.function(check_valid_fn))
  
  Y_true <- as.character(Y_true)
  N <- length(Y_true)
  num_samples <- round(miss_prop * N)
  
  category_counts <- table(Y_true)
  small_categories <- names(category_counts[category_counts <= small_n])
  
  missing_idx <- integer(0)
  
  # Small classes: mask ~miss_prop but keep at least one labeled if possible
  for (cod in small_categories) {
    cod_idx <- which(Y_true == cod)
    n_cod <- length(cod_idx)
    
    if (n_cod > 1) {
      n_miss <- min(round(miss_prop * n_cod), n_cod - 1L)
      if (n_miss > 0) {
        missing_idx <- c(missing_idx, sample(cod_idx, n_miss, replace = FALSE))
      }
    }
  }
  
  # Larger classes: sample remaining and enforce validity
  repeat {
    general_idx <- which(!Y_true %in% small_categories)
    
    need <- num_samples - length(missing_idx)
    need <- max(0L, need)
    
    general_missing <- if (need == 0L) integer(0) else sample(general_idx, need, replace = FALSE)
    final_missing <- c(missing_idx, general_missing)
    
    if (check_valid_fn(Y_true, final_missing)) break
  }
  
  final_missing
}


#' Label-shift split (no-overlap partition) via per-class Beta sampling
#'
#' Creates a labeled/unlabeled partition (no overlap) exhibiting label shift by:
#' sampling a per-class unlabeled proportion from \code{Beta(a,b)}, converting
#' to counts, then rescaling counts so the overall unlabeled size matches
#' \code{miss_prop * N}.
#'
#' @param Y_true Vector of true class labels (length N). Can be factor/character/integer.
#' @param miss_prop Proportion of observations to assign to the unlabeled set.
#' @param beta_shape Length-2 numeric giving \code{c(a,b)} for \code{Beta(a,b)} (default \code{c(0.2,0.2)}).
#'
#' @return A list with:
#' \describe{
#'   \item{labeled_idx}{Integer indices for labeled observations.}
#'   \item{unlabeled_idx}{Integer indices for unlabeled observations.}
#'   \item{pi_labeled}{Named numeric vector: empirical prevalence in labeled set.}
#'   \item{pi_unlabeled}{Named numeric vector: empirical prevalence in unlabeled set.}
#' }
#'
#' @keywords internal
split_idx_shift_beta <- function(Y_true, miss_prop, beta_shape = c(0.2, 0.2)) {
  stopifnot(miss_prop > 0, miss_prop < 1)
  stopifnot(length(beta_shape) == 2, all(beta_shape > 0))
  
  Y_true <- as.character(Y_true)
  cats <- sort(unique(Y_true))
  C <- length(cats)
  
  N <- length(Y_true)
  target_unlabeled <- round(miss_prop * N)
  
  # raw per-class unlabeled counts from Beta proportions
  unlab_raw <- numeric(C)
  for (i in seq_along(cats)) {
    idx <- which(Y_true == cats[i])
    p <- rbeta(1, beta_shape[1], beta_shape[2])
    unlab_raw[i] <- round(length(idx) * p)
  }
  
  # rescale to match target total
  scaling <- target_unlabeled / sum(unlab_raw)
  unlab_scaled <- round(unlab_raw * scaling)
  
  labeled_idx <- integer(0)
  unlabeled_idx <- integer(0)
  
  for (i in seq_along(cats)) {
    idx <- which(Y_true == cats[i])
    n_unlab <- min(length(idx), unlab_scaled[i])
    
    u <- if (n_unlab == 0) integer(0) else sample(idx, n_unlab, replace = FALSE)
    l <- setdiff(idx, u)
    
    labeled_idx <- c(labeled_idx, l)
    unlabeled_idx <- c(unlabeled_idx, u)
  }
  
  get_pi <- function(ii) {
    tab <- table(factor(Y_true[ii], levels = cats))
    prop <- as.numeric(tab / length(ii))
    names(prop) <- cats
    prop
  }
  
  list(
    labeled_idx   = labeled_idx,
    unlabeled_idx = unlabeled_idx,
    pi_labeled    = get_pi(labeled_idx),
    pi_unlabeled  = get_pi(unlabeled_idx)
  )
}


#' Label-shift split via Dirichlet + resampling (overlap allowed)
#'
#' Simulates label shift by drawing two prevalence vectors from a Dirichlet
#' distribution and then sampling labeled and unlabeled index sets via
#' multinomial class counts. Sampling is done with replacement within class,
#' so labeled/unlabeled sets may overlap and may contain duplicates.
#'
#' @param Y_true Vector of true class labels (length N). Can be factor/character/integer.
#' @param miss_prop Proportion of observations assigned to the labeled set size
#'   (mirrors the original implementation: labeled size = \code{miss_prop * N}).
#' @param dirichlet_alpha Positive scalar Dirichlet concentration parameter (default 1).
#' @param ensure_all_classes Logical; if TRUE, ensures every class appears at least once
#'   in the labeled set by adding any missing-class indices (default TRUE).
#'
#' @return A list with:
#' \describe{
#'   \item{labeled_idx}{Integer indices for labeled observations (may contain duplicates).}
#'   \item{unlabeled_idx}{Integer indices for unlabeled observations (may contain duplicates).}
#'   \item{pi_labeled}{Table-based prevalence in labeled set.}
#'   \item{pi_unlabeled}{Table-based prevalence in unlabeled set.}
#' }
#'
#' @keywords internal
split_idx_shift_dirichlet <- function(Y_true, miss_prop, dirichlet_alpha = 1, ensure_all_classes = TRUE) {
  stopifnot(miss_prop > 0, miss_prop < 1)
  stopifnot(dirichlet_alpha > 0)
  
  Y_true <- as.character(Y_true)
  cats <- unique(Y_true)
  C <- length(cats)
  N <- length(Y_true)
  
  # Dirichlet prevalence for labeled and unlabeled sets
  pi1 <- as.numeric(gtools::rdirichlet(1, rep(dirichlet_alpha, C)))
  pi2 <- as.numeric(gtools::rdirichlet(1, rep(dirichlet_alpha, C)))
  
  n_labeled <- round(miss_prop * N)
  n_unlabeled <- round((1 - miss_prop) * N)
  
  labeled_counts <- as.integer(rmultinom(1, size = n_labeled, prob = pi1))
  unlabeled_counts <- as.integer(rmultinom(1, size = n_unlabeled, prob = pi2))
  
  labeled_idx <- unlist(mapply(function(cat, n) {
    idx <- which(Y_true == cat)
    if (length(idx) == 0 || n == 0) return(integer(0))
    sample(idx, size = n, replace = TRUE)
  }, cats, labeled_counts))
  
  unlabeled_idx <- unlist(mapply(function(cat, n) {
    idx <- which(Y_true == cat)
    if (length(idx) == 0 || n == 0) return(integer(0))
    sample(idx, size = n, replace = TRUE)
  }, cats, unlabeled_counts))
  
  if (isTRUE(ensure_all_classes)) {
    cats_in_labeled <- unique(Y_true[labeled_idx])
    missing_cats <- setdiff(cats, cats_in_labeled)
    if (length(missing_cats) > 0) {
      add_idx <- unlist(lapply(missing_cats, function(cat) which(Y_true == cat)))
      labeled_idx <- c(labeled_idx, add_idx)
    }
  }
  
  pi_labeled <- table(Y_true[labeled_idx]) / length(labeled_idx)
  pi_unlabeled <- table(Y_true[unlabeled_idx]) / length(unlabeled_idx)
  
  list(
    labeled_idx   = labeled_idx,
    unlabeled_idx = unlabeled_idx,
    pi_labeled    = pi_labeled,
    pi_unlabeled  = pi_unlabeled
  )
}

#' Validity check for balanced missing-label splits
#'
#' Ensures that every class appearing in the missing (unlabeled) subset also
#' appears at least once in the observed (labeled) subset. This prevents any
#' class from being entirely unlabeled.
#'
#' @param Y_true Vector of true labels (length N).
#' @param missing_idx Integer indices treated as missing/unlabeled.
#'
#' @return Logical; TRUE if split is valid, FALSE otherwise.
#'
#' @keywords internal
check_valid_split <- function(Y_true, missing_idx) {
  Y_true <- as.character(Y_true)
  observed <- Y_true[-missing_idx]
  missing  <- Y_true[ missing_idx]
  
  observed_cats <- unique(observed[!is.na(observed)])
  missing_cats  <- unique(missing)
  
  all(missing_cats %in% observed_cats)
}