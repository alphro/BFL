#' Create label-missingness mask for BFL variants
#'
#' Generates label visibility patterns and index sets used by BFL model variants.
#' This function masks a proportion of labels (sets them to NA) and returns
#' indices needed for Base, Domain, Partial, and Mix configurations.
#'
#' @param Y_true True labels (length N). Factor/character/integer.
#' @param model One of "Base", "Domain", "Partial", "Mix".
#' @param miss_prop Proportion of observations to treat as unlabeled.
#'   Typical value is 0.80 for 20% labeled / 80% unlabeled.
#' @param mix_domain_prop_within_labeled For Mix only: proportion of labeled rows
#'   assigned to the domain subset (excluded from Stan).
#' @param seed Optional integer seed for reproducible splitting.
#' @param check_valid_fn Optional validity check for balanced splitting.
#'   Signature: function(Y_true, missing_idx) returning TRUE/FALSE.
#'
#' @return A list containing model, Y_obs, Y_true, missing_idx, stan_idx,
#'   and (for some models) labeled_idx, domain_idx, partial_idx.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' Y <- sample(c("A","B","C"), 30, replace = TRUE)
#'
#' # ---- Base: no labels observed; Stan uses all rows ----
#' s_base <- make_label_mask(Y, model = "Base")
#' sum(!is.na(s_base$Y_obs))           # 0
#' length(s_base$stan_idx)             # N
#'
#' # ---- Domain: Stan sees ONLY unlabeled rows ----
#' s_dom <- make_label_mask(Y, model = "Domain", miss_prop = 0.8, seed = 1)
#' length(s_dom$labeled_idx)           # ~20% of N
#' length(s_dom$missing_idx)           # ~80% of N
#' identical(s_dom$stan_idx, s_dom$missing_idx)  # TRUE
#'
#' # ---- Partial: Stan sees ALL rows (labeled + unlabeled) ----
#' s_par <- make_label_mask(Y, model = "Partial", miss_prop = 0.8, seed = 1)
#' length(s_par$stan_idx)              # N
#' sum(is.na(s_par$Y_obs))             # ~80% of N
#'
#' # ---- Mix: like Partial, but remove a labeled "domain" subset from Stan ----
#' s_mix <- make_label_mask(
#'   Y,
#'   model = "Mix",
#'   miss_prop = 0.8,
#'   mix_domain_prop_within_labeled = 0.5,
#'   seed = 1
#' )
#' # domain_idx excluded from Stan
#' intersect(s_mix$domain_idx, s_mix$stan_idx)   # integer(0)
#' # stan_idx = partial_idx + missing_idx
#' setequal(s_mix$stan_idx, c(s_mix$partial_idx, s_mix$missing_idx))  # TRUE
#' }
#'
#' @export
make_label_mask <- function(
    Y_true,
    model = c("Base","Domain","Partial","Mix"),
    miss_prop = 0.80,
    mix_domain_prop_within_labeled = 0.50,
    seed = NULL,
    check_valid_fn = NULL
) {

  if (anyNA(Y_true)) {
    stop("Y_true must not contain NA. Pass fully observed labels to make_label_mask().")
  }

  model <- match.arg(model)
  if (!is.null(seed)) set.seed(seed)

  N <- length(Y_true)
  Y_true_chr <- as.character(Y_true)

  if (is.null(check_valid_fn)) {
    check_valid_fn <- function(Y_true, missing_idx) {
      observed <- Y_true[-missing_idx]
      missing  <- Y_true[ missing_idx]
      all(unique(missing) %in% unique(observed))
    }
  }

  if (model == "Base") {
    return(list(
      model = model,
      Y_true = Y_true_chr,
      Y_obs = rep(NA_character_, N),
      missing_idx = seq_len(N),
      labeled_idx = integer(0),
      stan_idx = seq_len(N)
    ))
  }

  missing_idx <- split_idx_balanced(
    Y_true = Y_true_chr,
    miss_prop = miss_prop,
    check_valid_fn = check_valid_fn,
    small_n = 10
  )

  missing_idx <- sort(unique(missing_idx))
  labeled_idx <- sort(setdiff(seq_len(N), missing_idx))

  Y_obs <- Y_true_chr
  Y_obs[missing_idx] <- NA_character_

  if (model == "Domain") {
    return(list(
      model = model,
      Y_true = Y_true_chr,
      Y_obs = Y_obs,
      labeled_idx = labeled_idx,
      missing_idx = missing_idx,
      stan_idx = missing_idx
    ))
  }

  if (model == "Partial") {
    return(list(
      model = model,
      Y_true = Y_true_chr,
      Y_obs = Y_obs,
      labeled_idx = labeled_idx,
      missing_idx = missing_idx,
      stan_idx = seq_len(N)
    ))
  }

  if (model == "Mix") {

    nL <- length(labeled_idx)
    nD <- if (nL == 0L) 0L else max(1L, floor(mix_domain_prop_within_labeled * nL))

    domain_idx  <- if (nD == 0L) integer(0) else sort(sample(labeled_idx, nD))
    partial_idx <- sort(setdiff(labeled_idx, domain_idx))

    stan_idx <- sort(c(partial_idx, missing_idx))

    return(list(
      model = model,
      Y_true = Y_true_chr,
      Y_obs = Y_obs,
      labeled_idx = labeled_idx,
      domain_idx = domain_idx,
      partial_idx = partial_idx,
      missing_idx = missing_idx,
      stan_idx = stan_idx
    ))
  }
}
