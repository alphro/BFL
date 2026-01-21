filter_sparse_causes <- function(sim_data, threshold = 0) {

  # calculate frequency of each cause of death
  cause_freq <- table(sim_data$data.truth$Y.t)

  # identify causes with samples <= threshold
  sparse_causes <- names(cause_freq[cause_freq <= threshold])
  sparse_causes <- as.numeric(sparse_causes)

  # keep samples NOT in sparse causes
  keep_indices <- which(!(sim_data$data.truth$Y.t %in% sparse_causes))

  # filter the data
  sim_data$data.truth$Y.t <- sim_data$data.truth$Y.t[keep_indices]
  sim_data$data$Y.t       <- sim_data$data$Y.t[keep_indices]
  sim_data$data$X         <- sim_data$data$X[keep_indices, , drop = FALSE]
  sim_data$data$G         <- sim_data$data$G[keep_indices]

  # re-label causes to 1:C (continuous)
  unique_causes <- sort(unique(sim_data$data.truth$Y.t))
  new_labels <- seq_along(unique_causes)

  mapping_origin_to_new <- setNames(new_labels, unique_causes)
  mapping_new_to_origin <- setNames(unique_causes, new_labels)

  sim_data$data.truth$Y.t <- as.integer(mapping_origin_to_new[as.character(sim_data$data.truth$Y.t)])
  sim_data$data$Y.t       <- as.integer(mapping_origin_to_new[as.character(sim_data$data$Y.t)])

  # update C
  sim_data$C <- length(unique_causes)

  list(
    filtered_data = sim_data,
    mapping_origin_to_new = mapping_origin_to_new,
    mapping_new_to_origin = mapping_new_to_origin
  )
}


# checking valid_sample
check_valid_sample <- function(full_data, missing_indices) {
  observed_data <- full_data[-missing_indices]
  missing_data  <- full_data[missing_indices]

  observed_categories <- unique(observed_data[!is.na(observed_data)])
  missing_categories  <- unique(missing_data)

  all(missing_categories %in% observed_categories)
}


# generate missing Yt
generate_missing_Yt <- function(Y.t, miss_prop) {

  #so let's say missing proprtion is 1
  #then length(Y.t) is the entire length then
  num_samples <- round(miss_prop * length(Y.t))

  #for categorical_counts, here this is the count of each individual cod
  category_counts <- table(Y.t)
  #small categories is any which one where the category_counts is less than or equal to 10
  #these are indices
  small_categories <- as.numeric(names(category_counts[category_counts <= 10]))

  #missing indcies is empty for now
  missing_indices <- integer(0)

  # handle small categories: never remove all of any category
  for (cod in small_categories) {
    #lets track first small category (index 31)
    #this is grabbing global Y.t having cod
    cod_indices <- which(Y.t == cod)

    #so how many small samples are there
    n_cod <- length(cod_indices)

    if (n_cod > 1) {
      # this here is saying that im going to let the # that is missing equal to whicever is smaller
      # either miss_prop * n_cod
        # so larger miss_prop equates to larger missing # and smaller miss_prop equate to smaller missing #
        # n_cod-1 (so we don't want the exact amount? hm. interesting. I wonder why)
      num_miss <- min(round(miss_prop * n_cod), n_cod - 1)

      #here it's saying from those small samples, let's get the # of missing_samples
      #so let's say it was 5 missing cod_indic,es we are now only grabbing 4 for some reason.
      #if miss_prop is small, we are grabbing 1 for example?
      small_miss_indices <- sample(cod_indices, num_miss, replace = FALSE)
        # so here it is saying that it will take the cod_indices, and then it will sample
      # the number of missing, so here if miss_prop =0, then we are saying basially none is misisng
      # if missing_prop is one, we are saying sample all of the missing ones
      # so for some reason, we are not going to sample all of it
      missing_indices <- c(missing_indices, small_miss_indices)
      # and then here it is creating miss_indices,
      #here, it's saying we are going to add on the small_missing_indices for the sample
    }
  }

  repeat {
    #here, we are going to pick all of the samples of Y.t which are not in small_categories
    #so full matrix
    general_indices <- which(!(Y.t %in% small_categories))
    #and then we are saying that we have that the num_sample we are grabbing is the num of missing_prop * the length of Y.t
    #so smaller miss_prop means smaller remaining samples
    #we also subtrac length of missing indices.
    remaining <- num_samples - length(missing_indices) # so this function is essentially grabbing the # of indices wihtout small categories
    # and here, we are saying the missing_indices are really just controlling teh

    general_missing <- if (remaining > 0) {
      sample(general_indices, remaining, replace = FALSE)
    } else {
      integer(0)
    }

    final_missing_indices <- c(missing_indices, general_missing)

    if (check_valid_sample(Y.t, final_missing_indices)) break
  }

  final_missing_indices
}


###### WILL FILL OUT LATER #######
generate_train_test_Yt_unbalanced <- function(Y.t, miss_prop) {}
generate_missing_Yt_unbalanced <- function(Y.t, miss_prop) {}
#######


get_observed_sim_data <- function(sim_data) {
  obs_data <- list()
  obs_data$data.truth <- list()

  keep <- which(!is.na(sim_data$data$Y.t))

  obs_data$data.truth$Y.t <- sim_data$data$Y.t[keep]
  obs_data$data$X         <- sim_data$data$X[keep, , drop = FALSE]
  obs_data$data$G         <- rep(1L, length(obs_data$data.truth$Y.t))
  obs_data$NG             <- 1L

  unique_causes <- sort(unique(obs_data$data.truth$Y.t))
  obs_data$C <- length(unique_causes)

  list(obs_data = obs_data)
}


get_global_post_pred_Y <- function(posterior_phi,
                                   posterior_lambda,
                                   posterior_pi) {

  n_samples <- dim(posterior_lambda)[1]
  C <- dim(posterior_phi)[2]

  posterior_pred_Y_prob <- sapply(
    1:n_samples,
    function(s) {
      apply(posterior_phi, c(1, 2), function(phi_icm) sum(phi_icm * posterior_lambda[s, ])) *
        posterior_pi[s, ]
    },
    simplify = "array"
  )

  # normalize across C for each (sample, i)
  posterior_pred_Y_prob <- apply(posterior_pred_Y_prob, c(1, 3), function(x) x / sum(x))
  posterior_pred_Y <- apply(posterior_pred_Y_prob, c(2, 3), function(x) sample(1:C, 1, prob = x))

  list(
    posterior_pred_Y_prob = posterior_pred_Y_prob,
    posterior_pred_Y = t(posterior_pred_Y)
  )
}


get_sim_data_from_phmrc_by_location <- function(site, phmrc = NULL) {

  # If user didn't pass a data.frame, load toy CSV from inst/extdata/
  if (is.null(phmrc)) {
    phmrc_path <- system.file("extdata", "phmrc_clean.csv", package = "BFL")
    if (phmrc_path == "") stop("Example dataset not found in package: phmrc_clean.csv")
    phmrc <- read.csv(phmrc_path, stringsAsFactors = FALSE)
  }

  causes <- unique(phmrc$cause)

  phmrc_site <- phmrc[phmrc$site == site, , drop = FALSE]

  phmrc_site$cause <- factor(phmrc_site$cause, levels = causes)
  phmrc_site$Y <- as.integer(phmrc_site$cause)

  data <- list()
  data$cause_level <- causes
  data$cause <- phmrc_site$cause
  data$Y.t <- phmrc_site$Y
  data$X <- phmrc_site[, 1:168, drop = FALSE]
  data$G <- rep(1L, length(phmrc_site$Y))

  data.truth <- list()
  data.truth$Y.t <- phmrc_site$Y

  sim_data <- list()
  sim_data$data <- data
  sim_data$data.truth <- data.truth
  sim_data$C <- length(unique(phmrc_site$Y))
  sim_data$NG <- 1L

  sim_data
}


get_acc <- function(pred_Yt, true_Yt) {
  acc <- rep(NA_real_, ncol(pred_Yt))

  for (t in 1:ncol(pred_Yt)) {
    cause_frequencies <- table(pred_Yt[, t])
    most_probable_cause <- names(which.max(cause_frequencies))
    most_probable_cause_number <- as.integer(most_probable_cause)

    acc[t] <- ifelse(most_probable_cause_number == true_Yt[t], 1, 0)
  }

  mean(acc)
}


get_Y_pred <- function(pred_Yt) {
  Y <- rep(NA_integer_, ncol(pred_Yt))

  for (t in 1:ncol(pred_Yt)) {
    cause_frequencies <- table(pred_Yt[, t])
    most_probable_cause <- names(which.max(cause_frequencies))
    Y[t] <- as.integer(most_probable_cause)
  }

  Y
}


align_prevalence <- function(prevalence_A, mapping_A, prevalence_B, mapping_B) {
  causes_A <- unname(mapping_A)
  causes_B <- unname(mapping_B)

  all_causes <- sort(unique(c(causes_A, causes_B)))

  aligned_A <- setNames(rep(0, length(all_causes)), all_causes)
  aligned_B <- setNames(rep(0, length(all_causes)), all_causes)

  aligned_A[as.character(causes_A)] <- prevalence_A
  aligned_B[as.character(causes_B)] <- prevalence_B

  list(aligned_A = aligned_A, aligned_B = aligned_B)
}


CSMF_acc <- function(pi_hat, pi) {
  if (length(pi_hat) != length(pi) || sum(pi_hat) > 1.001 || sum(pi) > 1.001) {
    stop("Invalid")
  }

  csmf_accuracy <- 1 - (sum(abs(pi_hat - pi)) / (2 * (1 - min(pi))))
  csmf_accuracy
}


# split domain and partial; no need yet
split_domain_and_partial <- function(sim_data) {
  # TODO
}
