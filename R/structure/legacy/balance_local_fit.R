#' Execute balanced local fit
#'
#' @param test_site character
#' @param sites character vector
#' @param K integer
#' @param miss_prop numeric
#' @return list
#' @export
execute_balance_local_fit <- function(test_site, sites, K = 5, miss_prop) {

  sim_data_list <- sim_data_filtered_list <- list()
  LCVA_local_model_fit <- list()
  LCVA_local_model_test_obs_fit <- list()
  posterior_phi_full <- list()
  confusion_local_parafac_self_obs <- list()

  adjust_est_prev <- cause_pred <- list()

  CSMF_local_parafac_self_obs <- rep(NA, length(sites))
  ACC_local_parafac_self_obs <- rep(NA, length(sites))
  names(CSMF_local_parafac_self_obs) <- sites
  names(ACC_local_parafac_self_obs) <- sites

  # Step 1: Test site self fit
  sim_data_list[[test_site]] <- get_sim_data_from_phmrc_by_location(site = test_site)
  sim_data_filtered_list[[test_site]] <- filter_sparse_causes(sim_data_list[[test_site]], threshold = 0)

  missing_indices <- generate_missing_Yt(sim_data_filtered_list[[test_site]]$filtered_data$data$Y.t, miss_prop)
  sim_data_filtered_list[[test_site]]$filtered_data$data$Y.t[missing_indices] <- NA
  obs_data_test <- get_observed_sim_data(sim_data_filtered_list[[test_site]]$filtered_data)

  train_X <- as.matrix(obs_data_test$obs_data$data$X)
  train_Y <- obs_data_test$obs_data$data.truth$Y.t

  LCVA_local_model_test_obs_fit[[test_site]] <- LCVA.train(
    X = train_X, Y = train_Y,
    Domain = rep(1, length(train_Y)),
    K = K, model = "S",
    Nitr = 2000, thin = 2, seed = 12345, verbose = FALSE
  )

  test_X <- as.matrix(sim_data_filtered_list[[test_site]]$filtered_data$data$X[missing_indices, , drop = FALSE])

  out <- LCVA.pred(
    fit = LCVA_local_model_test_obs_fit[[test_site]],
    X_test = test_X, model = "C",
    Burn_in = 500, Nitr = 4000, verbose = FALSE
  )

  true_test_labels <- sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t[missing_indices]
  Y_pred <- get_assignment(out$Y_test, Burn_in = 2000)
  ACC_local_parafac_self_obs[test_site] <- mean(Y_pred == true_test_labels)

  conf_matrix <- caret::confusionMatrix(
    factor(Y_pred, levels = levels(as.factor(true_test_labels))),
    factor(true_test_labels, levels = levels(as.factor(true_test_labels)))
  )
  confusion_local_parafac_self_obs[[test_site]] <- conf_matrix

  posterior_prev_Y_test <- apply(out$pi_test[-c(1:2000), , drop = FALSE], 2, mean)
  full_truth <- sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t

  aligned_prev <- align_prevalence(
    posterior_prev_Y_test,
    sim_data_filtered_list[[test_site]]$mapping_new_to_origin,
    as.numeric(table(full_truth)) / length(full_truth),
    sim_data_filtered_list[[test_site]]$mapping_new_to_origin[as.numeric(names(table(full_truth)))]
  )

  CSMF_local_parafac_self_obs[test_site] <- CSMF_acc(aligned_prev$aligned_A, aligned_prev$aligned_B)

  cause_level <- sim_data_filtered_list[[test_site]]$filtered_data$data$cause_level
  adjust_est_prev[[test_site]] <- aligned_prev$aligned_A
  cause_ids <- as.integer(names(aligned_prev$aligned_A))
  names(adjust_est_prev[[test_site]]) <- cause_level[cause_ids]

  cause_pred[[test_site]] <- cause_level[
    as.numeric(sim_data_filtered_list[[test_site]]$mapping_new_to_origin[Y_pred])
  ]

  # Step 2: Fit across all other sites
  for (train_site in setdiff(sites, test_site)) {

    sim_data_list[[train_site]] <- get_sim_data_from_phmrc_by_location(site = train_site)
    sim_data_filtered_list[[train_site]] <- filter_sparse_causes(sim_data_list[[train_site]], threshold = 0)

    train_X <- as.matrix(sim_data_filtered_list[[train_site]]$filtered_data$data$X)
    train_Y <- sim_data_filtered_list[[train_site]]$filtered_data$data.truth$Y.t

    LCVA_local_model_fit[[train_site]] <- LCVA.train(
      X = train_X, Y = train_Y,
      Domain = rep(1, length(train_Y)),
      K = K, model = "S",
      Nitr = 2000, thin = 2, seed = 12345, verbose = FALSE
    )

    test_X <- as.matrix(sim_data_filtered_list[[test_site]]$filtered_data$data$X)

    out_full <- LCVA.pred(
      fit = LCVA_local_model_fit[[train_site]],
      X_test = test_X, model = "C",
      Burn_in = 500, Nitr = 4000, return_likelihood = TRUE, verbose = FALSE
    )

    posterior_phi_full[[train_site]] <- apply(out_full$x_given_y_prob[-c(1:2000), , , drop = FALSE], c(2, 3), mean)

    ############## Local prediction ################
    test_site_data <- sim_data_filtered_list[[test_site]]
    train_site_data <- sim_data_filtered_list[[train_site]]

    Y_test_partial <- as.numeric(test_site_data$mapping_new_to_origin[test_site_data$filtered_data$data$Y.t])
    Y_test_partial <- as.numeric(train_site_data$mapping_origin_to_new[Y_test_partial])

    test_X <- test_site_data$filtered_data$data$X
    test_X <- as.matrix(test_X)
    colnames(test_X) <- NULL
    rownames(test_X) <- NULL

    out <- LCVA.pred(
      fit = LCVA_local_model_fit[[train_site]],
      X_test = test_X, Y_test = Y_test_partial, model = "C",
      Burn_in = 500, Nitr = 4000, verbose = FALSE
    )

    # Compute ACC
    true_test_labels <- test_site_data$mapping_new_to_origin[
      test_site_data$filtered_data$data.truth$Y.t[missing_indices]
    ]

    Y_pred <- get_assignment(out$Y_test, Burn_in = 2000)
    Y_pred <- as.numeric(train_site_data$mapping_new_to_origin[Y_pred])

    ACC_local_parafac_self_obs[train_site] <- mean(Y_pred[missing_indices] == true_test_labels)

    # get confusion matrix
    Y_pred_factor <- factor(Y_pred[missing_indices], levels = levels(as.factor(true_test_labels)))
    true_test_labels <- factor(true_test_labels, levels = levels(as.factor(true_test_labels)))

    conf_matrix <- caret::confusionMatrix(Y_pred_factor, true_test_labels)
    confusion_local_parafac_self_obs[[train_site]] <- conf_matrix

    # Compute true test prevalence (for whole data)
    true_test_labels <- test_site_data$filtered_data$data.truth$Y.t
    true_test_label_distribution <- as.numeric(table(true_test_labels)) / length(true_test_labels)

    CSMF.est <- apply(out$pi_test[-c(1:2000), , drop = FALSE], 2, mean)

    # Align prevalence using the correct mapping
    aligned_prev <- align_prevalence(
      CSMF.est,
      train_site_data$mapping_new_to_origin,
      true_test_label_distribution,
      test_site_data$mapping_new_to_origin[as.numeric(names(table(true_test_labels)))]
    )

    # Compute CSMF accuracy for test site
    CSMF_local_parafac_self_obs[train_site] <- CSMF_acc(aligned_prev$aligned_A, aligned_prev$aligned_B)

    cause_level <- train_site_data$filtered_data$data$cause_level
    adjust_est_prev[[train_site]] <- aligned_prev$aligned_A
    cause_ids <- as.integer(names(aligned_prev$aligned_A))
    names(adjust_est_prev[[train_site]]) <- cause_level[cause_ids]

    cause_pred[[train_site]] <- cause_level[Y_pred]
  }

  list(
    csmf_local = adjust_est_prev,
    cause_pred_local = cause_pred,
    csmf_acc_local = CSMF_local_parafac_self_obs,
    acc_local = ACC_local_parafac_self_obs,
    confusion_local = confusion_local_parafac_self_obs,
    posterior_phi_full = posterior_phi_full,
    LCVA_local_model_test_obs_fit = LCVA_local_model_test_obs_fit,
    sim_data_filtered_list = sim_data_filtered_list
  )
}
