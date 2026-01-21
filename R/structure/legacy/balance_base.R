#' Execute BFL base model (balanced)
#'
#' @param test_site character
#' @param sites character vector
#' @param posterior_phi_full list
#' @param sim_data_filtered_list list
#' @return list(csmf, cause_pred, acc, conf_matrix, csmf_acc)
#' @export
execute_balance_base <- function(
    test_site,
    sites,
    posterior_phi_full,
    sim_data_filtered_list
) {

  # path to Stan program inside the installed package
  stan_path <- system.file("stan", "no_partial_labels.stan", package = "BFL")
  if (stan_path == "") stop("Stan file not found in package: inst/stan/no_partial_labels.stan")

  # compile Stan model
  stan_model <- rstan::stan_model(file = stan_path)

  cause_list <- list()
  for (i in seq_along(sites)) {
    train_site <- sites[i]
    cause_list[[train_site]] <- as.numeric(sim_data_filtered_list[[train_site]]$mapping_new_to_origin)
  }

  train_sites <- setdiff(sites, test_site)

  cause_lists_train <- lapply(train_sites, function(site) cause_list[[site]])
  C_max <- max(unlist(cause_lists_train))
  all_causes <- sort(unique(unlist(cause_lists_train)))
  num_causes <- length(all_causes)

  models <- seq_along(train_sites)
  site_to_model <- setNames(models, train_sites) # kept for parity (unused, like original)

  model_presence <- matrix(
    0,
    nrow = num_causes,
    ncol = length(train_sites),
    dimnames = list(all_causes, train_sites)
  )

  for (site in train_sites) {
    for (cause in cause_list[[site]]) {
      model_presence[as.character(cause), site] <- 1
    }
  }

  count <- apply(model_presence, 1, function(x) sum(x != 0))

  # kept for parity (unused later, like original)
  missing_indices <- which(is.na(sim_data_filtered_list[[test_site]]$filtered_data$data$Y.t))

  whole_X_i <- as.matrix(sim_data_filtered_list[[test_site]]$filtered_data$data$X)
  colnames(whole_X_i) <- NULL
  rownames(whole_X_i) <- NULL

  num_models <- length(train_sites)

  posterior_input_partial_i <- array(0, dim = c(nrow(whole_X_i), C_max, num_models))
  for (j in seq_along(train_sites)) {
    model_site <- train_sites[j]
    posterior_phi_j <- posterior_phi_full[[model_site]]
    for (i in seq_along(cause_list[[model_site]])) {
      cause_index <- cause_list[[model_site]][i]
      posterior_input_partial_i[, cause_index, j] <- posterior_phi_j[, i]
    }
  }

  stan_data <- list(
    N = nrow(whole_X_i),
    C_max = C_max,
    M = num_models,
    num_causes = num_causes,
    causes = all_causes,
    model_presence = model_presence,
    count = count,
    phi = posterior_input_partial_i
  )

  fit <- rstan::sampling(stan_model, data = stan_data, iter = 2000, chains = 4, refresh = 0)
  posterior_samples <- rstan::extract(fit)

  pred_res <- get_global_post_pred_Y(
    posterior_input_partial_i,
    posterior_samples$lambda,
    posterior_samples$pi
  )

  Y_pred <- pred_res$posterior_pred_Y

  true_test_labels <- as.numeric(sim_data_filtered_list[[test_site]]$mapping_new_to_origin[
    sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t
  ])

  acc <- get_acc(Y_pred, true_test_labels)

  Y_pred_factor <- factor(get_Y_pred(Y_pred), levels = levels(as.factor(true_test_labels)))
  true_test_labels_factor <- factor(true_test_labels, levels = levels(as.factor(true_test_labels)))
  conf_matrix <- caret::confusionMatrix(Y_pred_factor, true_test_labels_factor)

  true_test_labels <- as.numeric(sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t)
  true_test_label_distribution <- as.numeric(table(true_test_labels)) / length(true_test_labels)

  aligned_prev <- align_prevalence(
    colMeans(posterior_samples$pi),
    setNames(all_causes, all_causes),
    true_test_label_distribution,
    sim_data_filtered_list[[test_site]]$mapping_new_to_origin[as.numeric(names(table(true_test_labels)))]
  )

  csmf_acc <- CSMF_acc(aligned_prev$aligned_A, aligned_prev$aligned_B)

  cause_level <- sim_data_filtered_list[[test_site]]$filtered_data$data$cause_level
  adjust_est_prev <- aligned_prev$aligned_A
  names(adjust_est_prev) <- cause_level

  cause_pred <- cause_level[get_Y_pred(Y_pred)]

  list(
    csmf = adjust_est_prev,
    cause_pred = cause_pred,
    acc = acc,
    conf_matrix = conf_matrix,
    csmf_acc = csmf_acc
  )
}
