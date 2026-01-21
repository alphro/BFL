run_BFL_balance_model <- function(model_type,
                                  test_site,
                                  sites,
                                  sim_data_filtered_list,
                                  posterior_phi_full,
                                  LCVA_local_model_test_obs_fit = NULL) {

  # minimal: ignore model_type, always run base
  execute_balance_base(
    test_site = test_site,
    sites = sites,
    posterior_phi_full = posterior_phi_full,
    sim_data_filtered_list = sim_data_filtered_list
  )
}
