library(BFL)
library(rstan)
library(caret)
install.packages("remotes")
remotes::install_github("richardli/LCVA")
library(LCVA)

test_site <- "AP"  # example from the PDF
sites <- c("Mexico", "AP", "Bohol", "Dar", "Pemba", "UP")
K <- 5
miss_prop <- 0.8

## Local fit (balance) ----
local_res <- execute_balance_local_fit(
  test_site = test_site,
  sites = sites,
  K = K,
  miss_prop = miss_prop
)

# BFL base
base_res <- execute_balance_base(
  test_site = test_site,
  sites = sites,
  posterior_phi_full = local_res$posterior_phi_full,
  sim_data_filtered_list = local_res$sim_data_filtered_list
)

print(base_rec$csmf_acc)
print(base_rec$acc)
print(base_res$conf_matrix)
