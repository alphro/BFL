# BFL (Bayesian Federated Learning for Verbal Autopsy)

This repository contains a **prototype implementation (v0.01)** of Bayesian
Federated Learning (BFL) for verbal autopsy data.

The package is structured to clearly separate:
- user-facing APIs,
- core BFL mechanics,
- data alignment / validation,
- and legacy reference code.

---

## Package structure

```mermaid
mindmap
  root((BFL))
    DESCRIPTION
    NAMESPACE
    R
      api
        prep_BFL_local_LCVA.R
        run_BFL_global.R
        predict_BFL_local.R
      core
        bfl_prediction.R
        lcva_fit.R
        posterior_extraction.R
        bfl_stan_runner.R
      data_prep
        align_causes.R
        build_stan_data.R
        validate_inputs.R
      utils
        argmax.R
        safe_log.R
        softmax.R
      legacy
        assist_functions.R
        balance_base.R
        balance_local_fit.R
        run_BFL_model.R
    inst
      stan
        no_partial_labels.stan
      extdata
        phmrc_clean.csv
    vignettes
      01-bfl-data-format.Rmd
      02-phmrc-to-bfl.Rmd
    tests
      testthat
        test-toy-pipeline.R
```
