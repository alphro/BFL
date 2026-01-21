data {
  int<lower=1> N; // Number of observations
  int<lower=1> C_max; // Maximum index of causes of death
  int<lower=1> M; // Number of models
  int<lower=1> num_causes; // Actual number of causes (4 in your example)
  int<lower=1, upper=C_max> causes[num_causes]; // Actual causes list (1, 2, 3, 5)
  int<lower=0, upper=1> model_presence[num_causes, M]; // Models active for each cause, 0 or 1
  int<lower=0> count[num_causes]; // Number of models active for each cause
  real<lower=0, upper=1> phi[N, C_max, M]; // Probabilities for each cause and model
}

parameters {
  simplex[num_causes] pi; // Prior probabilities for Y categories
  vector[M] beta_raw[num_causes]; // Raw beta parameters for each cause and model
}

transformed parameters {
  simplex[M] lambda[num_causes]; // Softmax-transformed lambda parameters for each cause

  for (c in 1:num_causes) {
    vector[M] beta_effective;
    for (m in 1:M) {
      // Apply the indicator of model presence directly in the beta calculation
      // beta_effective[m] = beta_raw[c, m] * model_presence[c, m];
      if(model_presence[c, m] == 1){
          beta_effective[m] = beta_raw[c, m];
      }else{
        beta_effective[m] = -10;
      }
    }
    lambda[c] = softmax(beta_effective); // Apply softmax to ensure it's a simplex
  }
}

model {
  for (c in 1:num_causes) {
    for (m in 1:M) {
        beta_raw[c, m] ~ normal(0, 1);
    }
  }

  pi ~ dirichlet(rep_vector(1.0, num_causes)); // Uniform prior over actual causes

  for (i in 1:N) {
    real total_prob_i = 0;

    for (idx in 1:num_causes) {
      int c = causes[idx]; // Actual cause
      //real total_prob_i = 0;
      real inner_sum = 0;

      for (m in 1:M) {
        if (model_presence[idx, m] > 0) {
          inner_sum += phi[i, c, m] * lambda[idx][m];
          // total_prob_i += pi[idx] * phi[i, c, m] * lambda[idx][m];
        }
      }

       total_prob_i += pi[idx] * inner_sum;
      // target += log(total_prob_i); // Log-likelihood addition
    }

        // Only take the log and add to target if total_prob_i is positive
    if (total_prob_i > 0) {
      target += log(total_prob_i); // Add the log of this probability to the target (log-likelihood)
    } else {
      print("Warning: Non-positive total probability encountered.");
      target += negative_infinity(); // Handle log(0) case appropriately
    }
  }
}

