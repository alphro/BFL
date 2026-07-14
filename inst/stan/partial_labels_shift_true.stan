data {
  int<lower=1> N; // Number of observations
  int<lower=1> C_max; // Maximum index of causes of death
  int<lower=1> M; // Number of models
  int<lower=1> num_causes; // Actual number of causes (4 in your example)
  int<lower=1, upper=C_max> causes[num_causes]; // Actual causes list (1, 2, 3, 5)
  int<lower=0, upper=1> model_presence[num_causes, M]; // Models active for each cause, 0 or 1
  int<lower=0> count[num_causes]; // Number of models active for each cause
  real<lower=0, upper=1> phi[N, C_max, M]; // Probabilities for each cause and model
  int<lower=0, upper=1> Y_known[N];  // Indicator: 1 if Y_i is known, 0 otherwise

  // MIN FIX: use cause INDEX (1..num_causes) for labeled Y
  int<lower=1, upper=num_causes> Y_idx[N];  // Known category labels, only used where Y_known[i] == 1
}


parameters {
  simplex[num_causes] pi_O; // Prior probabilities for Y categories
  simplex[num_causes] pi;
  vector[M] beta_raw[num_causes]; // Raw beta parameters for each cause and model
  // vector[num_causes] alpha_raw; 
  // vector[num_causes] epsilon;
}

transformed parameters {
  simplex[M] lambda[num_causes]; // Softmax-transformed lambda parameters for each cause
  // simplex[num_causes] pi_O;
  // simplex[num_causes] pi_U;
  
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
  
   // pi_O = softmax(alpha_raw);
   //  pi_U = softmax(alpha_raw + epsilon);
}

model {
  // Prior distributions for pi and lambda
  pi ~ dirichlet(rep_vector(1.0, num_causes)); // Uniform prior over actual causes
  pi_O ~ dirichlet(rep_vector(1.0, num_causes)); // Uniform prior over actual causes
  for (c in 1:num_causes) {
    for (m in 1:M) {
        beta_raw[c, m] ~ normal(0, 1);
    }
    // alpha_raw[c] ~ normal(0, 1);
    // epsilon[c] ~ normal(0, 0.1);
  }
 

  for (i in 1:N) {
    if (Y_known[i] == 1) {
      // Use the known Y_i directly

      // MIN FIX: c1 is now a cause INDEX (1..num_causes), not a cause code
      int c1 = Y_idx[i];

      // MIN FIX: map cause index -> actual cause code for phi indexing
      int c_code = causes[c1];

      // int idx = which_first(to_array_1d(causes) == c);
      real inner_sum = 0;
      for (m in 1:M) {
        if (model_presence[c1, m] > 0) {
              inner_sum += phi[i, c_code, m] * lambda[c1][m];
        }
      }
      // Log of this particular configuration, because Y_i is known
      target += log(pi_O[c1] * inner_sum);
    } else {
      // Y_i is unknown, proceed as before
      real total_prob_i = 0;
      
      for (c in 1:num_causes) {
        // MIN FIX: map cause index -> actual cause code for phi indexing
        int c_code = causes[c];

        // int c = causes[idx]; // Actual cause
        //real total_prob_i = 0;
        real inner_sum = 0;
        
        for (m in 1:M) {
          if (model_presence[c, m] > 0) {
            inner_sum += phi[i, c_code, m] * lambda[c][m];
            // total_prob_i += pi[idx] * phi[i, c, m] * lambda[idx][m];
          }
        }
        
         total_prob_i += pi[c] * inner_sum;
        // target += log(total_prob_i); // Log-likelihood addition
      }
       target += log(total_prob_i);
          // Only take the log and add to target if total_prob_i is positive
     // if (total_prob_i > 0) {
      //  target += log(total_prob_i); // Add the log of this probability to the target (log-likelihood)
      //} else {
      //  print("Warning: Non-positive total probability encountered.");
      //  target += negative_infinity(); // Handle log(0) case appropriately
      //}
    }
  }
}