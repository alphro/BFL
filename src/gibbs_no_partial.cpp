#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Gibbs sampler for BFL no-partial-labels model
//
// Two prior options for lambda_k (model weights per cause k):
//
//   logistic_normal = FALSE  [default]
//     lam_k ~ Dirichlet(1,...,1)
//     Full conditional is conjugate: lam_k | Z,W ~ Dir(1 + v_{k,m})
//
//   logistic_normal = TRUE
//     beta_raw[k,m] ~ N(0,1),  lam_k = softmax(beta_eff[k])
//     where beta_eff[k,m] = beta_raw[k,m]  if model m is active for cause k
//                         = -10            otherwise  (mirrors Stan code)
//     Full conditional for active betas is NOT conjugate:
//       log p(beta_k | Z, W) proportional to
//         -0.5 * ||beta_active||^2
//         + sum_{m active} v_{k,m} * (beta_m - log_sum_exp(beta_eff))
//     A Metropolis-Hastings random-walk step is used for this block.
//     All other steps (Z, W, pi) remain conjugate/exact.
//
// ---------------------------------------------------------------------------

// ---- Utilities -------------------------------------------------------------

static NumericVector rdirichlet_cpp(const NumericVector& alpha) {
  int K = alpha.size();
  NumericVector x(K);
  double s = 0.0;
  for (int k = 0; k < K; k++) {
    x[k] = (alpha[k] > 0.0) ? R::rgamma(alpha[k], 1.0) : 0.0;
    s += x[k];
  }
  if (s <= 0.0) { std::fill(x.begin(), x.end(), 1.0 / K); return x; }
  for (int k = 0; k < K; k++) x[k] /= s;
  return x;
}

// Categorical sample from unnormalized probabilities; uniform fallback.
static int rcategorical_cpp(const NumericVector& p) {
  int n = p.size();
  double s = 0.0;
  for (int k = 0; k < n; k++) s += p[k];
  if (s <= 0.0) return (int)(R::unif_rand() * n);
  double u = R::unif_rand() * s;
  double cum = 0.0;
  for (int k = 0; k < n - 1; k++) {
    cum += p[k];
    if (u < cum) return k;
  }
  return n - 1;
}

// Numerically stable log-sum-exp over a std::vector.
static double log_sum_exp(const std::vector<double>& x) {
  double mx = *std::max_element(x.begin(), x.end());
  double s = 0.0;
  for (double v : x) s += std::exp(v - mx);
  return mx + std::log(s);
}

// Softmax of beta_eff -> lambda, storing into the provided output vector.
static void softmax_into(const std::vector<double>& beta_eff,
                         NumericVector& lam, int M) {
  double lse = log_sum_exp(beta_eff);
  for (int m = 0; m < M; m++)
    lam[m] = std::exp(beta_eff[m] - lse);
}

// Log-posterior for beta_active[k] given count vector v_k (length = |active|).
//
//   log p ∝ -0.5 * ||beta_active||^2
//           + sum_{j} v_k[j] * (beta_active[j] - log_sum_exp(beta_eff))
//
// beta_eff is the full M-length vector (inactive entries = -10).
static double log_post_beta(const std::vector<double>& beta_active,
                             const std::vector<int>&    active_ms,
                             const std::vector<double>& beta_eff_full,
                             const std::vector<double>& v_k) {
  double lse = log_sum_exp(beta_eff_full);
  double lp  = 0.0;
  for (int j = 0; j < (int)active_ms.size(); j++) {
    lp += -0.5 * beta_active[j] * beta_active[j];      // N(0,1) prior
    lp += v_k[j] * (beta_active[j] - lse);             // multinomial likelihood
  }
  return lp;
}

// ---- Main exported function ------------------------------------------------

//' Gibbs sampler for BFL no-partial-labels model
//'
//' Run one MCMC chain. Call from R once per chain with a different seed.
//'
//' @param phi            Flattened phi array (column-major, dim N x C_max x M).
//' @param model_presence num_causes x M integer matrix; 1 = model active.
//' @param causes_r       1-indexed cause codes, length num_causes.
//' @param N,C_max,M,num_causes Array dimensions.
//' @param n_iter         Total iterations (warmup + sampling).
//' @param n_warmup       Iterations to discard.
//' @param logistic_normal If TRUE, use N(0,1)->softmax prior with MH step.
//' @param mh_scale       Random-walk MH step-size (only used when logistic_normal=TRUE).
//'
//' @return List: pi [n_keep x num_causes], lambda [n_keep*num_causes*M],
//'   mh_accept_rate [num_causes] (NA when logistic_normal=FALSE).
//' @keywords internal
// [[Rcpp::export]]
List gibbs_no_partial_cpp(
    const NumericVector& phi,
    const IntegerMatrix& model_presence,
    const IntegerVector& causes_r,
    int N, int C_max, int M, int num_causes,
    int n_iter, int n_warmup,
    bool logistic_normal = false,
    double mh_scale      = 0.5
) {
  int n_keep = n_iter - n_warmup;

  // Active-model index and reverse lookup (all 0-based)
  std::vector<std::vector<int>> active(num_causes);
  std::vector<std::vector<int>> m_to_j(num_causes, std::vector<int>(M, -1));
  for (int k = 0; k < num_causes; k++)
    for (int m = 0; m < M; m++)
      if (model_presence(k, m) == 1) {
        m_to_j[k][m] = (int)active[k].size();
        active[k].push_back(m);
      }

  std::vector<int> causes(num_causes);
  for (int k = 0; k < num_causes; k++) causes[k] = causes_r[k] - 1;

  auto phi_at = [&](int i, int c, int m) -> double {
    return phi[i + c * N + m * N * C_max];
  };

  // ---- Initialize ----
  NumericVector pi(num_causes, 1.0 / num_causes);

  // lambda[k][m]; inactive m always 0 (Dirichlet path) or ~0 (LN path)
  std::vector<NumericVector> lambda(num_causes);
  for (int k = 0; k < num_causes; k++) {
    lambda[k] = NumericVector(M, 0.0);
    int nact = (int)active[k].size();
    if (nact > 0) {
      double v = 1.0 / nact;
      for (int m : active[k]) lambda[k][m] = v;
    }
  }

  // beta_raw per cause (only used when logistic_normal = true)
  // beta_eff[k][m] = beta_raw for active m, -10 for inactive
  std::vector<std::vector<double>> beta_eff(num_causes,
                                            std::vector<double>(M, -10.0));
  if (logistic_normal) {
    for (int k = 0; k < num_causes; k++)
      for (int m : active[k])
        beta_eff[k][m] = 0.0;  // start at 0 -> uniform lambda
  }

  IntegerVector Z(N, 0), W(N, 0);

  // ---- Output ----
  NumericMatrix pi_out(n_keep, num_causes);
  NumericVector lambda_out((R_xlen_t)n_keep * num_causes * M, 0.0);

  // MH diagnostics
  std::vector<int> mh_accept(num_causes, 0);
  std::vector<int> mh_total(num_causes, 0);

  NumericVector prob_z(num_causes);

  // ---- Main Gibbs loop ----
  for (int iter = 0; iter < n_iter; iter++) {

    // Step 1 & 2: sample Z[i] then W[i]
    for (int i = 0; i < N; i++) {
      for (int k = 0; k < num_causes; k++) {
        int c = causes[k];
        double s = 0.0;
        for (int m : active[k])
          s += phi_at(i, c, m) * lambda[k][m];
        prob_z[k] = pi[k] * s;
      }
      Z[i] = rcategorical_cpp(prob_z);

      int k = Z[i], c = causes[k];
      int nact = (int)active[k].size();
      NumericVector prob_w(nact);
      for (int j = 0; j < nact; j++) {
        int m = active[k][j];
        prob_w[j] = phi_at(i, c, m) * lambda[k][m];
      }
      W[i] = active[k][rcategorical_cpp(prob_w)];
    }

    // Step 3: sample pi | Z ~ Dir(1 + n_k)
    NumericVector alpha_pi(num_causes, 1.0);
    for (int i = 0; i < N; i++) alpha_pi[Z[i]] += 1.0;
    pi = rdirichlet_cpp(alpha_pi);

    // Step 4: update lambda[k]
    for (int k = 0; k < num_causes; k++) {
      int nact = (int)active[k].size();
      if (nact == 0) continue;

      // Count model assignments within cause k: v_k[j] = #{i: Z_i=k, W_i=active[k][j]}
      std::vector<double> v_k(nact, 0.0);
      for (int i = 0; i < N; i++)
        if (Z[i] == k) {
          int j = m_to_j[k][W[i]];
          if (j >= 0) v_k[j] += 1.0;
        }

      if (!logistic_normal) {
        // ---- Conjugate Dirichlet update ----
        NumericVector alpha_lam(nact);
        for (int j = 0; j < nact; j++) alpha_lam[j] = 1.0 + v_k[j];
        NumericVector lam_new = rdirichlet_cpp(alpha_lam);
        for (int m = 0; m < M; m++) lambda[k][m] = 0.0;
        for (int j = 0; j < nact; j++)
          lambda[k][active[k][j]] = lam_new[j];

      } else {
        // ---- MH random-walk step for beta_raw (active components only) ----
        //
        // Propose: beta_prop[active] = beta_curr[active] + mh_scale * N(0, I)
        // Inactive entries stay at -10; proposal is symmetric so Hastings ratio = 1.
        //
        // log p(beta | v_k) = -0.5 * ||beta_active||^2
        //                   + sum_j v_k[j] * (beta_active[j] - log_sum_exp(beta_eff))

        std::vector<double> beta_curr_active(nact);
        for (int j = 0; j < nact; j++)
          beta_curr_active[j] = beta_eff[k][active[k][j]];

        double lp_curr = log_post_beta(beta_curr_active, active[k],
                                       beta_eff[k], v_k);

        // Draw proposal
        std::vector<double> beta_prop_active(nact);
        std::vector<double> beta_eff_prop = beta_eff[k];  // copy; inactive stay -10
        for (int j = 0; j < nact; j++) {
          beta_prop_active[j] = beta_curr_active[j] + mh_scale * R::norm_rand();
          beta_eff_prop[active[k][j]] = beta_prop_active[j];
        }

        double lp_prop = log_post_beta(beta_prop_active, active[k],
                                       beta_eff_prop, v_k);

        mh_total[k]++;
        if (std::log(R::unif_rand()) < lp_prop - lp_curr) {
          // Accept
          beta_eff[k] = beta_eff_prop;
          mh_accept[k]++;
        }

        // Update lambda from (possibly new) beta_eff
        softmax_into(beta_eff[k], lambda[k], M);
      }
    }

    // Store post-warmup draws
    if (iter >= n_warmup) {
      int s = iter - n_warmup;
      for (int k = 0; k < num_causes; k++)
        pi_out(s, k) = pi[k];
      for (int k = 0; k < num_causes; k++)
        for (int m = 0; m < M; m++)
          lambda_out[s + (long)k * n_keep + (long)m * n_keep * num_causes] =
              lambda[k][m];
    }
  }

  // MH acceptance rates (NA when Dirichlet path used)
  NumericVector mh_accept_rate(num_causes, NA_REAL);
  if (logistic_normal)
    for (int k = 0; k < num_causes; k++)
      mh_accept_rate[k] = mh_total[k] > 0
          ? (double)mh_accept[k] / mh_total[k] : NA_REAL;

  return List::create(
    Named("pi")             = pi_out,
    Named("lambda")         = lambda_out,
    Named("mh_accept_rate") = mh_accept_rate
  );
}
