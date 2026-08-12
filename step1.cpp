// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List sample_Y_tilde(const mat& X, const mat& b, const vec& cC_inv_xi,
                    IntegerVector Z, NumericVector sej, NumericVector spj,
                    List ind_list, List pool_list,
                    IntegerVector Y_tilde, int t, int N,
                    IntegerVector success_counts) {
  
  if (t < 2) throw std::invalid_argument("t must be >= 2");
  int success_count = 0;
  IntegerVector output_Y_tilde = clone(Y_tilde); // Clone to avoid modifying input
  
  for (int i = 0; i < N; i++) {
    IntegerVector ind_j = ind_list[i]; // Pool indices (1-based) for individual i+1
    if (ind_j.size() == 0) continue; // Skip if no pools
    
    // Linear predictor: X[i,] %*% b[t-1,] + cC_inv_xi[i]
    double dot_product = dot(X.row(i), b.row(t - 2).t());
    double linear_predictor = dot_product + cC_inv_xi[i];
    double eta = 1.0 / (1.0 + std::exp(-linear_predictor));
    
    // p1: eta * product over all pools in ind_j
    double prod_p1 = 1.0;
    for (int k = 0; k < ind_j.size(); k++) {
      int idx = ind_j[k] - 1; // 0-based pool index
      prod_p1 *= std::pow(sej[idx], Z[idx]) * std::pow(1.0 - sej[idx], 1 - Z[idx]);
    }
    double p1_val = eta * prod_p1;
    
    // p0 computation
    IntegerVector pool_members = pool_list[ind_j[0] - 1];
    int sij = 0;
    for (int m = 0; m < pool_members.size(); m++) {
      int pm = pool_members[m] - 1; // 0-based individual index
      sij += output_Y_tilde[pm];
    }
    sij -= output_Y_tilde[i];
    
    double prod_p0 = 1.0;
    for (int k = 0; k < ind_j.size(); k++) {
      int idx = ind_j[k] - 1;
      double sens_term = std::pow(sej[idx], Z[idx]) * std::pow(1.0 - sej[idx], 1 - Z[idx]);
      double spec_term = std::pow(1.0 - spj[idx], Z[idx]) * std::pow(spj[idx], 1 - Z[idx]);
      int sij_k = (k == 0) ? sij : 0;
      prod_p0 *= std::pow(sens_term, static_cast<int>(sij_k > 0)) * std::pow(spec_term, static_cast<int>(!(sij_k > 0)));
    }
    double p0_val = (1.0 - eta) * prod_p0;
    
    double total = p1_val + p0_val;
    double prob;
    if (total > 0.0) {
      prob = p1_val / total;
      success_count++;
    } else {
      prob = 0.5;
    }
    output_Y_tilde[i] = R::rbinom(1, prob);
  }
  
  // Return a list with named elements
  return List::create(_["Y_tilde"] = output_Y_tilde,
                      _["success_count"] = success_count);
}