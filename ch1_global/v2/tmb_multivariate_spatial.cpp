#include <TMB.hpp>

// Multivariate Gaussian trait model with:
// - shared environmental fixed effects for all traits
// - trait-specific species random intercepts
// - low-rank shared spatial factors on an SPDE mesh
//
template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(y);
  DATA_IVECTOR(obs_index);
  DATA_IVECTOR(trait_index);
  DATA_IVECTOR(species_index);
  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(A);
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  DATA_INTEGER(n_traits);
  DATA_INTEGER(n_factors);

  PARAMETER_MATRIX(beta);               // predictors x traits
  PARAMETER_MATRIX(species_re);         // species x traits
  PARAMETER_VECTOR(log_sigma_species);  // trait-specific species SD
  PARAMETER_VECTOR(log_sigma_obs);      // trait-specific residual SD
  PARAMETER_MATRIX(omega);              // mesh vertices x factors
  PARAMETER_MATRIX(loadings);           // traits x factors; upper triangle mapped fixed at 0 in R
  PARAMETER_VECTOR(log_kappa);          // factor-specific spatial scale
  PARAMETER_VECTOR(log_tau);            // factor-specific spatial precision scale

  Type nll = 0.0;
  const int n_species = species_re.rows();
  const int n_mesh = omega.rows();

  for (int t = 0; t < n_traits; ++t) {
    Type sd = exp(log_sigma_species(t));
    for (int s = 0; s < n_species; ++s) {
      nll -= dnorm(species_re(s, t), Type(0.0), sd, true);
    }
  }

  matrix<Type> spatial_at_obs(A.rows(), n_factors);
  for (int k = 0; k < n_factors; ++k) {
    Type kappa = exp(log_kappa(k));
    Type tau = exp(log_tau(k));
    Eigen::SparseMatrix<Type> Q = pow(tau, Type(2.0)) *
      (pow(kappa, Type(4.0)) * M0 + Type(2.0) * pow(kappa, Type(2.0)) * M1 + M2);
    vector<Type> field(n_mesh);
    for (int j = 0; j < n_mesh; ++j) field(j) = omega(j, k);
    nll += density::GMRF(Q)(field);
    spatial_at_obs.col(k) = A * field;
  }

  vector<Type> pred(y.size());
  for (int q = 0; q < y.size(); ++q) {
    int i = obs_index(q);
    int t = trait_index(q);
    Type mu = species_re(species_index(i), t);
    for (int p = 0; p < X.cols(); ++p) mu += X(i, p) * beta(p, t);
    for (int k = 0; k < n_factors; ++k) mu += spatial_at_obs(i, k) * loadings(t, k);
    Type sd = exp(log_sigma_obs(t));
    pred(q) = mu;
    nll -= dnorm(y(q), mu, sd, true);
  }

  ADREPORT(beta);
  ADREPORT(loadings);
  ADREPORT(log_sigma_species);
  ADREPORT(log_sigma_obs);
  ADREPORT(log_kappa);
  ADREPORT(log_tau);
  REPORT(pred);
  REPORT(spatial_at_obs);
  return nll;
}
