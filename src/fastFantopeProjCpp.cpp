#include <RcppArmadillo.h>
using namespace Rcpp;

// Note: to use eigs_sym, need to
// edit include/armadillo_bits/config.hpp and uncomment ARMA_USE_ARPACK
// which is in /Library/Frameworks/R.framework/Versions/3.2/Resources/library/RcppArmadillo

int simplexCpp(arma::vec& eigval, int ndim){
  int rank = 3;

  return rank;
}

// [[Rcpp::export]]
arma::sp_mat fastFantopeProjCpp(arma::sp_mat& S, int& ndim) {
  int rank;
  arma::vec eigval;
  arma::vec eigval_;
  arma::mat eigvec;

  // eigen values
  arma::eigs_sym(eigval,S, S.n_rows);
  // arma::eig_sym(eigval,S);

  // threshold theta
  rank = simplexCpp(eigval, ndim);

  // eigen vectors
  arma::eigs_sym(eigval_, eigvec, S, rank);
  // arma::eig_sym(eigval_, eigvec, S);

  // reconstruct
  S = (
    eigvec *
    diagmat(eigval.subvec(eigval.n_elem - rank, eigval.n_elem - 1)) *
    eigvec.t()
  );

  return S;
}



