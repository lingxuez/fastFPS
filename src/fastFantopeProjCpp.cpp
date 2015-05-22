#include <RcppArmadillo.h>
using namespace Rcpp;

int simplexCpp(arma::vec& eigval, int ndim){
  int rank = 3;

  return rank;
}

// [[Rcpp::export]]
arma::mat fastFantopeProjCpp(arma::mat& S, int& ndim) {
  int rank;
  arma::vec eigval;
  arma::vec eigval_;
  arma::mat eigvec;

  // eigen values
  // arma::eigs_sym(eigval,S, S.n_rows);
  arma::eig_sym(eigval,S);

  // threshold theta
  rank = simplexCpp(eigval, ndim);

  // eigen vectors
  // arma::eigs_sym(eigval_, eigvec, S, rank);
  arma::eig_sym(eigval_, eigvec, S);

  // reconstruct
  S = (
    eigvec *
    diagmat(eigval.subvec(eigval.n_elem - rank, eigval.n_elem - 1)) *
    eigvec.t()
  );

  return S;
}



