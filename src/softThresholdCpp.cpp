#include <RcppArmadillo.h>
using namespace Rcpp;

// Proximal operator for \lambda |x|_1 (from Vu's package)
struct SoftThresholdOp
{
  SoftThresholdOp(const double& z) : z(z) {}
  const double operator() (const double& x) const {
    return ((x > 0) - (x < 0)) * std::max(0.0, std::abs(x) - z);
  }

private:
    const double z;
};

// [[Rcpp::export]]
arma::mat softThresholdCpp(arma::mat& x, double& lambda) {
  x.transform( SoftThresholdOp(lambda) );
  return x;
}

