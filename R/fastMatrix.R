#' A function that reconstructs matrix using eigen values and vectors
#' @param Din diagonal matrix
#' @param Vin the eigenvector matrix
#' @return recon, the reconstructed matrix Vmat %*% Dmat %*% t(Vmat)

# reconSVDcpp <- cxxfunction(signature(Din="numeric", Vin="numeric"),
#                              plugin="RcppArmadillo",
#                              body = '
#         mat D = as<mat>(Din);
#         mat V = as<mat>(Vin);
#         mat recon = V * D * trans(V);
#         return wrap(recon);
#         ',
#                              includes='
#         using namespace Rcpp;
#         using namespace arma;
#         ')

reconSVDcpp <- cppFunction(depends="RcppArmadillo",
        code = '
        arma::mat reconSVD(arma::mat & D, arma::mat & V){
        arma::mat recon = V * D * trans(V);
        return recon;
        }
        ')


