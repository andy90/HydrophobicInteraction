#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector structure_factorC(NumericMatrix corrs) {
  int ncol = corrs.ncol();
  int nrow = corrs.nrow();
  
  int nrow_pair = nrow*(nrow-1)/2;
  NumericMatrix dis_corrs(nrow_pair, ncol);
  
  int ij = 0;
  for (int i=0; i < (nrow-1); i++){
    for (int j = (i+1); j< nrow; j++){
      for (int k=0; k < ncol; k++){
        dis_corrs(ij, k) = corrs(i, k) - corrs(j, k);
      }
      ij = ij + 1;
    }
  }
  
  double k0 = 2*M_PI;
  int nk = 50;
  NumericVector sk(50);
  for (int ik = 0; ik < nk; ik++){
    for (int i = 0; i < nrow_pair; i++){
      for (int j = 0; j < ncol; j++){
        sk(ik) += 2 * cos((ik+1) * k0 * dis_corrs(i, j)) /3;
      }
    }
    sk(ik) += nrow;
  }
  
  return sk;
}


// [[Rcpp::export]]
NumericVector pp_longitudinalC(NumericMatrix corrs, NumericMatrix dipoles) {
  int ncol = corrs.ncol();
  int nrow = corrs.nrow();
  
  int nrow_pair = nrow*(nrow-1)/2;
  NumericMatrix dis_corrs(nrow_pair, ncol);
  NumericMatrix mult_dipoles(nrow_pair, ncol);
  
  int ij = 0;
  for (int i=0; i < (nrow-1); i++){
    for (int j = (i+1); j< nrow; j++){
      for (int k=0; k < ncol; k++){
        dis_corrs(ij, k) = corrs(i, k) - corrs(j, k);
        mult_dipoles(ij, k) = dipoles(i, k) * dipoles(j, k);
      }
      ij = ij + 1;
    }
  }
  
  double total_dipolesquare  = 0;
  for (int i=0; i<nrow; i++){
    for (int j=0; j<ncol; j++){
      total_dipolesquare += dipoles(i, j) * dipoles(i, j); 
    }
  }
  
  double k0 = 2*M_PI;
  int nk = 50;
  NumericVector sk(50);
  for (int ik = 0; ik < nk; ik++){
    for (int i = 0; i < nrow_pair; i++){
      for (int j = 0; j < ncol; j++){
        sk(ik) += 2 * cos((ik+1) * k0 * dis_corrs(i, j)) * mult_dipoles(i, j) /3;
      }
    }
    sk(ik) += total_dipolesquare/3;
  }
  
  return sk;
}


/*** R
*/
