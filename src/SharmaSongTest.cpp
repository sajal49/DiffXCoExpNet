// Created by Sajal Kumar
// Copyright (c) NMSU Song lab

#include "methods.h"
#include <armadillo>

// Sharma-Song test implementation
// Given 'k' contingency tables representing 'k' experimental conditions
// the Sharma-Song test returns the amount of second order difference between the tables
ldouble SharmaSongTest(const std::vector<frame<int > > & tables, ldouble & pvalue, double & estimate){

  int k = tables.size(), total_sum=0;
  int temp;
  ldouble stat = 0;
  // Get independent standard normal variables vectors and sample sizes for each table
  vec<int> ssizes(k);
  std::vector<ublas_vec<double > > emat = get_e_matrix(tables, ssizes);
  size_t df=0;

  // check to see if all ssizes are non-zero
  for(int i=0; i<k; i++){
    if(ssizes[i] == 0){
      std::cerr<<"No table can have 0 sample size. Exiting!"<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // pool emat
  ublas_vec<double> pooled_emat(emat[0].size(), 0);
  for(int i=0; i<k; i++){
    for(int j=0; j<emat[0].size(); j++){
      pooled_emat(j) += emat[i](j);
    }
  }

  // scaling coefficients
  ublas_vec<double> b(k, 0);
  // sum of all ssizes
  double sum_ssizes = 0;
  for(int i=0; i<k; i++){
    sum_ssizes+= sqrt((double)ssizes[i]);
    total_sum+=ssizes[i];
  }

  for(int i=0; i<k; i++){
    b(i) = sqrt((double)ssizes[i])/(sum_ssizes);
  }

  // calculate emat's deviation from pooled emat
  ublas_frame<double> U(k, emat[0].size(), 0);
  for(int i=0; i<k; i++){
    for(int j=0; j<emat[0].size(); j++){
      U(i,j) = emat[i](j) - (b(i) * pooled_emat(j));
    }
  }

  // get covariance matrix for b
  ublas_frame<double> C = covariance_matrix(b);

  // find eigen values for C
  // convert C to arma::mat
  arma::mat arma_C = arma::zeros<arma::mat>(C.size1(), C.size2());
  for(int i=0; i<C.size1(); i++){
    for(int j=0; j<C.size2(); j++){
      arma_C(i,j) = C(i,j);
    }
  }

  // get eigen values
  arma::vec eigval;
  arma::mat eigvec;

  // rank of matrix C is k-1 (proved)
  int R = k-1;

  // covariance matrices should be symmetric and hermitian
  arma::eig_sym(eigval, eigvec, arma_C);

  // Find non-zero eigen values
  vec<int> nz_eindex;
  for(int i=0; i<R; i++){
    temp = eigval.n_elem - 1 - i;
    nz_eindex.push_back(temp);
  }

  // non zero eigen vector(s)
  ublas_frame<double> S(eigvec.n_rows, nz_eindex.size(), 0);
  for(int i=0; i<eigvec.n_rows; i++){
    for(int j=0; j<nz_eindex.size(); j++){
      S(i,j) = eigvec(i,nz_eindex[j]);
    }
  }

  // diagonal matrix with 1/non-zero eigen values
  ublas_vec<double> one_over_nz_eig(nz_eindex.size(),0);
  for(int i=0; i<nz_eindex.size(); i++){
    one_over_nz_eig(i) = 1/sqrt(eigval(nz_eindex[i]));
  }

  ublas_frame<double> Z = diagonal_mat(one_over_nz_eig, nz_eindex.size());

  // compute Mahanalobis distance
  ublas_frame<double> dist_mat = ublas::prod(ublas_frame<double>(ublas::prod(ublas::trans(U), S)),
                                          Z);
  for(int i=0; i<dist_mat.size1(); i++){
    for(int j=0; j<dist_mat.size2(); j++){
      stat += dist_mat(i,j) * dist_mat(i,j);
    }
  }

  // degrees of freedom
  df = R * pooled_emat.size();

  // effect size
  estimate = sqrt(stat/(df*total_sum));

  // chi-squared p-value
  boost::math::chi_squared sharma_song_dist(df);
  if(isnan(stat)){
    stat = 0;
  }
  pvalue = boost::math::cdf(boost::math::complement(sharma_song_dist, stat));
  return stat;

}


// given tables under different experimental conditions, obtain standard normal vectors
std::vector<ublas_vec<double > > get_e_matrix(const std::vector<frame<int > > & tables, vec<int> & ssizes){

  int k = tables.size(); // number of experimental conditions
  std::vector<ublas_vec<double> > e_list(k);

  // a strong assumption here is that all tables[i] should have same dimension
  int nrows = tables[0].size();
  int ncols = tables[0][0].size();

  // temporary data-structures needed
  frame<double> temp_expec(nrows, vec<double>(ncols, 0));
  ublas_frame<double> temp_A(nrows, ncols, 0);
  vec<int> temp_rowsum(nrows);
  vec<int> temp_colsum(ncols);
  ublas_vec<double> temp_prdot(nrows);
  ublas_vec<double> temp_pcdot(ncols);
  ublas_vec<double> ei_vector((nrows-1) * (ncols-1));
  int temp_n;

  for(int i=0; i<k; i++){

    // get dims for each table
    DimSums(tables[i], temp_rowsum, temp_colsum, temp_n);

    // calculate expected for each table
    temp_expec = pchisq_expec(tables[i], temp_rowsum, temp_colsum, temp_n);

    // calculate deviation from expected table
    for(int j=0; j<nrows; j++){
      for(int l=0; l<ncols; l++){
        if(temp_expec[j][l]!=0){
          temp_A(j,l) = (tables[i][j][l] - temp_expec[j][l])/(sqrt(temp_expec[j][l]));
        } else {
          temp_A(j,l) = 0;
        }
      }
    }

    // calculate prdot and pcdot
    for(int j=0; j<nrows; j++){
      temp_prdot(j) = (double)temp_rowsum[j] / temp_n;
    }

    for(int j=0; j<ncols; j++){
      temp_pcdot(j) = (double)temp_colsum[j] / temp_n;
    }

    // calculate V and W
    ublas_frame<double> temp_V = helmert_transform(temp_prdot);
    ublas_frame<double> temp_W = helmert_transform(temp_pcdot);

    // product of row and column helmert matrix with normalized sampled data
    ublas_frame<double> ei = boost::numeric::ublas::prod(ublas_frame<double>(boost::numeric::ublas::prod(temp_V, temp_A)),
                                                         trans(temp_W));
    // flatten ei to ei_vector

    for(int j=1; j<ei.size1(); j++){
      for(int l=1; l<ei.size2(); l++){
        ei_vector(((j-1)*(ei.size2()-1))+(l-1)) = ei(j,l);
      }
    }

    e_list[i] = ei_vector;
    ssizes[i] = temp_n;

  }

  return e_list;
}

// Given a vector of prob., return its helmert transform
ublas_frame<double> helmert_transform(const ublas_vec<double> & p){

  int n = p.size();
  ublas_frame<double> hel_mat(n, n, 0);

  // 1st row is sqrt(p)
  for(int j=0; j<hel_mat.size2(); j++){
    hel_mat(0,j) = sqrt(p(j));
  }

  // if p has only 1 probability return a 1x1 hel_mat
  if(n < 2){
    return hel_mat;
  }

  // summmation of p[i - 1]
  double sump_iminus1 = p(0);
  // summation of p[i]
  double sump_i = 0;

  for(int i=1; i<hel_mat.size1(); i++){
    sump_i = sump_iminus1 + p(i);

    // diagonal elements
    if(sump_iminus1 == 0){
      hel_mat(i,i) = 0;
    } else {
      hel_mat(i,i) = -sqrt(sump_iminus1/sump_i);
    }

    for(int j=0; j<i; j++){
      // other elements
      if((p(i)*p(j)) == 0){
        hel_mat(i,j) = 0;
      } else {
        hel_mat(i,j) = sqrt((p(i)*p(j))/(sump_i*sump_iminus1));
      }
    }

    // update sump_iminus1
    sump_iminus1 = sump_i;

  }

  return hel_mat;
}


// Given a vector of coefficients, return the covariance matrix
ublas_frame<double> covariance_matrix(const ublas_vec<double> & b){

  // number of experimental conditions
  int k = b.size();

  // create intermediate matrices
  ublas_frame<double> jk(k, k, 1);

  ublas_vec<double> veck(1, k);
  ublas_frame<double> diag_k = diagonal_mat(veck);

  ublas_frame<double> diag_b = diagonal_mat(b);

  // pre-computing k * b^T * b
  ublas_frame<double> kbtb(b.size(), b.size(), 0);
  for(int i=0; i<b.size(); i++){
    for(int j=0; j<b.size(); j++){
      kbtb(i,j) = k*b(i)*b(j);
    }
  }

  // covariance matrix
  ublas_frame<double> cov_mat = (diag_k) -
                                boost::numeric::ublas::prod(jk, diag_b) -
                                boost::numeric::ublas::prod(diag_b, jk) +
                                kbtb;

  return cov_mat;
}



// Given a vector of size 'n', construct a nxn matrix
// with values of vector on diagonals
// If the size of vector is 1 with value 'm', construct a mxm identity matrix
ublas_frame<double> diagonal_mat(const ublas_vec<double> & n, int nrows){

  // check if the size of vector is equal to 1
  if(n.size() == 1 && nrows == -1){
    int num = (int)n(0);

    if(num == 0){
      boost::numeric::ublas::zero_matrix<double> emp(0,0);
      return emp;
    } else {
      // create a n[0]xn[0] identity matrix
      boost::numeric::ublas::identity_matrix<double> diag_mat(n(0));
      return diag_mat;
    }
  } else if(n.size() > 1) { // if the size of vector is greater than 1

    // create 0 matrix
    ublas_frame<double> diag_mat(n.size(), n.size(), 0);

    for(int i=0; i<n.size(); i++){
      diag_mat(i,i) = n(i);
    }
    return diag_mat;
  } else { // if size of vector is still 1 but nrows and ncols are not -1

    // create 0 matrix
    ublas_frame<double> diag_mat(nrows, nrows, 0);

    for(int i=0; i<nrows; i++){
      diag_mat(i,i) = n(0);
    }
    return diag_mat;

  }

}
