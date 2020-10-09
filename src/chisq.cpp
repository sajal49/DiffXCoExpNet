#include "methods.h"

double chisq(const frame<int> & O, ldouble & p_value, double & estimate){

  double chsq = 0.0;

  if (O.size() == 0) { // if number of rows is 0 then return 0.
    return chsq;
  } else if(O[0].size() == 0) { // if number of columns is 0 then return 0.
    return chsq;
  }

  size_t nrows = O.size();  // number of rows
  size_t ncols = O[0].size();  // number of columns

  // get dim sums
  vec<int> rowsums(nrows, 0);
  vec<int> colsums(ncols, 0);
  
  int n = 0;
  DimSums(O, rowsums, colsums, n);

  if(n == 0){ // if no samples are present then return 0.
    return chsq;
  }

  // get expected
  frame<double> eij = pchisq_expec(O, rowsums, colsums, n);

  size_t df = 0; // degrees of freedom

  // compute chisq
  for (size_t i=0; i<nrows; ++i) {
    // Expected count for cell (i,j):
    for (size_t j=0; j<ncols; ++j) {
      if(eij[i][j] > 0){
        chsq += (O[i][j] - eij[i][j]) * (O[i][j] - eij[i][j]) / eij[i][j];
      }
    }
  }

  if(chsq < 0){  // It is possible to get a really small negative value
    chsq = 0;
  }

  // compute estimate -- Cramer's V
  double maxchsq = 0.0;
  if(nrows < ncols){
    maxchsq = (n * (nrows - 1));
  } else {
    maxchsq = (n * (ncols - 1));
  }

  if(maxchsq > 0) {
    estimate = std::sqrt(std::abs(chsq) / maxchsq);
  } else {
    estimate = 0;
  }

  // compute p-value
  df = (nrows - 1) * (ncols - 1);
  boost::math::chi_squared chidist(df);
  p_value = boost::math::cdf(boost::math::complement(chidist, chsq));

  return chsq;

}
