#include <Rcpp.h>
#include "methods.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericMatrix Discrete_Coexpnet(IntegerMatrix children, IntegerMatrix parents, StringVector c_names, StringVector p_names,
                                String method, int nthreads)
{

  // prepare data
  vec<std::string> _c_names = as<std::vector<std::string > >(c_names);
  vec<std::string> _p_names = as<std::vector<std::string > >(p_names);
  std::string _method = method;

  // prepare parents
  frame<int> _parents;
  vec<int> temp;
  IntegerVector ptemp;
  for(int i=0; i<parents.nrow(); i++){
    ptemp = parents(i,_);
    temp = as<std::vector<int > >(ptemp);
    _parents.push_back(temp);
  }

  // prepare children
  frame<int> _children;
  IntegerVector ctemp;
  for(int i=0; i<children.nrow(); i++){
    ctemp = children(i,_);
    temp = as<std::vector<int > >(ctemp);
    _children.push_back(temp);
  }

  // Get stats
  frame<mydouble> _stat_collect= discrete_coexpnet(_children, _parents, _c_names, _p_names, _method, nthreads);

  // convert to numeric matrix
  NumericMatrix stat_collect(_stat_collect[0].size(), 4);
  NumericVector sttemp;
  for(int i=0; i<4; i++){
    sttemp = wrap(_stat_collect[i]);
    stat_collect(_,i) = sttemp;
  }

  return stat_collect;

}

// //' @export
// // [[Rcpp::export]]
// NumericMatrix Discrete_Coexpnet(List children, List parents, StringVector c_names, StringVector p_names,
//                                 String method, int nthreads){
//
//   // prepare data
//   vec<std::string> _c_names = as<std::vector<std::string > >(c_names);
//   vec<std::string> _p_names = as<std::vector<std::string > >(p_names);
//   std::string _method = method;
//
//   // prepare parents and children
//   frame<int> _parents;
//   frame<int> _children
//   vec<int> temp;
//   IntegerVector ptemp;
//   IntegerVector ctemp;
//
//   for(int i=0; i<parents.nrow(); i++){
//     ptemp = as<IntegerVector>(parents[i]);
//     ctemp = as<IntegerVector>(children[i]);
//
//     temp = as<std::vector<int > >(ptemp);
//     _parents.push_back(temp);
//
//     temp = as<std::vector<int > >(ctemp);
//     _children.push_back(temp);
//   }
//
//   // Get stats
//   frame<mydouble> _stat_collect= discrete_coexpnet(_children, _parents, _c_names, _p_names, _method, nthreads);
//
//   // convert to numeric matrix
//   NumericMatrix stat_collect(_stat_collect[0].size(), 4);
//   NumericVector sttemp;
//   for(int i=0; i<4; i++){
//     sttemp = wrap(_stat_collect[i]);
//     stat_collect(_,i) = sttemp;
//   }
//
//   return stat_collect;
//
// }

//' @export
// [[Rcpp::export]]
NumericMatrix tableRcpp(NumericVector x, NumericVector y, int xlevels=-1, int ylevels=-1) {

  //Convert numeric vector to standard vector
  std::vector<int> x_vec = as< std::vector<int> >(x);
  std::vector<int> y_vec = as< std::vector<int> >(y);

  // Get table from CPP
  frame<int> tablecpp = tableCpp(x_vec, y_vec, xlevels, ylevels);

  // Convert to numeric matrix
  NumericMatrix table(tablecpp.size(), tablecpp[0].size());

  for(int i=0; i<tablecpp.size(); i++){
    for(int j=0; j<tablecpp[0].size(); j++){
      table(i,j) = tablecpp[i][j];
    }
  }

  return table;
}

//' @export
// [[Rcpp::export]]
NumericMatrix helmert_trans(NumericVector p){

  // convert to standard vector
  ublas_vec<double> _p(p.size());
  for(int i=0; i<p.size(); i++){
    _p(i) = p(i);
  }

  // Get result
  ublas_frame<double> _hel_mat = helmert_transform(_p);

  // convert to numeric matrix
  NumericMatrix hel_mat(_hel_mat.size1(), _hel_mat.size2());

  for(int i=0; i<_hel_mat.size1(); i++){
    for(int j=0; j<_hel_mat.size2(); j++){
      hel_mat(i,j) = _hel_mat(i,j);
    }
  }

  return hel_mat;
}

//' @export
// [[Rcpp::export]]
NumericMatrix covariance_mat(NumericVector b){

  // convert to standard vector
  ublas_vec<double> _b(b.size());
  for(int i=0; i<b.size(); i++){
    _b(i) = b(i);
  }

  // Get result
  ublas_frame<double> _cov_mat = covariance_matrix(_b);

  // convert to numeric matrix
  NumericMatrix cov_mat(_cov_mat.size1(), _cov_mat.size2());

  for(int i=0; i<_cov_mat.size1(); i++){
    for(int j=0; j<_cov_mat.size2(); j++){
      cov_mat(i,j) = _cov_mat(i,j);
    }
  }

  return cov_mat;
}


//' @export
// [[Rcpp::export]]
List get_e_mat(List tables, int nrows, int ncols){

  List t(tables.size());

  // convert List tables to vector of frame<int>
  std::vector<frame<int > > _tables(tables.size(), frame<int>(nrows, vec<int>(ncols, 0)));

  // Generic temp variable
  GenericVector temp;

  for(int i=0; i<tables.size(); i++){
    temp = tables(i);
    for(int j=0; j<nrows; j++){
      for(int k=0; k<ncols; k++){
        _tables[i][j][k] = temp((j*ncols)+k);
      }
    }
  }

  // prepare ssizes
  std::vector<int> _ssizes(tables.size());

  std::vector<ublas_vec<double > > e_mat = get_e_matrix(_tables, _ssizes);

  for(int i=0; i<e_mat.size(); i++){
    NumericVector temp2(((nrows-1)*(ncols-1)));
    for(int j=0; j<((nrows-1)*(ncols-1)); j++){
      temp2[j] = e_mat[i](j);
    }
    t[i] = temp2;
  }

  return t;

}


//' @export
// [[Rcpp::export]]
void SharmaSongTestRcpp(List tables, int nrows, int ncols){

  List t(tables.size());

  // convert List tables to vector of frame<int>
  std::vector<frame<int > > _tables(tables.size(), frame<int>(nrows, vec<int>(ncols, 0)));

  // Generic temp variable
  GenericVector temp;

  for(int i=0; i<tables.size(); i++){
    temp = tables(i);
    for(int j=0; j<nrows; j++){
      for(int k=0; k<ncols; k++){
        _tables[i][j][k] = temp((k*nrows)+j);
      }
    }
  }

  mydouble pval, estimate;
  mydouble stat = SharmaSongTest(_tables, pval, estimate);

  Rcout<<"Statistic:"<<stat<<std::endl;
  Rcout<<"P-value:"<<pval<<std::endl;

}

//' @export
// [[Rcpp::export]]
void ConservedTestRcpp(List tables, int nrows, int ncols){

  List t(tables.size());

  // convert List tables to vector of frame<int>
  std::vector<frame<int > > _tables(tables.size(), frame<int>(nrows, vec<int>(ncols, 0)));

  // Generic temp variable
  GenericVector temp;

  for(int i=0; i<tables.size(); i++){
    temp = tables(i);
    for(int j=0; j<nrows; j++){
      for(int k=0; k<ncols; k++){
        _tables[i][j][k] = temp((k*nrows)+j);
      }
    }
  }

  mydouble pval, estimate;
  mydouble stat = ConservedTest(_tables, pval, estimate);

  Rcout<<"Statistic:"<<stat<<std::endl;
  Rcout<<"P-value:"<<pval<<std::endl;

}

//' @export
// [[Rcpp::export]]
NumericMatrix Discrete_Diffcoexpnet(IntegerMatrix exp1, IntegerMatrix exp2, IntegerVector parents, IntegerVector children,
                                    IntegerVector plevels, IntegerVector clevels, int nthreads,
                                    std::string method="differential")
{

  // prepare data
  vec<int> _parents = as<vec<int > >(parents);
  vec<int> _children = as<vec<int > >(children);
  vec<int> _plevels = as<vec<int > >(plevels);
  vec<int> _clevels = as<vec<int > >(clevels);

  frame<int> _exp1;
  vec<int> temp;
  IntegerVector ptemp;
  for(int i=0; i<exp1.nrow(); i++){
    ptemp = exp1(i,_);
    temp = as<vec<int > >(ptemp);
    _exp1.push_back(temp);
  }

  frame<int> _exp2;
  for(int i=0; i<exp2.nrow(); i++){
    ptemp = exp2(i,_);
    temp = as<vec<int > >(ptemp);
    _exp2.push_back(temp);
  }

  // Get stats
  frame<mydouble> _stat_collect= discrete_diffcoexpnet(_exp1, _exp2, _parents, _children, _plevels, _clevels,
                                                       nthreads, method);

  // convert to numeric matrix
  NumericMatrix stat_collect(_stat_collect[0].size(), 4);
  NumericVector sttemp;
  for(int i=0; i<4; i++){
    sttemp = wrap(_stat_collect[i]);
    stat_collect(_,i) = sttemp;
  }

  return stat_collect;

}
