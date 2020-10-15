#include <Rcpp.h>
#include "methods.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
DataFrame Discrete_Coexpnet(IntegerMatrix children, IntegerMatrix parents, StringVector c_names, StringVector p_names,
                            int nthreads)
{

  // prepare data
  vec<std::string> _c_names = as<std::vector<std::string > >(c_names);
  vec<std::string> _p_names = as<std::vector<std::string > >(p_names);
  
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
  vec<stat_collector> _stat_collect= discrete_coexpnet(_children, _parents, _c_names, _p_names, nthreads);

  // convert stats to Dataframe
  IntegerVector pid, cid;
  StringVector pn, cn;
  NumericVector pvalue, estimate;
  
  for(int i=0; i<_stat_collect.size(); i++){
    pid(i) = _stat_collect[i].pindex;
    cid(i) = _stat_collect[i].cindex;
    pn(i) = _stat_collect[i].pname;
    cn(i) = _stat_collect[i].cname;
    pvalue(i) = _stat_collect[i].pvalue;
    estimate(i) = _stat_collect[i].estimate;
  }
  
  DataFrame stat_collect = DataFrame::create(Named("PID")=pid,
                                             Named("CID")=cid,
                                             Named("PNAMES")=pn,
                                             Named("CNAMES")=cn,
                                             Named("PVALUE")=pvalue,
                                             Named("ESTIMATE")=estimate);
  // 
  // NumericMatrix stat_collect(_stat_collect[0].size(), 4);
  // NumericVector sttemp;
  // for(int i=0; i<4; i++){
  //   sttemp = wrap(_stat_collect[i]);
  //   stat_collect(_,i) = sttemp;
  // }

  return stat_collect;

}

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

  ldouble pval;
  double estimate;
  double stat = SharmaSongTest(_tables, pval, estimate);

  Rcout<<"Statistic:"<<stat<<std::endl;
  Rcout<<"P-value:"<<pval<<std::endl;

}

// //' @export
// // [[Rcpp::export]]
// NumericMatrix Discrete_Diffcoexpnet(IntegerMatrix exp1, IntegerMatrix exp2, IntegerVector parents, IntegerVector children,
//                                     IntegerVector plevels, IntegerVector clevels, int nthreads,
//                                     std::string method="differential")
// {
// 
//   // prepare data
//   vec<int> _parents = as<vec<int > >(parents);
//   vec<int> _children = as<vec<int > >(children);
//   vec<int> _plevels = as<vec<int > >(plevels);
//   vec<int> _clevels = as<vec<int > >(clevels);
// 
//   frame<int> _exp1;
//   vec<int> temp;
//   IntegerVector ptemp;
//   for(int i=0; i<exp1.nrow(); i++){
//     ptemp = exp1(i,_);
//     temp = as<vec<int > >(ptemp);
//     _exp1.push_back(temp);
//   }
// 
//   frame<int> _exp2;
//   for(int i=0; i<exp2.nrow(); i++){
//     ptemp = exp2(i,_);
//     temp = as<vec<int > >(ptemp);
//     _exp2.push_back(temp);
//   }
// 
//   // Get stats
//   frame<mydouble> _stat_collect= discrete_diffcoexpnet(_exp1, _exp2, _parents, _children, _plevels, _clevels,
//                                                        nthreads, method);
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
