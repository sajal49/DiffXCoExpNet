#include <Rcpp.h>
#include "methods.h"
using namespace Rcpp;

// helper functions
DataFrame ConvertResult(const vec<stat_collector> & _stat_collect, bool identical);

// [[Rcpp::export]]
DataFrame Discrete_Coexpnet(IntegerMatrix children, IntegerMatrix parents, StringVector c_names, 
                            StringVector p_names, int nthreads, bool identical)
{

  // prepare names
  vec<std::string> _c_names = as<std::vector<std::string > >(c_names);
  vec<std::string> _p_names = as<std::vector<std::string > >(p_names);
  
  // prepare parents
  frame<int> _parents;
  vec<int> temp;
  IntegerVector ptemp;
  for(int i=0; i<parents.ncol(); i++){
    ptemp = parents(_,i);
    temp = as<std::vector<int > >(ptemp);
    _parents.push_back(temp);
  }

  // prepare children
  frame<int> _children;
  IntegerVector ctemp;
  for(int i=0; i<children.ncol(); i++){
    ctemp = children(_,i);
    temp = as<std::vector<int > >(ctemp);
    _children.push_back(temp);
  }

  // Get stats
  vec<stat_collector> _stat_collect= discrete_coexpnet(_children, _parents, _c_names, _p_names, nthreads);
  
  // Get converted result
  DataFrame stat_collect = ConvertResult(_stat_collect, identical);
  
  return stat_collect;
}

// [[Rcpp::export]]
DataFrame Coexpnet(NumericMatrix children, NumericMatrix parents, StringVector c_names, 
                   StringVector p_names, IntegerVector k, int nthreads, bool identical)
{
  
  // prepare names
  vec<std::string> _c_names = as<std::vector<std::string > >(c_names);
  vec<std::string> _p_names = as<std::vector<std::string > >(p_names);
  
  // prepare k
  vec<int> _k = as<std::vector<int > >(k);
  
  // prepare parents
  frame<double> _parents;
  vec<double> temp;
  NumericVector ptemp;
  for(int i=0; i<parents.ncol(); i++){
    ptemp = parents(_,i);
    temp = as<std::vector<double > >(ptemp);
    _parents.push_back(temp);
  }
  
  // prepare children
  frame<double> _children;
  NumericVector ctemp;
  for(int i=0; i<children.ncol(); i++){
    ctemp = children(_,i);
    temp = as<std::vector<double > >(ctemp);
    _children.push_back(temp);
  }
  
  // Get stats
  vec<stat_collector> _stat_collect= coexpnet(_children, _parents, _c_names, _p_names, _k, nthreads);
  
  // Get converted result
  DataFrame stat_collect = ConvertResult(_stat_collect, identical);
  
  return stat_collect;
}


// [[Rcpp::export]]
DataFrame Discrete_DiffCoexpnet(IntegerMatrix exp_matr, int n_conditions, IntegerVector conditions, 
                                IntegerMatrix indices, StringVector g_names, IntegerVector g_levels, 
                                int nthreads)
{
  
  // prepare gene names
  vec<std::string> _gnames = as<std::vector<std::string > >(g_names);
  
  // prepare per sample condition vector
  vec<int> _conds = as<std::vector<int > >(conditions);
  
  // prepare maximum levels for each gene
  vec<int> _glevels = as<std::vector<int > >(g_levels);
  
  // prepare expression matrix
  frame<int> _exp_matr;
  vec<int> temp;
  IntegerVector etemp;
  for(int i=0; i<exp_matr.ncol(); i++){
    etemp = exp_matr(_,i);
    temp = as<std::vector<int > >(etemp);
    _exp_matr.push_back(temp);
  }
  
  // prepare interaction indices
  frame<int> _indices;
  for(int i=0; i<indices.ncol(); i++){
    etemp = indices(_,i);
    temp = as<std::vector<int > >(etemp);
    _indices.push_back(temp);
  }
  
  // Get stats
  vec<stat_collector> _stat_collect= discrete_diffcoexpnet(_exp_matr, n_conditions, _conds, _indices, 
                                                           _gnames, _glevels, nthreads);
  
  // Get converted result
  DataFrame stat_collect = ConvertResult(_stat_collect, false);
  
  return stat_collect;
}

// [[Rcpp::export]]
DataFrame DiffCoexpnet(NumericMatrix exp_matr, int n_conditions, IntegerVector conditions, 
                       IntegerMatrix indices, StringVector g_names, IntegerVector k, 
                       int nthreads)
{
  
  // prepare gene names
  vec<std::string> _gnames = as<std::vector<std::string > >(g_names);
  
  // prepare per sample condition vector
  vec<int> _conds = as<std::vector<int > >(conditions);
  
  // prepare k
  vec<int> _k = as<std::vector<int > >(k);
  
  // prepare expression matrix
  frame<double> _exp_matr;
  vec<double> temp;
  NumericVector etemp;
  for(int i=0; i<exp_matr.ncol(); i++){
    etemp = exp_matr(_,i);
    temp = as<std::vector<double > >(etemp);
    _exp_matr.push_back(temp);
  }
  
  // prepare interaction indices
  frame<int> _indices;
  vec<int> vtemp;
  IntegerVector itemp;
  for(int i=0; i<indices.ncol(); i++){
    itemp = indices(_,i);
    vtemp = as<std::vector<int > >(itemp);
    _indices.push_back(vtemp);
  }
  
  // Get stats
  vec<stat_collector> _stat_collect= diffcoexpnet(_exp_matr, n_conditions, _conds, _indices, 
                                                  _gnames, _k, nthreads);
  
  // Get converted result
  DataFrame stat_collect = ConvertResult(_stat_collect, false);
  
  return stat_collect;
}


DataFrame ConvertResult(const vec<stat_collector> & _stat_collect, bool identical){
  
  // convert stats to Dataframe
  IntegerVector pid(_stat_collect.size()), cid(_stat_collect.size());
  StringVector pn(_stat_collect.size()), cn(_stat_collect.size());
  NumericVector pvalue(_stat_collect.size()), estimate(_stat_collect.size());
  
  int r=0, c=0;
  int nparents = 0;
  
  // if identical calculate the number of parents
  nparents = std::sqrt(_stat_collect.size());
  
  for(int i=0; i<_stat_collect.size(); i++){
    
    
    if(!identical || c > r){
      
      pid(i) = _stat_collect[i].pindex;
      cid(i) = _stat_collect[i].cindex;
      pn(i) = _stat_collect[i].pname;
      cn(i) = _stat_collect[i].cname;
      pvalue(i) = _stat_collect[i].pvalue;
      estimate(i) = _stat_collect[i].estimate;
      
    } else {
      
      pid(i) = NA_INTEGER;
      cid(i) = NA_INTEGER;
      pn(i) = NA_STRING;
      cn(i) = NA_STRING;
      pvalue(i) = NA_REAL;
      estimate(i) = NA_REAL;
      
    }
    
    if(identical){
      
      if((c % (nparents-1)) == 0 && c != 0){
        r = r + 1;
        c = 0;
      } else {
        c = c + 1;
      }
      
    }
  }
  
  // return named dataframe
  DataFrame stat_collect = DataFrame::create(Named("PID")=pid,
                                             Named("CID")=cid,
                                             Named("PNAMES")=pn,
                                             Named("CNAMES")=cn,
                                             Named("PVALUE")=pvalue,
                                             Named("ESTIMATE")=estimate);
  return(stat_collect);
  

}


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

// [[Rcpp::export]]
List SharmaSongTestRcpp(List tables, int nrows, int ncols){

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

  // return variables
  ldouble pval;
  double estimate;
  size_t df;
  double stat = SharmaSongTest(_tables, pval, estimate, df);
  
  return(List::create(Rcpp::Named("statistics")=stat,
                      Rcpp::Named("p.value")=pval,
                      Rcpp::Named("estimate")=estimate,
                      Rcpp::Named("parameters")=(int)df));

}