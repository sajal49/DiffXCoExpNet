//
// Created by Sajal Kumar.
//
// Copyright (c) NMSU Song lab

#ifndef DIFFXTABLES_COEXPNET_H
#define DIFFXTABLES_COEXPNET_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <numeric>
#include <stdlib.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
using std::vector;

// data structures
// #if defined _WIN32
// typedef double mydouble;
// //typedef long double mydouble;
// #else
typedef long double ldouble;
//#endif

struct stat_collector{
  int pindex; // parent index
  int cindex; // child index
  std::string pname; // parent name
  std::string cname; // child name
  ldouble pvalue; // pvalue
  ldouble estimate; // effect size
};

// Generic frame: vector of vector
template <typename T>
using frame = vector<vector<T> >;

// Generic vector
template <typename T>
using vec = vector<T>;

// Generic ublas matrix
template <typename T>
using ublas_frame = ublas::matrix<T>;

// Generic ublas vector
template <typename T>
using ublas_vec = ublas::vector<T>;

// utility functions
// Get rowsum, colsum and totalsum for a given table
void DimSums(const frame<int> & O, vec<int> & rowsums, vec<int> & colsums, int & totsum);

// Get pooled table from a list of tables
frame<int> PoolTables(const vec<frame<int > > & obs);

// Get Pearson's chi-square expected table
frame<mydouble> pchisq_expec(const frame<int> & obs, const vec<int> & rowsums, const vec<int> & colsums,
                             int totsum);

// a faster table function
frame<int> tableCpp(vector<int> x_vec, vector<int> y_vec, int xlevels, int ylevels);

// function to print tables
template <typename T>
void print_tableCpp(frame<T> table){
  // For internal use

  // column names
  std::cout<<"   ";
  for(int j=0; j<table[0].size(); j++){
    std::cout<<j<<" ";
  }
  std::cout<<std::endl;

  // table
  for(int i=0; i<table.size(); i++){
    // rownames
    std::cout<<i<<": ";
    for(int j=0; j<table[0].size(); j++){
      std::cout<<table[i][j]<<" ";
    }
    std::cout<<std::endl;
  }

  std::cout<<std::endl;

}

// Pearson's chi-squared statistic for co-expression network
double chisq(const vector<vector<int> > & O, ldouble & p_value, double & estimate);


// Slave process -- Discrete co-expression network that discretizes continuous variables on the fly
void coexpnet_thread(const frame<double> & children, 
                     const frame<double> & parents, 
                     const vec<std::string> & c_names,
                     vec<int> c_indices, 
                     const vec<std::string> & p_names,
                     const vec<int> & k,
                     vec<stat_collector> & sc);

// Master process -- Discrete co-expression network that discretizes continuous variables on the fly
vec<stat_collector> coexpnet(const frame<double> & children, 
                             const frame<double> & parents, 
                             const vec<std::string> & c_names,
                             const vec<std::string> & p_names,
                             const vec<int> & k,
                             int nthreads);

// Discrete Coexpression network
void discrete_coexpnet_thread(const frame<int> & children, const frame<int> & parents, const vec<std::string> & c_names,
                              vec<int> c_indices, const vec<std::string> & p_names, std::string method,
                              frame<mydouble> & stat_collect);

frame<mydouble> discrete_coexpnet(const frame<int> & children, const frame<int> & parents, const vec<std::string> & c_names,
                                  const vec<std::string> & p_names, std::string method, int nthreads);


// Sharma-Song Test of 2nd order differential expression
// returns standard normal vectors
std::vector<ublas_vec<double > > get_e_matrix(const std::vector<frame<int > > & tables, vec<int> & ssizes);

// helmert transform of probabilities
ublas_frame<double> helmert_transform(const ublas_vec<double> & p);

// get nxn identity matrix
ublas_frame<double> diagonal_mat(const ublas_vec<double> & n, int nrows=-1);

// get covariance matrix
ublas_frame<double> covariance_matrix(const ublas_vec<double> & b);

// Sharma-Song test
mydouble SharmaSongTest(const std::vector<frame<int > > & tables, mydouble & pvalue, mydouble & estimate);

// Conserved test
mydouble ConservedTest(const std::vector<frame<int > > & tables, mydouble & pvalue, mydouble & estimate);

//Discrete Differential Coexpression Network
frame<mydouble> discrete_diffcoexpnet(const frame<int> & exp1, const frame<int> & exp2, const vec<int> & parents,
                                      const vec<int> & children, const vec<int> & plevels, const vec<int> & clevels,
                                      int nthreads, std::string method);

void discrete_diffcoexpnet_thread(const frame<int> & exp1, const frame<int> & exp2, const vec<int> & parents,
                                  const vec<int> & children, const vec<int> & plevels, const vec<int> & clevels,
                                  const vec<int> & int_index, frame<mydouble> & stat_collect, std::string method);
