// Created by Sajal Kumar
// Copyright (c) NMSU Song lab

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <numeric>
#include <cstdlib>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// Namespace
namespace ublas = boost::numeric::ublas;
using std::vector;

// Data

// long double
typedef long double ldouble;

// a structure to collect results from coexpnet and diffcoexpnet
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

// Program

// Functions for method_utility.cpp
// Get rowsum, colsum and totalsum for a given table
void DimSums(const frame<int> & O, vec<int> & rowsums, vec<int> & colsums, int & totsum);

// Get pooled table from a list of tables
frame<int> PoolTables(const vec<frame<int > > & obs);

// Get Pearson's chi-square expected table
frame<double> pchisq_expec(const frame<int> & obs, const vec<int> & rowsums, const vec<int> & colsums,
                           int totsum);

// method to fill a container frame row wise, given parent and child
void create_pair_data(frame<double> & container, const vec<double> & dim1, const vec<double> & dim2);

// Container filler for stat_collector structure
void fill_stat_collector(stat_collector & sc, int pindex, int cindex, const std::string & pname,
                         const std::string & cname, ldouble pvalue, double estimate);

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

// transpose a frame
template <typename T>
frame<T> transpose_frame(frame<T> f){

    // For internal use
    frame<T> tf(f[0].size(), vec<T>(f.size(), 0));

    for(int i=0; i<f[0].size(); i++){
        for(int j=0; j<f.size(); j++){
            tf[i][j] = f[j][i];
        }
    }

    return tf;
}

// segment expression matrix by experimental conditions
template <typename T>
vec<frame<T>> Segment_expr_matr_by_conds(const frame<T> & expr_matr, int exp_conds, const vec<int> & conds){

    vec<frame<T>> cond_exp_matr(exp_conds);

    // transpose expr_matr
    frame<T> expr_matr_t = transpose_frame(expr_matr);

    // prepare matrix segments
    for(int i=0; i<conds.size(); i++)
        cond_exp_matr[conds[i]].push_back(expr_matr_t[i]);

    // transpose all matrix segments
    for(int i=0; i<exp_conds; i++)
        cond_exp_matr[i] = transpose_frame(cond_exp_matr[i]);

    return cond_exp_matr;

}

// Append a vector to another vector
template <typename T>
vec<T> append_vec(const vec<T> & vec1, const vec<T> & vec2){
    vec<T> app_vec(vec1.size() + vec2.size(), 0);
    for(int i=0; i<app_vec.size(); i++) {
        if (i < vec1.size()) {
            app_vec[i] = vec1[i];
        } else {
            app_vec[i] = vec2[i - vec1.size()];
        }
    }
    return app_vec;
}

// join expression matrix by experimental conditions
template <typename T>
frame<T> Join_expr_matr_by_conds(const vec<frame<T>> & cond_expr_matr, int exp_conds, int sm){

    frame<T> exp_matr(cond_expr_matr[0].size(), vec<T>(sm, 0));
    vec<T> temp_vec; // temporary vector

    for(int v=0; v<exp_matr.size(); v++){

        // append variable over conditions
        for(int i=1; i<exp_conds; i++){
            if(i == 1){
                temp_vec = append_vec(cond_expr_matr[0][v], cond_expr_matr[1][v]);
            } else {
                temp_vec = append_vec(temp_vec, cond_expr_matr[i][v]);
            }
        }

        exp_matr[v] = temp_vec;
        temp_vec.clear();
    }

    return exp_matr;
}

// Standardize data (scale and shift)
void StandardScalerTr(frame<double> & data);

// Use Kmeans on a n-dimensional data to get cluster assignments
vec<int> GetClusterAssignments(const arma::mat & data, int cluster);

// Apply GridOnCluster on _data using the cluster assignments
frame<int> ApplyGOC(const vec<int> & assignments, const frame<double> & _data, int cluster, int min_levels=2);

// Discretize _data given grid_lines
frame<int> DiscretizeData(const frame<double> & grid_lines, const frame<double> & _data);

// Get mean silhouette score
double GetSilMeanScore(const frame<double> sil_scores);

// Functions for chisq.cpp
// Pearson's chi-squared statistic for co-expression network
double chisq(const vector<vector<int> > & O, ldouble & p_value, double & estimate);

// Functions for coexpnet.cpp
// Worker process -- Discrete co-expression network that discretizes continuous variables on the fly
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

// Worker process -- Discrete co-expression network
void discrete_coexpnet_thread(const frame<int> & children,
                              const frame<int> & parents,
                              const vec<std::string> & c_names,
                              vec<int> c_indices,
                              const vec<std::string> & p_names,
                              vec<stat_collector> & sc);

// Master process -- Discrete co-expression network
vec<stat_collector> discrete_coexpnet(const frame<int> & children,
                                      const frame<int> & parents,
                                      const vec<std::string> & c_names,
                                      const vec<std::string> & p_names,
                                      int nthreads);

// Functions for SharmaSongTest.cpp
// returns standard normal vectors
std::vector<ublas_vec<double > > get_e_matrix(const std::vector<frame<int > > & tables, vec<int> & ssizes);

// helmert transform of probabilities
ublas_frame<double> helmert_transform(const ublas_vec<double> & p);

// get nxn identity matrix
ublas_frame<double> diagonal_mat(const ublas_vec<double> & n, int nrows=-1);

// get covariance matrix
ublas_frame<double> covariance_matrix(const ublas_vec<double> & b);

// Sharma-Song test
ldouble SharmaSongTest(const std::vector<frame<int > > & tables, ldouble & pvalue, double & estimate);

// Functions for diffcoexpnet.cpp
// Worker process -- Discrete differential coexpression network
void discrete_diffcoexpnet_thread(const frame<int> & exp_matr,
                                  int exp_conds,
                                  const vec<int> & conds,
                                  const frame<int> & indices,
                                  const vec<std::string> & gnames,
                                  const vec<int> & glevels,
                                  const vec<int> & int_index,
                                  vec<stat_collector> & sc);

// Master process -- Discrete differential coexpression network
vec<stat_collector> discrete_diffcoexpnet(const frame<int> & exp_matr,
                                          int exp_conds,
                                          const vec<int> & conds,
                                          const frame<int> & indices,
                                          const vec<std::string> & gnames,
                                          const vec<int> & glevels,
                                          int nthreads);

// Worker process -- differential coexpression network with on the fly discretization
void diffcoexpnet_thread(const frame<double> & exp_matr,
                         int exp_conds,
                         const vec<int> & conds,
                         const frame<int> & indices,
                         const vec<std::string> & gnames,
                         const vec<int> & k,
                         const vec<int> & int_index,
                         vec<stat_collector> & sc);

// Master process -- differential coexpression network with on the fly discretization
vec<stat_collector> diffcoexpnet(const frame<double> & exp_matr,
                                 int exp_conds,
                                 const vec<int> & conds,
                                 const frame<int> & indices,
                                 const vec<std::string> & gnames,
                                 const vec<int> & k,
                                 int nthreads);