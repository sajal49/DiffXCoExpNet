// Created by Sajal Kumar
// Copyright (c) NMSU Song lab

#include "methods.h"
#include <boost/thread.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <RcppArmadillo.h>

/*
 * Slave -- Parallel discrete coexpression network builder that works with continuous values.
 * 
 * Parameters
 * 
 * children: A continuous sample x gene matrix representing 
 * discrete levels of gene expression at given samples. This matrix
 * contains all child genes.
 * 
 * parents: A continuous sample x gene matrix representing 
 * discrete levels of gene expression at given samples. This matrix
 * contains all parent genes.
 * 
 * c_names: a vector of strings containing gene names for the children matrix.
 * 
 * c_indices: a vector of integers representing the memory location in 'sc' that the 
 * slave process should work on. 
 * 
 * p_names: a vector of strings containing gene names for the parents matrix.
 * 
 * k: a vector of integers that can contain the following:
 * i) k with size two: the minimum and maximum number of discrete levels allowed.
 * ii) k with size one: the maximum number of discrete levels allowed (minimum would be assumed to be 2).
 * iii) k with size n (where n is nparents x nchildren): the number of discrete levels to be used for each
 * nparents x nchildren interactions.
 *
 * sc: An instance of stat_collector to gather results.
 * 
 */

void coexpnet_thread(const frame<double> & children, 
                     const frame<double> & parents, 
                     const vec<std::string> & c_names,
                     vec<int> c_indices, 
                     const vec<std::string> & p_names,
                     const vec<int> & k,
                     vec<stat_collector> & sc){
  
  // dims for input parents
  int nparents = parents.size();
  int nchildren_t = c_indices.size();
  
  // temporary table
  frame<int> temp_tab;
  
  // temporary frame
  frame<double> temp_frame;
  
  // temporary variables for kmeans
  size_t clusters;
  arma::Row<size_t> assignments;
  
  // temporary variables for stat collection
  mydouble stat, p_value, estimate, effind;
  
  for(int c=0; c<nchildren_t; c++) {
    
    // collect result for all parents except which c_names[c_indices[c]] == p_names[i]
    for(int i=0; i<nparents; i++){
      
      // effective index
      effind = (nparents * c_indices[c]) + i;
      
      if(c_names[c_indices[c]] == p_names[i]){
        
        // Not applicable -- store worst result
        fill_stat_collector(sc[effind], i, c_indices[c], p_names[i], c_names[c_indices[c]], 1, 0);
        
      } else {
        
        // if size of k is 2 then it specifies the maximum and minimum range.
        if(k.size() == 2){
          
          // TO-DO  
        
        } else if(k.size() == 1){ // if size of k is 1 then it specifies the maximum range
          
          // TO-DO
          
        } else if(k.size() == nparents * children.size()){ 
          // if size of k is nparents * nchildren, then it specifies the exact k for each combination
          
          temp_frame.push_back(parents[i]);
          temp_frame.push_back(children[c_indices[c]])
          
          // create an arma matrix consisting data
          arma::mat data(temp_frame);
          
          KMeans<> k(500);
          k.Cluster(data, clusters, assignments);
          
        }
        
        // make a contingency table
        temp_tab = tableCpp(parents[i], children[c_indices[c]], -1, -1);
        
        // get stats
        stat = chisq(temp_tab, p_value, estimate);
        
        // store stats
        fill_stat_collector(sc[effind], i, c_indices[c], p_names[i], c_names[c_indices[c]], p_value, estimate);
        
      }
      
    } // parents
    
  } // children
  
}


/*
 * Master -- Parallel discrete coexpression network builder that works with continuous values.
 * 
 * Parameters
 * 
 * children: A continuous sample x gene matrix representing 
 * discrete levels of gene expression at given samples. This matrix
 * contains all child genes.
 * 
 * parents: A continuous sample x gene matrix representing 
 * discrete levels of gene expression at given samples. This matrix
 * contains all parent genes.
 * 
 * c_names: a vector of strings containing gene names for the children matrix.
 * 
 * p_names: a vector of strings containing gene names for the parents matrix.
 * 
 * k: a vector of integers that can contain the following:
 * i) k with size two: the minimum and maximum number of discrete levels allowed.
 * ii) k with size one: the maximum number of discrete levels allowed (minimum would be assumed to be 2).
 * iii) k with size n (where n is nparents x nchildren): the number of discrete levels to be used for each
 * nparents x nchildren interactions.
 * 
 * nthreads: an integer representing the numbers of thread to be used.
 * 
 */

vec<stat_collector> coexpnet(const frame<ldouble> & children, 
                             const frame<ldouble> & parents, 
                             const vec<std::string> & c_names,
                             const vec<std::string> & p_names,
                             const vec<int> & k,
                             int nthreads=-1) {
  
  // fix nthreads
  int true_threads = boost::thread::hardware_concurrency()-1;
  
  if(nthreads > true_threads) {
    nthreads = true_threads;
  }
  
  if(nthreads < 2) {
    std::cout <<"WARNING: No multithreading enabled!"<<std::endl;
  }
  
  // create threads
  boost::thread_group threads;
  
  // dims for input
  int nchildren = children.size();
  int nparents = parents.size();
  
  // a vector of stat_collector data structure
  vec<stat_collector> sc(nparents*nchildren);
  
  // load balancing on number of children for multi-threads
  
  frame<int> c_indices;
  vec<int> temp_indices;
  int eq_child = (int)ceil((double)nchildren/nthreads); // total children for each thread
  
  for(int child_index=0; child_index!=nchildren; child_index++){
    if((child_index+1) % eq_child != 0){ // prepare child indices for each thread
      temp_indices.push_back(child_index);
    } else {
      temp_indices.push_back(child_index);
      c_indices.push_back(temp_indices); // move on to next thread
      temp_indices.clear();
    }
  }
  if(temp_indices.size()!=0){
    c_indices.push_back(temp_indices); // residual child indices
  }
  
  // Create threads for each set of children
  for(int i=0; i<c_indices.size(); i++){
    
    if(nthreads >= 2){ // multiple threads
      threads.create_thread(boost::bind(&coexpnet_thread,
                                        boost::cref(children),
                                        boost::cref(parents),
                                        boost::cref(c_names),
                                        boost::cref(c_indices[i]),
                                        boost::cref(p_names),
                                        boost::cref(k),
                                        boost::ref(sc)));
    } else { // single thread
      coexpnet_thread(children, parents, c_names, c_indices[i], p_names, k, sc);
    }
    
  }
  
  // join threads
  if(nthreads >= 2){
    threads.join_all();
  }
  
  return sc;
}

/*
 * Slave -- Parallel discrete coexpression network builder
 * Parameters
 * 
 * children: A discrete sample x gene matrix representing 
 * discrete levels of gene expression at given samples. This matrix
 * contains all child genes.
 * 
 * parents: A discrete sample x gene matrix representing 
 * discrete levels of gene expression at given samples. This matrix
 * contains all parent genes.
 * 
 * c_names: a vector of strings containing gene names for the children matrix
 * 
 * c_indices: a vector of integers representing the memory location in 'sc' that the 
 * slave process should work on. 
 * 
 * p_names: a vector of strings containing gene names for the parents matrix
 * 
 * sc: An instance of stat_collector to gather results.
 * 
 */

void discrete_coexpnet_thread(const frame<int> & children, 
                              const frame<int> & parents, 
                              const vec<std::string> & c_names,
                              vec<int> c_indices, 
                              const vec<std::string> & p_names, 
                              vec<stat_collector> & sc){

  // dims for input parents
  int nparents = parents.size();
  int nchildren_t = c_indices.size();

  // temporary table
  frame<int> temp_tab;

  // temporary variables for stat collection
  mydouble stat, p_value, estimate, effind;

  for(int c=0; c<nchildren_t; c++) {

   // collect result for all parents except which c_names[c_indices[c]] == p_names[i]
   for(int i=0; i<nparents; i++){

     // effective index
     effind = (nparents * c_indices[c]) + i;

     if(c_names[c_indices[c]] == p_names[i]){

       // Not applicable -- store worst result
       fill_stat_collector(sc[effind], i, c_indices[c], p_names[i], c_names[c_indices[c]], 1, 0);

     } else {

       // make a contingency table
       temp_tab = tableCpp(parents[i], children[c_indices[c]], -1, -1);

       // get stats
       stat = chisq(temp_tab, p_value, estimate);

       // store stats
       fill_stat_collector(sc[effind], i, c_indices[c], p_names[i], c_names[c_indices[c]], p_value, estimate);
       
      }

    } // parents

  } // children

}


/*
 * Master -- Parallel discrete coexpression network builder
 * Parameters
 * 
 * children: A discrete sample x gene matrix representing 
 * discrete levels of gene expression at given samples. This matrix
 * contains all child genes.
 * 
 * parents: A discrete sample x gene matrix representing 
 * discrete levels of gene expression at given samples. This matrix
 * contains all parent genes.
 * 
 * c_names: a vector of strings containing gene names for the children matrix
 * 
 * p_names: a vector of strings containing gene names for the parents matrix
 * 
 * nthreads: an integer representing the numbers of thread to be used.
 * 
 */

vec<stat_collector> discrete_coexpnet(const frame<int> & children, 
                                      const frame<int> & parents, 
                                      const vec<std::string> & c_names,
                                      const vec<std::string> & p_names, 
                                      int nthreads=-1) {

  // fix nthreads
  int true_threads = boost::thread::hardware_concurrency()-1;

  if(nthreads > true_threads) {
     nthreads = true_threads;
  }

  if(nthreads < 2) {
    std::cout <<"WARNING: No multithreading enabled!"<<std::endl;
  }

  // create threads
  boost::thread_group threads;

  // dims for input
  int nchildren = children.size();
  int nparents = parents.size();

  // a vector of stat_collector data structure
  vec<stat_collector> sc(nparents*nchildren);

  // load balancing on number of children for multi-threads

  frame<int> c_indices;
  vec<int> temp_indices;
  int eq_child = (int)ceil((double)nchildren/nthreads); // total children for each thread

  for(int child_index=0; child_index!=nchildren; child_index++){
    if((child_index+1) % eq_child != 0){ // prepare child indices for each thread
      temp_indices.push_back(child_index);
    } else {
      temp_indices.push_back(child_index);
      c_indices.push_back(temp_indices); // move on to next thread
      temp_indices.clear();
    }
  }
  if(temp_indices.size()!=0){
    c_indices.push_back(temp_indices); // residual child indices
  }

  // Create threads for each set of children
  for(int i=0; i<c_indices.size(); i++){

    if(nthreads >= 2){ // multiple threads
      threads.create_thread(boost::bind(&discrete_coexpnet_thread,
                                        boost::cref(children),
                                        boost::cref(parents),
                                        boost::cref(c_names),
                                        boost::cref(c_indices[i]),
                                        boost::cref(p_names),
                                        boost::ref(sc)));
    } else { // single thread
      discrete_coexpnet_thread(children, parents, c_names, c_indices[i], p_names, sc);
    }

  }

  // join threads
  if(nthreads >= 2){
    threads.join_all();
  }

  return sc;
}
