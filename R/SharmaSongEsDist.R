# SharmaSongEsDist.R
#
# Created by Sajal Kumar
# Copyright (c) NMSU Song lab

#' Sharma-Song effect size distribution generator
#'
#' This generator simulates differential interactions (alternate population) and computes their Sharma-Song effect
#' size. The resultant distribution can be used to threshold effect size as well as calculate the statistical
#' power.
#'
#' @import Rcpp
#' @import doParallel
#' @importFrom FunChisq simulate_tables
#' @importFrom parallel parLapply makeForkCluster stopCluster detectCores
#' @importFrom stats p.adjust mad sd
#'
#' @param N a vector of sample size for each experimental conditions.
#' 
#' @param n_conditions number of experimental conditions.
#' 
#' @param nthreads number of parallel threads to use.
#' 
#' @param num number of simulations to consider (default:100000)
#'
#' @param d the maximum dimension (row and column) for simulating tables.
#'
#' @return
#'
#' A data frame containing the SharmaSong test statistic, p-value and effect size for all \code{num} 
#' simulations.
#' 
#' @author 
#' 
#' Ruby Sharma, Sajal Kumar and Mingzhou Song
#' 
#' @seealso
#'
#' See \pkg{DiffXTables} for computing the SharmaSong test statistic for n tables.
#'
#' @export

SharmaSongEsDist = function(N, n_conditions, nthreads, d, num=100000){
  
  # get simulated tables
  tb_list = SimulateDiffInteractions(N = N, n_conditions = n_conditions, nthreads = nthreads, 
                                     num = num, d = d)
  
  # get sharma-song test statistics
  shsong_res = ComputeSharmaSong(tb_list)
  
  # return sharma-song stats
  return(shsong_res)
}

# Uses simulate_tables from FunChisq to generate differential interactions by independently generating
# tables that differ within a specific margin.
#
# N : a vector of sample size for each experimental conditions.
# n_conditions : number of experimental conditions.
# nthreads : number of parallel threads to use.
# num : number of simulations to consider.
# d : if specified, the tables would have dimensions upto d.
# r : if specified, the tables would have exactly r rows.
# s : if specified, the tables would have exactly s columns.
#
SimulateDiffInteractions = function(N, n_conditions, nthreads, num, 
                                    d=NULL, r=NULL, s=NULL){
  
  # check for ncores available
  
  if(nthreads > detectCores() - 2){
    nthreads = detectCores() - 2
  }
  
  if(nthreads <= 1){
    
    # stores tables
    tb_list = vector("list", length = num)
    
    # simulation
    for(i in 1:num){
      
      # determine r and s
      if(is.null(r) && is.null(s) & !is.null(d)){
        
        if(d <= 2){
          r1 = 2
          s1 = 2
        } else {
          r1 = sample(2:d, 1)
          s1 = sample(2:d, 1)  
        }
      
      }
      
      # keep generating tables independently until they are differ within a tolerance
      
      difference = matrix(0, nrow = r1, ncol = s1)
      
      while(sum(difference) < .Machine$double.eps){
        
        # collect tables for simulation i
        tbs = vector("list", n_conditions)
        
        for(ni in 1:n_conditions){
          
          # generate table
          
          if(N[ni] < r1*s1){
            ttypes = sample(c("functional","independent"), 1)
          } else {
            ttypes = sample(c("functional","independent","dependent.non.functional"), 1)
          }
          tbs[[ni]] = simulate_tables(n = N[ni], nrow = r1, ncol = s1, type = ttypes)$sample.list[[1]]
          
          # compute running difference
          if(sum(difference) == 0){
            difference = (tbs[[ni]] - GetExpected(tbs[[ni]]))/N[ni]
          } else {
            difference = abs(difference - ((tbs[[ni]] - GetExpected(tbs[[ni]]))/N[ni]))
          }
            
        }
        
      }
      
      # store tables for simulation i
      tb_list[[i]] = tbs  
      
    }
    
  } else {
    
    # prepare cluster for parallel
    
    cl = makeForkCluster(nnodes = nthreads)
    
    tb_list = parLapply(cl = cl, X = 1:num, fun = function(i){
      
      # set seed to avoid simulating the exact tables (simulate_tables is not thread safe)
      set.seed(i)
      
      # determine r and s
      if(is.null(r) && is.null(s) & !is.null(d)){
        
        if(d <= 2){
          r1 = 2
          s1 = 2
        } else {
          r1 = sample(2:d, 1)
          s1 = sample(2:d, 1)  
        }
        
      }
      
      # keep generating tables independently until they are differ within a tolerance
      
      difference = matrix(0, nrow = r1, ncol = s1)
      
      while(sum(difference) < .Machine$double.eps){
        
        # collect tables for simulation i
        tbs = vector("list", n_conditions)
        
        for(ni in 1:n_conditions){
          
          # generate table
          
          if(N[ni] < r1*s1){
            ttypes = sample(c("functional","independent"), 1)
          } else {
            ttypes = sample(c("functional","independent","dependent.non.functional"), 1)
          }
          tbs[[ni]] = simulate_tables(n = N[ni], nrow = r1, ncol = s1, type = ttypes)$sample.list[[1]]
          
          # compute running difference
          if(sum(difference) == 0){
            difference = (tbs[[ni]] - GetExpected(tbs[[ni]]))/N[ni]
          } else {
            difference = abs(difference - ((tbs[[ni]] - GetExpected(tbs[[ni]]))/N[ni]))  
          }
        }
        
      }
      
      return(tbs)
      
    })
    
    stopCluster(cl)
    
  }
  
  return(tb_list)
}

# returns the Pearson's chi-squared test based expected values for a given table
# tb : a contingency table
GetExpected = function(tb){
  
  # calculate row, column and total sums
  rsum = rowSums(tb)
  csum = colSums(tb)
  tsum = sum(rsum)
  
  # calculate expected
  expected = (rsum %*% t(csum))/tsum
  
  return(expected)
}


# Computes sharma-song test statistic, pvalue, estimate and degrees of freedom for a list of 
# differential tables
# tb_list : a list of differential tables (list of lists of tables).
ComputeSharmaSong = function(tb_list){
  
  # prepare result container
  result = matrix(0, nrow=length(tb_list), ncol=5)
  colnames(result) = c("INDEX", "STATISTIC", "PVALUE", "ESTIMATE", "DF")
  
  for(i in 1:length(tb_list)){
    
    # get table dimensions
    nr = nrow(tb_list[[i]][[1]])
    nc = ncol(tb_list[[i]][[1]])
    
    # get sharma song result
    stat = SharmaSongTestRcpp(tb_list[[i]], nr, nc)
    
    result[i, 1:5] = c(i, stat$statistic, stat$p.value, stat$estimate, stat$parameters)
  }
  
  # return
  return(result)
 
}
