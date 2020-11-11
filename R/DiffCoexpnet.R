# Coexpnet.R
#
# Created by Sajal Kumar
# Copyright (c) NMSU Song lab

#' A parallel differential co-expression network builder
#'
#' Builds a differential co-expression network based on the Sharma-Song test. 
#' The builder computes a Sharma-Song chi-squared p-value and an effect size for each pair of gene,
#' under different experimental conditions.The builder provides options to discretize the expression 
#' profiles by either Ckmeans.1d.dp, that is an optimal univariate clustering algorithm, or 
#' GridOnCluster, that discretizes genes involved in co-expression patterns across the conditions, 
#' based on their joint distribution. 
#'
#' @import Rcpp
#' @importFrom stats p.adjust mad sd
#' @useDynLib DiffXCoExpNet
#'
#' @param exp_matr a sample (rows) x gene (columns) continuous expression profile matrix representing 
#' gene expression at given samples. This is a combined expression profile across multiple experimental
#' conditions.
#' 
#' @param n_conditions an integer specifying the number of unique experimental conditions.
#' 
#' @param conditions a vector of the same size as the number of samples in \code{expr_matr}, containing (1 based)
#' experimental condition ids for each sample.
#' 
#' @param indices a matrix of all gene interactions to study. The matrix containins the indices (1 based) 
#' of genes in \code{expr_matr}. The first column represents the parent genes and the second column represents 
#' child genes. 
#' 
#' @param g_names a vector containing names of all genes in \code{expr_matr}.
#' 
#' @param method a string specifying which discretization method to use. The options are:
#' 'univariate' that will choose \code{Ckmeans.1d.dp} or 'joint' that will choose
#' \code{GridOnCluster}.
#'
#' @param k a range of clusters to be found on \code{exp_matr}. The default minimum is 2 and the 
#' maximum is automatically computed to require atleast five points per cluster per experimental condition. 
#' If the \code{method} picked is \code{'joint'}, then a \code{Silhouette} score determines the number of 
#' \code{k}; if the \code{method} picked  is \code{'univariate'}, then a \code{BIC} score determines the 
#' number of \code{k}. 
#' 
#' @param nthreads an integer specifying the number of parallel threads to be used. If not specified,
#' a non-parallel version of the code would be executed.
#'
#'
#' @return
#'
#' A data frame containing information on all possible parent x child differential co-expression pairs. 
#' The following items are reported about each gene pair:
#' \item{PID}{the one(1)-based column index in \code{expr_matr} matrix of the parent gene.}
#' 
#' \item{CID}{the one(1)-based column index in \code{expr_matr} matrix of the child gene.}
#'
#' \item{PNAMES}{the name of the parent gene.}
#' 
#' \item{CNAMES}{the name of the child gene.}
#' 
#' \item{PVALUE}{the adjusted sharma-song chi-squared p-value, depicting the amount of interaction 
#' heterogeneity in the gene pair across the experimental conditions.}
#' 
#' \item{ESTIMATE}{the sharma-song effect size, depicting the amount of interaction heterogeneity in
#' the gene pair across the experimental conditions.}
#' 
#' @author 
#' 
#' Ruby Sharma, Sajal Kumar and Mingzhou Song
#' 
#' @seealso
#'
#' See \pkg{Ckmeans.1d.dp} for discretizing univariate continuous data.
#' See \pkg{GridOnClusters} for joint discretization of continuous data.
#'
#' @export

Build_DiffCoexpnet = function(exp_matr, 
                              n_conditions,
                              conditions,
                              indices,
                              g_names,
                              method = c('univariate', 'joint'), 
                              k=NULL,
                              nthreads=NULL){
  
  # check if length of g_names is equal to the number of genes in exp_matr
  if(length(g_names) != ncol(exp_matr)){
    stop("The length of g_names should be equal to the number of genes in exp_matr.")
  }
  
  # check if length of conditions is equal to the number of samples in exp_matr
  if(length(conditions) != nrow(exp_matr)){
    stop("The length of g_names should be equal to the number of samples in exp_matr.")
  }
  
  # check if the number of unique conditions is equal to n_conditions
  if(length(unique(conditions)) != n_conditions){
    stop("The number of unique conditions should be equal to n_conditions.")
  }
  
  # check if the condition ids are correct
  if(!all(sort(unique(conditions)) == c(1:n_conditions))){
    stop("The condition ids should be 1 based and have maximum n_condition levels.")
  }
  
  # check the range of indices
  if(min(indices) < 1 || max(indices) > ncol(exp_matr)){
    stop("Indices should be 1 based and within ncol(exp_matr) (inclusive).")
  }
  
  method = match.arg(method)
  
  if(method == 'univariate'){
    
    # separate data by conditions
    cond_exp_matr = vector("list", length = n_conditions)
    
    for(i in 1:n_conditions){
      cond_exp_matr[[i]] = exp_matr[conditions == i, ]
    }
    
    # scale and shift data
    for(i in 1:n_conditions){
      cond_exp_matr[[i]] = StdScaledTr(cond_exp_matr[[i]])
    }
    
    # merge data
    tr_exp_matr = matrix(0, nrow=nrow(exp_matr), ncol=ncol(exp_matr))
    j = 1
    for(i in 1:n_conditions){
      tr_exp_matr[j : (j + nrow(cond_exp_matr[[i]]) - 1), ] = cond_exp_matr[[i]]
      j = j + nrow(cond_exp_matr[[i]])
    }
    
    # Discretize tr_exp_matr using Ckmeans.1d.dp
    dis_exp_matr = DiscretizeData(data = tr_exp_matr, ncores = nthreads, n_cond = n_conditions)
    
    # find number of levels for all genes
    glevel = apply(dis_exp_matr, 2, max)
    
    # offset by 1
    dis_exp_matr = dis_exp_matr - 1
    conditions = conditions - 1
    indices = indices - 1
    
    # build differential coexpression network (not filtration needed)
    diffcoexpnet_result = Discrete_DiffCoexpnet(exp_matr = dis_exp_matr, 
                                                n_conditions = n_conditions,
                                                conditions = conditions,
                                                indices = indices, 
                                                g_names = g_names,
                                                g_levels = glevel,
                                                nthreads = nthreads)
    
  } else {
    
    # calculate maximum possible cluster
    max_possible = ceiling(sqrt((nrow(exp_matr)/5))/n_conditions)
    
    # check if user has provided k
    if(!is.null(k)){
      
      # check if the length of k is 2
      if(length(k) == 2){
        
        # check if the range is valid
        if(k[1] < 2 || k[2] < 2 || k[1] > max(max_possible,2) || k[2] > max(max_possible,2)){
          stop(paste0("Invalid range of supplied in k, the allowed range is between 2 to ", 
                      max(max_possible,2),"."))
        }
        
      } else {
        
        stop("Length of k should be 2. The first and second element should denote the minimum 
             and maximum number of clusters to consider, respectively.")
        
      }   
    } else { # determine the range of k
      
      k = c(2, max(max_possible, 2))
      if(k[2] == 2){
        k = 2
      }
      
    }
    
    # offset by 1
    conditions = conditions - 1
    indices = indices - 1
    
    # build differential coexpression network (no filtration needed)
    diffcoexpnet_result = DiffCoexpnet(exp_matr = exp_matr,
                                       n_conditions = n_conditions,
                                       conditions = conditions,
                                       indices = indices,
                                       g_names = g_names,
                                       k = k,
                                       nthreads = nthreads)
  }
  
  # adjust p-values
  diffcoexpnet_result$PVALUE = p.adjust(diffcoexpnet_result$PVALUE, method = "BH")
  
  # return differential coexpression network
  return(diffcoexpnet_result)
}

