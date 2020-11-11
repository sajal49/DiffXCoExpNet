# Coexpnet.R
#
# Created by Sajal Kumar
# Copyright (c) NMSU Song lab

#' A parallel co-expression network builder
#'
#' Builds a co-expression network based on the Pearson's chi-square test. 
#' The builder computes a chi-squared p-value and a cramer's V effect size for each pair of gene.
#' The builder provides options to discretize the expression profiles by either Ckmeans.1d.dp,
#' that is an optimal univariate clustering algorithm, or GridOnCluster, that discretizes genes
#' involved in a co-expression pattern based on their joint distribution. 
#'
#' @import Rcpp
#' @importFrom stats p.adjust mad sd
#' @useDynLib DiffXCoExpNet
#'
#' @param parent_expr a sample (rows) x gene (columns) continuous expression profile matrix containing 
#' the the parent genes. This matrix also represents one half of all desirable co-expression
#' gene pairs.
#' 
#' @param child_expr a sample (rows) x gene (columns) continuous expression profile matrix containing 
#' the the child genes. This matrix also represents the other half of all desirable co-expression
#' gene pairs. The parent_expr and child_expr can be the same expression matrix if
#' a complete graph between all genes in a dataset is required.
#' 
#' @param c_names a vector containing names of the child genes.
#' 
#' @param p_names a vector containing names of the parent genes. The p_names and c_names
#' vector can be the same if a parent_expr and child_expr are equivalent.
#'
#' @param method a string specifying which discretization method to use. The options are:
#' 'univariate' that will choose \code{Ckmeans.1d.dp} or 'joint' that will choose
#' \code{GridOnCluster}.
#'
#' @param k a range of clusters to be found on \code{parent_expr} and \code{child_expr}. 
#' The default minimum is 2 and the maximum is automatically computed to require atleast 
#' five points per cluster. If the \code{method} picked is \code{'joint'}, then a \code{Silhouette} 
#' score determines the number of \code{k}; if the \code{method} picked is \code{'univariate'}, 
#' then a \code{BIC} score determines the number of \code{k}. 
#' 
#' @param nthreads an integer specifying the number of parallel threads to be used. If not specified,
#' a non-parallel version of the code would be executed.
#'
#'
#' @return
#'
#' A data frame containing information on all possible parent x child co-expression pairs. 
#' The following items are reported about each gene pair:
#' \item{PID}{the one(1)-based column index in \code{parent_expr} matrix of the parent gene.}
#' 
#' \item{CID}{the one(1)-based column index in \code{child_expr} matrix of the child gene.}
#'
#' \item{PNAMES}{the name of the parent gene.}
#' 
#' \item{CNAMES}{the name of the child gene.}
#' 
#' \item{PVALUE}{the adjusted chi-squared p-value, depicting the statistical dependency between the gene pair.}
#' 
#' \item{ESTIMATE}{the cramer's v effect size, depicting the strength of the relationship between the 
#' gene pair.}
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

Build_Coexpnet = function(parent_expr, 
                          child_expr, 
                          c_names, 
                          p_names,
                          method = c('univariate', 'joint'), 
                          k=NULL,
                          nthreads=NULL){
  
  # check if length of p_names is equal to the number of genes in parent_expr
  if(length(p_names) != ncol(parent_expr)){
    stop("The length of p_names should be equal to the number of genes in parent_expr.")
  }
  
  # check if length of c_names is equal to the number of genes in child_expr
  if(length(c_names) != ncol(child_expr)){
    stop("The length of c_names should be equal to the number of genes in child_expr.")
  }
  
  method = match.arg(method)
  
  # check if parent_expr and child_expr are equivalent
  if(identical(parent_expr, child_expr)){
    pc_expr_ident = TRUE
  } else {
    pc_expr_ident = FALSE
  }
  
  if(method == 'univariate'){
    
    # Discretize parent_expr using Ckmeans.1d.dp
    dis_parent_expr = DiscretizeData(parent_expr, nthreads)
    
    # if parent_expr and child_expr are equivalent then no need to discretize again
    if(pc_expr_ident){
      
      dis_child_expr = dis_parent_expr
    
    } else {
      
      # Discretize child_expr using Ckmeans.1d.dp
      dis_child_expr = DiscretizeData(child_expr, nthreads)
      
    }
    
    # find if there are any genes with only 1 level (they will removed)
    plevel = apply(dis_parent_expr, 2, max)
    remp = plevel == 1
    if (any(remp)){
      
      parent_expr = parent_expr[ , !remp]
      p_names = p_names[!remp]
      dis_parent_expr = dis_parent_expr[ , !remp]
      
      # child_expr will be treated similarly if identical parent_expr
      if(pc_expr_ident){
        
        child_expr = child_expr[ , !remp]
        c_names = c_names[!remp]
        dis_child_expr = dis_child_expr[ , !remp]
        
      }
    }
    
    # offset by 1
    dis_parent_expr = dis_parent_expr - 1
    dis_child_expr = dis_child_expr - 1
    
    if(pc_expr_ident){
      
      # build coexpression network with filtration
      coexpnet_result = Discrete_Coexpnet(children = dis_child_expr, 
                                          parents = dis_parent_expr,
                                          c_names = c_names,
                                          p_names = p_names,
                                          nthreads = nthreads,
                                          identical = TRUE)
      
    } else {
      
      # build coexpression network without filtration
      coexpnet_result = Discrete_Coexpnet(children = dis_child_expr, 
                                          parents = dis_parent_expr,
                                          c_names = c_names,
                                          p_names = p_names,
                                          nthreads = nthreads,
                                          identical = FALSE)
      
    }
    
  } else {
    
    # calculate maximum possible cluster
    max_possible = ceiling(sqrt((nrow(parent_expr)/5)))
    
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
    
    if(pc_expr_ident){
      
      # build coexpression network with filtration
      coexpnet_result = Coexpnet(children = child_expr,
                                 parents = parent_expr,
                                 c_names = c_names,
                                 p_names = p_names,
                                 k = k,
                                 nthreads = nthreads,
                                 identical = TRUE)
      
    } else {
      
      # build coexpression network without filtration
      coexpnet_result = Coexpnet(children = child_expr,
                                 parents = parent_expr,
                                 c_names = c_names,
                                 p_names = p_names,
                                 k = k,
                                 nthreads = nthreads,
                                 identical = FALSE)
      
    }
    
  }
  
  # remove NA rows
  rem_entries = is.na(coexpnet_result$PID)
  
  if(any(rem_entries)){
    coexpnet_result = coexpnet_result[!rem_entries, ]
  }
  
  # adjust p-values
  coexpnet_result$PVALUE = p.adjust(coexpnet_result$PVALUE, method = "BH")
  
  # return coexpression network
  return(coexpnet_result)
}

