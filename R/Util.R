# Util.R
#
# Created by Ruby Sharma
# Modified by Sajal Kumar
# Copyright (c) NMSU Song lab

# Set of utility functions

#' Parallel discretization of each column in 'data'
#'
#' @import doParallel
#' @importFrom parallel parLapply makeForkCluster stopCluster detectCores
#' @importFrom Ckmeans.1d.dp Ckmeans.1d.dp
#' 
#' @param data A samples x genes normalized matrix.
#' 
#' @param ncores Number of parallel threads to use.
#' 
#' @param levls Number of levels for each column in 'data' (optional)
#' 
#' @param n_cond  Number of conditions (groups) samples of 'data' represent (optional)
#' 
#' @return 
#'
#' A matrix containing the discretized \code{data}.
#' 
#' @author 
#' 
#' Ruby Sharma, Sajal Kumar and Mingzhou Song
#' 
#' @seealso
#'
#' See \pkg{Ckmeans.1d.dp} for discretizing univariate continuous data.
#' 
#' @export
DiscretizeData = function(data, 
                          ncores,
                          levls = -1, 
                          n_cond = 1)
{
  
  # check levls
  if(levls != -1){
    len_levls = length(levls)
    if(nrow(data) != length(levls)){
      stop("Length of 'levls' should be equal to ncol(data).")
    }
  }
  
  # if only one core is to be used
  if(ncores==1){
    
    discreet_data = lapply(1:ncol(data), function(x){
      dt = as.numeric(data[,x])
      if(levls == -1){
        undt = length(unique(dt))
        if(undt < 10){
          levl = 2
        } else {
          levl = ceiling((sqrt(undt/5)/n_cond))
          levl = max(2, levl) # cannot be smaller than 2
        }  
      } else {
        levl = levls[x]
      }
      dis = Ckmeans.1d.dp(x=dt,k=c(2,levl))$cluster
      return(dis)
    })
    
    return(discreet_data)
  }
  
  # true number of cores in the machine
  true_cores = detectCores() - 2
  
  # adjust requested cores by the available number of cores
  if(true_cores<ncores)
    ncores = true_cores
  
  # prepare fork clusters
  clust = makeForkCluster(nnodes = ncores)
  
  # discretize data
  discreet_data = parLapply(cl=clust, 1:ncol(data), function(x){
    dt = as.numeric(data[,x])
    if(levls == -1){
      undt = length(unique(dt))
      if(undt < 10){
        levl = 2
      } else {
        levl = ceiling((sqrt(undt/5)/n_cond))
        levl = max(2, levl) # cannot be smaller than 2
      }  
    } else {
      levl = levls[x]
    }
    dis = Ckmeans.1d.dp(x=dt,k=c(2,levl))$cluster
    return(dis)
  })
  
  stopCluster(clust)
  
  # convert list to matrix
  dexpr = ListToDataMatrix(discreet_data)
  
  return(dexpr)
}


# Convert list obtained from DiscretizeData to matrix
# data_list : Result obtained from DiscretizeData
ListToDataMatrix = function(data_list){
  
  new_data = matrix(0,nrow=length(data_list[[1]]),ncol=length(data_list))
  for(i in 1:ncol(new_data))
    new_data[,i] = data_list[[i]]
  return(new_data)
}


# Scale and shift data
# data : A samples x genes normalized matrix
StdScaledTr = function(data){
  
  # get mean
  dmn = apply(data, 2, mean)
  
  # get sd
  dsd = apply(data, 2, sd)
  
  # transform
  for(i in 1:ncol(data)){
    data[,i] = (data[,i] - dmn[i]) / dsd[i]
  }
  
  return(data)
}

#' Returns indices of low variance genes using MAD
#'
#' Employs Median Absolute Deviation (MAD) to return genes expression profiles with zero MAD or their
#' MAD among the bottom 5%. 
#' 
#' @importFrom parallel parSapply makeForkCluster detectCores stopCluster
#'
#' @param data a sample (rows) x gene (columns) continuous expression profile matrix.
#' 
#' @param nthreads an integer specifying the number of parallel threads to be used. If not specified,
#' a non-parallel version of the code would be executed.
#' 
#' @return
#'
#' A logical vector of size ncol(data) containing whether a particular gene show be kept or removed.
#' 
#' @author 
#' 
#' Ruby Sharma, Sajal Kumar and Mingzhou Song
#' 
#' @export
FilterGenesByMAD = function(data, nthreads=-1)
{
  
  # get mad for each gene
  if(nthreads == -1){
    
    mad_collect = sapply(X = c(1:ncol(data)), function(x){
      return(mad(as.numeric(data[,x])))
    })
      
  } else {
    
    # true number of cores in the machine
    true_cores = detectCores() - 2
    
    # adjust requested cores by the available number of cores
    if(true_cores<nthreads)
      nthreads = true_cores
    
    cl = makeForkCluster(nnodes = nthreads)
    
    mad_collect = parSapply(cl = cl, X = c(1:ncol(data)), function(x){
      return(mad(as.numeric(data[,x])))
    })
    
    stopCluster(cl)
    
  }
  
  # remove genes with MAD in the bottom 5% 
  # only considering genes with non-zero MAD
  
  rem_index = mad_collect == 0
  
  nz_mad_collect = sort(unique(mad_collect[!rem_index]))
  cutoff = nz_mad_collect[round(0.05*length(nz_mad_collect))]
  
  rem_index = rem_index | mad_collect <= cutoff
  
  # return indices
  return(rem_index)
  
}

#' Returns common interaction to be considered for differential co-expression network.
#' 
#' Given a list of coexpression networks, this function returns the union of common significant
#' interactions across all co-expression networks, to be considered for differential co-expression
#' network
#' 
#' @param coexpnet_list a list of 2 or more coexpression networks, as returned by \code{Build_Coexpnet()}.
#' 
#' @return
#'
#' A matrix of all gene interactions to study for \code{Build_DiffCoexpnet()}. It containins the indices 
#' (1 based) of parent and child genes.
#' 
#' @author 
#' 
#' Ruby Sharma, Sajal Kumar and Mingzhou Song
#' 
#' @export
GetInteractionIndices = function(coexpnet_list){
  
  # number of experiments
  n_conditions = length(coexpnet_list)
  
  # form interaction pairs
  parents = c()
  for(i in 1:n_conditions){
    parents = c(parents, coexpnet_list[[i]]$PID)
  }
  
  children = c()
  for(i in 1:n_conditions){
    children = c(children, coexpnet_list[[i]]$CID)
  }
  
  interactions = paste0(parents, ",", children)
  interactions = unique(interactions)  
  
  # children and parent sets
  children = as.numeric(gsub(".*,","",interactions)) + 1
  parents = as.numeric(gsub(",.*","",interactions)) + 1
  
  inter_indices = cbind(parents, children)
  
  return(inter_indices)
  
}


#' Plot differential expression across multiple conditions.
#' 
#' Given results from differential coexpression, expression profiles and a vector containing conditions id,
#' this function plots the differential expression across multiple conditions in the same plot.
#' 
#' @importFrom stats loess predict
#' @importFrom ggplot2 alpha
#' @importFrom graphics par layout plot points axis lines legend
#' 
#' @param diffcoexp_res the differential coexpression network, as returned by \code{Build_DiffCoexpnet()}.
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
#' @param g_names a vector containing names of all genes in \code{expr_matr}.
#' 
#' @param clrs a vector of size \code{n_conditions} specifying what colors to use for each experimental 
#' condition
#' 
#' @param cond_cat a vector of size \code{n_conditions} specifying the string category for each experimental
#' condition
#' 
#' @param exp_title the title for the differential experiment
#'
#' 
#' @author 
#' 
#' Ruby Sharma, Sajal Kumar and Mingzhou Song
#' 
#' @export
PlotDiffCoExpNet = function(diffcoexp_res, exp_matr, n_conditions, conditions, g_names, clrs,
                            cond_cat, exp_title){
  
  # unique conditions
  uncn = unique(conditions)
  
  # The length of uncn and n_conditions should be the same
  if(length(uncn) != n_conditions){
    stop("n_conditions and should be equivalent to the number of unique elements in conditions.")
  }
  
  # create a list by condition
  exp_matr_bycond = lapply(1:n_conditions, function(i){
    return(exp_matr[conditions == uncn[i],,drop=FALSE])
  })
  
  # prepare for PlotInter
  
  # prepare categorical conds
  conds = c()
  conds_tab = table(conditions)
  for(i in 1:n_conditions){
    conds = c(conds, rep(cond_cat[i], conds_tab[i]))
  }
  PlotInter(diffcoexp_res, exp_matr, rep(21, n_conditions), clrs, alpha(clrs, 0.5), clrs,
            conds, exp_title, clrs)
}

## Plot top differential interactions 
PlotInter = function(DiffResult, ExpMat, Pointshape, Pointcol, Linecol, Bgcol,
                     Conds, Title, Cols){
  
  par(lwd=2)
  par(mar=c(6,6,6,9)+0.1)
  for(i in 1:nrow(DiffResult))
  {
    UnCond = unique(Conds)
    st = DiffResult[i, ]
    c.id = as.numeric(st[1])+1
    p.id = as.numeric(st[2])+1
    
    
    par(xpd=FALSE)
    layout(matrix(c(1,1,1,2,2,2,2,2,2), nrow = 3, ncol = 3, byrow = TRUE))
    
    par(mar=c(2,2,2,1))
    plot(1:5, 1:5,xaxt="n",yaxt="n",pch="",ylab="",xlab="", main="", sub="", bty= "n")
    
    legend("center", title = Title,
           legend = unique(Conds), pt.bg=Cols, pch= Pointshape,
           col = Pointcol, lty=c(1,1), lwd=rep(10,length(unique(Conds))),
           box.lty=0, cex = 2.5, horiz = TRUE, inset=c(0,3),
           xpd=TRUE)
    
    P = unlist(st[3])
    C = unlist(st[4])
    
    p = as.numeric(ExpMat[, p.id])
    c = as.numeric(ExpMat[, c.id])
    
    par(mar=c(6,6,2.4,1))
    plot(0, type="n", xlim=c(min(p), max(p)), ylim=c(min(c), max(c)+0.2),
         xlab=paste0(P),
         ylab=paste0(C),
         main = paste0("\nP = ",format(as.numeric(st[5]), digits=2), 
                       "\nE = ", format(as.numeric(st[6]), digits=4)), 
         cex.lab=2.5, cex.axis=2, cex.main=2, bty="n", xaxt = "n",yaxt = "n")
    
    
    for(i in 1:length(UnCond)){
      p = as.numeric(ExpMat[Conds==UnCond[i], p.id])
      c = as.numeric(ExpMat[Conds==UnCond[i], c.id])
      PlotPoint(p, c ,Pointshape[i],Bgcol[i],Pointcol[i], Linecol[i])
    }
    
    axis(side=1,   lwd = 4, cex.axis = 1.5, font.axis = 2)
    axis(side=2,   lwd = 4, cex.axis = 1.5, font.axis = 2)
    
  }
}

## Plot points and the loess fit line 
PlotPoint = function(p, c, ptsh, bgcol, ptcol, linecol){
  
  points(c~p, cex = 2.5, pch = ptsh, bg = bgcol, col = ptcol,lwd=4)
  model = predict(loess(jitter(c,amount = 0.05,factor = 0.01)~
                          jitter(p,amount = 0.05,factor = 0.01)),
                  se = TRUE)
  y = model$fit
  x = p
  ord_x = order(x)
  
  y = y[order(x)]
  x = x[order(x)]
  
  lines(y~x,col=linecol, lty=1, lwd=10)
}
