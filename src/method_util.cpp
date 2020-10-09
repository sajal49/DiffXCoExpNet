#include "methods.h"

frame<int> PoolTables(const vec<frame<int > > & obs){

  // pooled table
  frame<int> pooled_table(obs[0].size(), vector<int>(obs[0][0].size(), 0));

  // cell-wise addition of all observed tables
  for(int i=0; i<obs.size(); i++){ //experimental conditions
    for(int j=0; j<obs[0].size(); j++){ // rows in each table
      for(int k=0; k<obs[0][0].size(); k++){ // columns in each table
        pooled_table[j][k] += obs[i][j][k];
      }
    }
  }

  return pooled_table;
}

void DimSums(const frame<int> & O, vec<int> & rowsums, vec<int> & colsums, int & totsum){

  // row, col and total sum
  totsum = 0;

  // clean rowsums
  for(size_t i=0; i<rowsums.size(); i++){
    rowsums[i] = 0;
  }

  // clean colsums
  for(size_t i=0; i<colsums.size(); i++){
    colsums[i] = 0;
  }

  for (size_t i=0; i<O.size(); i++) {
    for (size_t j=0; j<O[i].size(); j++) {
      totsum += O[i][j];
      colsums[j] += O[i][j];
      rowsums[i] += O[i][j];
    }
  }

  return;
}

frame<double> pchisq_expec(const frame<int> & obs, const vec<int> & rowsums, const vec<int> & colsums,
                           int totsum){

  // get dims
  int rows = obs.size();
  int cols = obs[0].size();

  // get expected
  frame<double> expec(rows, vec<double>(cols, 0.0));

  for(size_t i=0; i<rows; i++){
    for(size_t j=0; j<cols; j++){
      expec[i][j] = ((double)rowsums[i]*(double)colsums[j])/totsum;
    }
  }

  return expec;
}


void fill_stat_collector(stat_collector & sc, int pindex, int cindex, std::string pname, 
                         std::string cname, ldouble pvalue, ldouble estimate){
  
  sc.pindex = pindex;
  sc.cindex = cindex;
  sc.pname = pname;
  sc.cname = cname;
  sc.pvalue = pvalue;
  sc.estimate = estimate;
  
}

vec<int> GetClusterAssignments(const vector<vector<double>> & _data, int cluster, int epoch=-1){
  
  // if epoch == -1 then use default epoch
  if(epoch == -1){
    KMeans<> kobj();  
  } else { // otherwise specify the number of epoch
    KMeans<> kobj(epoch);
  }
  
  // create an arma matrix from _data
  arma::mat data(_data);
  
  // prepare a container for cluster labels
  arma::Row<size_t> assignments;
  
  // apply kmeans
  kobj.Cluster(data, (size_t)cluster, assignments);
  
  // return the assignments in correct format
  return arma::conv_to<vector<int>>>::from(assignments);
}


frame<int> ApplyGOC(const vec<int> & assignments, const frame<double> & _data, int cluster, int min_levels=2){
  
  
  // prepare clusters
  Cluster clust_obj(cluster, assignments, _data);
  
  // get grid lines (with minimum levels 2)
  frame<double> grid_lines = Find_Grid(clust_obj, min_levels);
  
  // discretize each variable
  
}


discretize_data = function(data, gridlines){
  
# use gridlines for each dimension
  for(i in 1:ncol(data)){
    
    if(length(unique(gridlines))==0){
      discr = rep(1, nrow(data))
    }else{
      discr = rep(length(gridlines[[i]])+1, nrow(data))
      gridlines[[i]] = sort(gridlines[[i]])
      for(j in 1:length(gridlines[[i]])){ # determine discretization levels
        if(j == 1) {
          discr[data[,i] < gridlines[[i]][j]] = 1
        } else {
          discr[data[,i] < gridlines[[i]][j] & data[,i] >= gridlines[[i]][j-1]] = j
        }
      }
    }
    data[,i] = discr
  }
  
  return(data)
}

frame<int> DiscretizeData(const frame<double> & grid_lines, const frame<double> & _data){
  
  // discretize each dimension using grid_lines
  frame<int> _Ddata(_data.size(), vec<int>(_data[0].size(), 0));
  vec<int> temp_vec;
  for(int j=0; j<_data.size(); j++){
    
    if(std::unique(grid_lines[j]).size() > 0){
      
      std::sort(grid_lines[j].beg(), grid_lines[j].end());
      
      for(int )
      
    }
    
  }
  
}


frame<int> tableCpp(vector<int> x_vec, vector<int> y_vec, int xlevels, int ylevels) {

  int min_val_x = INT_MAX, min_val_y = INT_MAX, max_val_x = INT_MIN, max_val_y = INT_MIN;

  // Find min and max for both vectors
  for(int i=0; i<x_vec.size(); i++){
    if(x_vec[i] < min_val_x)
      min_val_x = x_vec[i];

    if(x_vec[i] > max_val_x)
      max_val_x = x_vec[i];

    if(y_vec[i] < min_val_y)
      min_val_y = y_vec[i];

    if(y_vec[i] > max_val_y)
      max_val_y = y_vec[i];
  }


  //If xlevels and ylevels are not provided and the smallest element is greater than or equal to 1
  // then subtract min_val to offset all elements to 0.

  if(xlevels == -1 && min_val_x >= 1)
  {
    for(int i=0; i<x_vec.size(); i++){
      x_vec[i]-=min_val_x;
    }
    max_val_x-=min_val_x;
    min_val_x = 0;
  }

  if(ylevels == -1 && min_val_y >= 1)
  {
    for(int i=0; i<y_vec.size(); i++){
      y_vec[i]-=min_val_y;
    }
    max_val_y-=min_val_y;
    min_val_y=0;
  }

  // setup xlevels and ylevels if not provided
  if(xlevels == -1 && ylevels == -1){
    xlevels = max_val_x+1;
    ylevels = max_val_y+1;
  } else if(xlevels == -1){
    xlevels = max_val_x+1;
  } else if(ylevels == -1){
    ylevels = max_val_y+1;
  }

  // ASSUMPTION : tables are already offset by 1 if xlevels and ylevels are provided

  vector<vector<int > > table(xlevels, vector<int>(ylevels,0));

  for(int i=0; i<x_vec.size(); i++)
    table[x_vec[i]][y_vec[i]]++;

  return table;
}