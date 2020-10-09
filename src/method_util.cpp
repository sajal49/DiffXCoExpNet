#include "methods.h"

frame<int> PoolTables(const vec<frame<int > > & obs){

  // pooled table
  frame<int> pooled_table(obs[0].size(), std::vector<int>(obs[0][0].size(), 0));

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

frame<ldouble> pchisq_expec(const frame<int> & obs, const vec<int> & rowsums, const vec<int> & colsums,
                             int totsum){

  // get dims
  int rows = obs.size();
  int cols = obs[0].size();

  // get expected
  frame<ldouble> expec(rows, vec<ldouble>(cols, 0.0));

  for(size_t i=0; i<rows; i++){
    for(size_t j=0; j<cols; j++){
      expec[i][j] = ((ldouble)rowsums[i]*(ldouble)colsums[j])/totsum;
    }
  }

  return expec;
}

frame<ldouble> fchisq_expec(const frame<int> & obs, const vec<int> & rowsums){

  // get dims
  int rows = obs.size();
  int cols = obs[0].size();

  // get expected
  frame<ldouble> expec(rows, vec<ldouble>(cols, 0.0));

  for(size_t i=0; i<rows; i++){
    for(size_t j=0; j<cols; j++){
      expec[i][j] = ((ldouble)rowsums[i]/cols);
    }
  }

  return expec;
}

frame<ldouble> cpchisq_expec(const frame<int> & pooled_obs, const vec<int> & pooled_rowsums,
                              const vec<int> & pooled_colsums, int pooled_sum, int t_sum){

  // get dims
  int rows = pooled_obs.size();
  int cols = pooled_obs[0].size();

  // get expected
  frame<ldouble> expec(rows, vec<ldouble>(cols, 0.0));

  for(size_t i=0; i<rows; i++){
    for(size_t j=0; j<cols; j++){
      expec[i][j] = ((ldouble)pooled_rowsums[i]*(ldouble)pooled_colsums[j]*(ldouble)t_sum)/pooled_sum;
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

frame<int> tableCpp(std::vector<int> x_vec, std::vector<int> y_vec, int xlevels, int ylevels) {

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

  std::vector<std::vector<int > > table(xlevels, std::vector<int>(ylevels,0));

  for(int i=0; i<x_vec.size(); i++)
    table[x_vec[i]][y_vec[i]]++;

  return table;
}