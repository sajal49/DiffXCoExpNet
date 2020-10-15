// Created by Sajal Kumar
// Copyright (c) NMSU Song lab

#include "methods.h"
#include "Clusters.h"
#include "Joint_Grid.h"

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

// Container filler for stat_collector structure
void fill_stat_collector(stat_collector & sc, int pindex, int cindex, const std::string & pname,
                         const std::string & cname, ldouble pvalue, double estimate){

    sc.pindex = pindex;
    sc.cindex = cindex;
    sc.pname = pname;
    sc.cname = cname;
    sc.pvalue = pvalue;
    sc.estimate = estimate;

}

// Use Kmeans on a n-dimensional data to get cluster assignments
vec<int> GetClusterAssignments(const arma::mat & data, int cluster){

    // initialize a kmeans object (with 100 epochs)
    mlpack::kmeans::KMeans<> kobj(100);

    // prepare a container for cluster labels
    arma::Row<size_t> assignments;

    // apply kmeans
    kobj.Cluster(data, (size_t)cluster, assignments);

    // return the assignments in correct format
    return arma::conv_to< vector<int> >::from(assignments);
}

void create_pair_data(frame<double> & container, const vec<double> & dim1, const vec<double> & dim2){

    // fill container row-wise
    container.clear();
    vec<double> temp_vec;
    for(int i=0; i<dim1.size(); i++){

        temp_vec.push_back(dim1[i]);
        temp_vec.push_back(dim2[i]);
        container.push_back(temp_vec);
        temp_vec.clear();
    }

}

// Apply GridOnCluster on _data using the cluster assignments
frame<int> ApplyGOC(const vec<int> & assignments, const frame<double> & _data, int cluster, int min_levels){

    // prepare clusters
    Cluster clust_obj(cluster, assignments, _data);

    // get grid lines (with minimum levels 2)
    frame<double> grid_lines = Find_Grid(clust_obj, min_levels);

    // discretize each variable
    frame<int> Ddata = DiscretizeData(grid_lines, _data);

    return Ddata;

}

// Discretize _data given grid_lines
frame<int> DiscretizeData(const frame<double> & grid_lines, const frame<double> & _data){

    // discretize each dimension using grid_lines
    frame<int> Ddata(_data[0].size(), vec<int>(_data.size(), 0));
    double start1, end1;
    bool _terminate;

    // temporary vector for gridlines
    vec<double> temp_gridlines;

    for(int j=0; j<_data[0].size(); j++){

        temp_gridlines = grid_lines[j];

        std::sort(temp_gridlines.begin(), temp_gridlines.end());
        start1 = std::numeric_limits<double>::min();

        _terminate = false;

        for(int gl=0; gl<temp_gridlines.size(); gl++){

            end1 = temp_gridlines[gl];

            if(_terminate){ // if _terminate is true, no need to look further
                break;
            }

            if(end1 == std::numeric_limits<double>::max()){
                // no need to proceed further, except the current instance
                _terminate = true;
            }

            for(int x=0; x<_data.size(); x++){
                if(_data[x][j] >= start1 && _data[x][j] <= end1){
                    Ddata[j][x] = gl;
                }
            }
            start1 = end1;
        } // grid_lines
    } // dimension

    return Ddata;
}

// Standardize data (scale and shift)
void StandardScalerTr(frame<double> & data){

    double mean = 0;
    double sd = 0;
    frame<double> std_data(data.size(), vec<double>(data[0].size(), 0));

    for(int i=0; i<data.size(); i++){

        // get mean
        for(int j=0; j<data[i].size(); j++){
            mean += data[i][j];
        }
        mean = ( mean / (double) data[i].size());

        // get sd (sample)
        for(int k=0; k<data[i].size(); k++){
            sd += pow((data[i][k] - mean), 2);
        }
        sd = sqrt(sd)/(sqrt(data[i].size()-1));

        // standardize
        for(int l=0; l<data[i].size(); l++){
            data[i][l] = (data[i][l] - mean)/sd;
        }
    }
}

// Get mean silhouette score
double GetSilMeanScore(const frame<double> sil_scores){

    double sil_mean = 0;

    for(int i=0; i<sil_scores.size(); i++){
        sil_mean += sil_scores[i][4];
    }

    return (sil_mean / sil_scores.size());
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