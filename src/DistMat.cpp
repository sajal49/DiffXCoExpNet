// Created by Ruby Sharma
// Modified by Sajal Kumar
// Copyright (c) NMSU Song lab

#include "Silhouette.h"
using namespace std;

//template <typename T>
//
//// Find distance between two vectors
//double VecDistance(const vector<T>& Firstpt, const vector<T>& Secondpt){
//    std::vector<double> vecd;
//    double sum =0.0;
//
//    // Finding distance using std transform
//    std::transform(Firstpt.begin(),
//                   Firstpt.end(),
//                   Secondpt.begin(),
//                   std::back_inserter(vecd),
//                   [](T point1, T point2){return pow((point1 - point2),2);});
//
//    sum = std::accumulate(vecd.begin(), vecd.end(), 0.0);
//    return sqrt(sum);
//
//}

template <typename T>
// Find distance between two vectors
double VecDistance(const vector<T>& Firstpt, const vector<T>& Secondpt){

    double ed = 0;
    for(int i=0; i<Firstpt.size(); i++){
        ed += pow((Firstpt[i] - Secondpt[i]),2);
    }
    return sqrt(ed);

}

// Find the distance between samples points in multidimension
vector<vector<double>>  DistMat (const vector<vector<double>> & data){

    int rows = data.size();

    vector<vector<double>> DistanceMat(rows, vector<double> (rows));

    for(size_t i = 0; i<rows; i++){
        for(size_t nr = i; nr<rows; nr++){
            DistanceMat[nr][i] = VecDistance(data[i],data[nr]);
        }
    }
    return   DistanceMat;
}

//// Print matrix
//void printMat(const vector<vector<double>>& t){
//
//    int rows = t.size();
//    int cols = t[0].size();
//    for(int i = 0; i <rows; i++  ){
//        for(int j  = 0; j<cols; j++){
//            cout<<t[i][j]<<" ";
//        }
//
//        printf("\n");
//    }
//}
   
