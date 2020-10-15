
#include <cmath>
#include <iterator>
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <functional>
#include <cctype>
#include <numeric>
#include <fstream>
#include <sstream>





using namespace std;

// Function to find the distance between two vectors.
template <typename T>
double VecDistance(const vector<T>& Firstpt,const vector<T>& Secondpt);

// Function to create a distance matrix 
vector<vector<double>>  DistMat (const vector<vector<double>> & data);

// Function to print the matrix
void printMat(const vector<vector<double>>& t);

vector<vector<double>> silhouette_score(const vector<int>& clustering, const vector<vector<double>> & DisMat);