// Created by Ruby Sharma
// Copyright (c) NMSU Song lab
#include "Silhouette.h"

vector<vector<double>> silhouette_score( const vector<int>& clustering, const vector<vector<double>>& DisMat ){

    int n = clustering.size();
    int k = 0, r = 0, c = 0, cli = 0, clj = 0, ct = 0;
    double avg = 0.0, b = 0.0, a = 0.0, div, icd = 0.0, ocd = 0.0,
            m_dis = MAXFLOAT, minCl, ai = 0.0;


    vector<int> clustid = clustering;
    vector<vector<double>> ClusScore(n, vector<double> (5));


    // remove consecutive (adjacent) duplicates and obtain unique clusters
    sort(clustid.begin(), clustid.end());
    auto uni = unique(clustid.begin(), clustid.end());
    clustid.erase(uni , clustid.end());

    // total number of clusters
    k = clustid.size();
    vector<int> count(k);

    // Get the sample count in each clusters
    for(int j = 0; j<n; j++)
    {
        int cl = clustering[j];
        count[cl]++;
    }

    // Find the neighbour for each clusters and populate ClusScore matrix
    for(int i = 0; i<n ; i++) {
        
       vector<double> bi(k, 0);
       cli = clustering[i];
       ai = 0.0;
       ClusScore[i][0] = cli;

       for (int j = 0; j < n; j++) {
           clj = clustering[j];
           if (i == j) {
               ai = ai + 0.0;

           } else {
               if (i > j) {
                   r = i;
                   c = j;
               } else {
                   r = j;
                   c = i;
               }
               if (cli == clj) {

                   icd = DisMat[r][c];
                   ai = ai+icd;

               } else {
                   ocd = DisMat[r][c];
                   bi[clj] = bi[clj] + ocd;
               }

           }

       }

       // Find the average distance all the clusters from a sample point
       // Obtain neighbour of each sample point
       //for(double lm : bi){
       //    cout<<lm;
       //}
       //cout<<"\n";
       //cout<<"i"<<i<<"\n";
       vector<double> avg(k,0);
        m_dis = MAXFLOAT, minCl = 0.0;
       for(int l = 0; l<k; l++) {

           if (l != cli) {
               avg[l] = bi[l] / count[l];
               //cout<<"avg "<<avg[l]<<"\n";
               //cout<<"c: "<<count[l]<<"\n";
               //cout<<"b: "<<bi[l]<<"\n";

               if (m_dis > avg[l]) {
                      m_dis = avg[l];
                      minCl = l;

                  }
              }
           //cout<<"mds"<<m_dis<<"\n";
           ClusScore[i][1] = minCl;
           ClusScore[i][2] = m_dis;
       }

       ct = count[cli]-1;
       if(ai!=0) {
           ClusScore[i][3] = ai / ct;
       }
       else{
           ClusScore[i][3] =  0;
       }
    }

     // Find the silhouette coefficient for each sample point
     for(int i = 0; i<n ; i++){
         cli = clustering[i]-1;
         b = ClusScore[i][2];
         a = ClusScore[i][3];

         if(b>a){
            div = b;
         }else
            div = a;
         if(b != a && count[cli]>1){
             ClusScore[i][4] = (b-a)/div;
         }else{
             ClusScore[i][4] = 0.0;
         }
     }

     return ClusScore;

}

