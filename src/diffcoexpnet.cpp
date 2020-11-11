// Created by Sajal Kumar
// Copyright (c) NMSU Song lab

#include "methods.h"
#include "Silhouette.h"
#include <boost/thread.hpp>
#include <armadillo>

/*
 * Worker -- Parallel discrete differential coexpression network builder
 * Parameters
 *
 * exp_matr: A discrete sample x gene matrix representing
 * discrete levels of gene expression at given samples.
 *
 * exp_conds: Total number of experimental conditions.
 *
 * conds: A vector with the same size as the number of samples in exp_matr, representing the
 * experimental condition id (starting 0) of each sample.
 *
 * indices: A matrix containing the indices of interaction pairs to be extracted from exp_matr. The first column
 * contains parent indices and the second column contains child indices.
 *
 * gnames: A vector of strings containing gene names for all genes in exp_matr
 *
 * glevels: A vector containing the discretization levels for all genes in exp_matr
 *
 * int_index: A vector of thread specific interaction indices, specifying how much work the current slave
 * thread would do.
 *
 * sc: An instance of stat_collector to gather results.
 *
 */
void discrete_diffcoexpnet_thread(const frame<int> & exp_matr,
                                  int exp_conds,
                                  const vec<int> & conds,
                                  const frame<int> & indices,
                                  const vec<std::string> & gnames,
                                  const vec<int> & glevels,
                                  const vec<int> & int_index,
                                  vec<stat_collector> & sc)
{
    // to store temporary tables
    std::vector<frame<int > > tables(exp_conds);

    // separate exp_matr into conditions
    vec<frame<int>> cond_exp_matr = Segment_expr_matr_by_conds(exp_matr, exp_conds, conds);

    // temporary variables for stat collection
    double estimate=0, stat=0;
    ldouble pvalue=1;
    size_t df;

    for(int i=0; i<int_index.size(); i++){

        // construct tables by experimental conditions
        for(int cnd = 0; cnd < exp_conds; cnd++){
            tables[cnd] = tableCpp(cond_exp_matr[cnd][indices[0][int_index[i]]],
                                   cond_exp_matr[cnd][indices[1][int_index[i]]],
                                   glevels[indices[0][int_index[i]]],
                                   glevels[indices[1][int_index[i]]]);
        }

        // Apply SharmaSong on tables
        stat = SharmaSongTest(tables, pvalue, estimate, df);

        // update results
        fill_stat_collector(sc[int_index[i]], indices[0][int_index[i]], indices[1][int_index[i]],
                            gnames[indices[0][int_index[i]]], gnames[indices[1][int_index[i]]], pvalue,
                            estimate);

    }
}

/*
 * Master -- Parallel discrete differential coexpression network builder
 * Parameters
 *
 * exp_matr: A discrete sample x gene matrix representing
 * discrete levels of gene expression at given samples.
 *
 * exp_conds: Total number of experimental conditions.
 *
 * conds: A vector with the same size as the number of samples in exp_matr, representing the
 * experimental condition id (starting 0) of each sample.
 *
 * indices: A matrix containing the indices of interaction pairs to be extracted from exp_matr. The first column
 * contains parent indices and the second column contains child indices.
 *
 * gnames: A vector of strings containing gene names for all genes in exp_matr
 *
 * glevels: A vector containing the discretization levels for all genes in exp_matr
 *
 * nthreads: an integer representing the numbers of thread to be used.
 *
 */
vec<stat_collector> discrete_diffcoexpnet(const frame<int> & exp_matr,
                                          int exp_conds,
                                          const vec<int> & conds,
                                          const frame<int> & indices,
                                          const vec<std::string> & gnames,
                                          const vec<int> & glevels,
                                          int nthreads=-1)
{
    // fix nthreads
    int true_threads = (int)boost::thread::hardware_concurrency()-1;

    if(nthreads > true_threads) {
        nthreads = true_threads;
    }

    // if(nthreads < 2) {
    //     std::cout <<"WARNING: No multithreading enabled!"<<std::endl;
    // }

    // create threads
    boost::thread_group threads;

    // // a vector of stat_collector data structure
    vec<stat_collector> sc((int)indices[0].size());

    // load balancing on number of interactions for multi-threads

    frame<int> int_index;
    vec<int> temp_indices;
    int eq_int = (int)ceil((double)indices[0].size()/nthreads); // total interaction for each thread

    for(int intid=0; intid!=indices[0].size(); intid++){
        if((intid+1) % eq_int != 0){ // prepare interaction indices for each thread
            temp_indices.push_back(intid);
        } else {
            temp_indices.push_back(intid);
            int_index.push_back(temp_indices); // move on to next thread
            temp_indices.clear(); //vec<int>().swap(temp_indices); // clear temp_indices
        }
    }
    if(temp_indices.size()!=0){
        int_index.push_back(temp_indices); // residual interaction indices
    }

    // Create threads for each set of interactions
    for(int i=0; i<int_index.size(); i++){

        if(nthreads >= 2){ // multiple threads
            threads.create_thread(boost::bind(&discrete_diffcoexpnet_thread,
                                              boost::cref(exp_matr),
                                              exp_conds,
                                              boost::cref(conds),
                                              boost::cref(indices),
                                              boost::cref(gnames),
                                              boost::cref(glevels),
                                              boost::cref(int_index[i]),
                                              boost::ref(sc)));
        } else { // single thread
            discrete_diffcoexpnet_thread(exp_matr, exp_conds, conds, indices, gnames,
                                         glevels, int_index[i], sc);
        }

    }

    // join threads
    if(nthreads >= 2){
        threads.join_all();
    }

    return sc;
}

/*
 * Worker -- Parallel differential coexpression network builder that discretizes on the fly
 * Parameters
 *
 * exp_matr: A continuous sample x gene matrix representing
 * discrete levels of gene expression at given samples.
 *
 * exp_conds: Total number of experimental conditions.
 *
 * conds: A vector with the same size as the number of samples in exp_matr, representing the
 * experimental condition id (starting 0) of each sample.
 *
 * indices: A matrix containing the indices of interaction pairs to be extracted from exp_matr. The first column
 * contains parent indices and the second column contains child indices.
 *
 * gnames: A vector of strings containing gene names for all genes in exp_matr
 *
 * int_index: A vector of thread specific interaction indices, specifying how much work the current slave
 * thread would do.
 *
 * k: a vector of integers that can contain the following:
 * i) k with size two: the minimum and maximum number of discrete levels allowed.
 * ii) k with size one: the maximum number of discrete levels allowed (minimum would be assumed to be 2).
 * iii) k with size n (where n is indices[0].size()): the number of discrete levels to be used for each
 * interaction.
 *
 * sc: An instance of stat_collector to gather results.
 *
 */
void diffcoexpnet_thread(const frame<double> & exp_matr,
                         int exp_conds,
                         const vec<int> & conds,
                         const frame<int> & indices,
                         const vec<std::string> & gnames,
                         const vec<int> & k,
                         const vec<int> & int_index,
                         vec<stat_collector> & sc)
{
    // to store temporary tables
    std::vector<frame<int > > tables(exp_conds);

    // temporary variables for stat collection
    double estimate=0, stat=0;
    ldouble pvalue=1;
    size_t df;
    int temp_k_min, temp_k_max, unq;

    // temporary arma mat for all pairs
    arma::mat temp_pair_data(2, exp_matr[0].size());

    // temporary frame to store gene pairs
    frame<double> temp_RData;

    // temporary frame to store discretized data
    frame<int> temp_Ddata;

    // temporary vector of frame to store segmented discretized pair
    vec<frame<int>> temp_cond_Ddata;

    // temporary vector to store cluster assignments
    vec<int> temp_classn;

    // temporary frame to store distance matrix
    frame<double> dist_mat;

    // temporary variables to store silhouette scores (and best k)
    frame<double> sil_scores;
    double mean_sil_score;
    int best_k;
    double best_sil_score;

    // separate exp_matr into conditions
    vec<frame<double>> cond_exp_matr = Segment_expr_matr_by_conds(exp_matr, exp_conds, conds);

    // scale and shift data
    for(int cnd = 0; cnd < exp_conds; cnd++){
        StandardScalerTr(cond_exp_matr[cnd]);
    }

    // merge scaled data
    frame<double> std_exp_matr = Join_expr_matr_by_conds(cond_exp_matr, exp_conds, exp_matr[0].size());

    for(int i=0; i<int_index.size(); i++){
       
        // create data
        create_pair_data(temp_RData, std_exp_matr[indices[0][int_index[i]]],
                         std_exp_matr[indices[1][int_index[i]]]);

        // find the number of unique samples
        unq = FindUniqueSamples(temp_RData);

        // if there is only one unique sample, no need to proceed
        if( unq == 1 ){

            fill_stat_collector(sc[int_index[i]], indices[0][int_index[i]], indices[1][int_index[i]],
                                gnames[indices[0][int_index[i]]], gnames[indices[1][int_index[i]]], 1,
                                0);
        } else {

            for(int mi=0; mi<temp_RData[0].size(); mi++){
                for(int mj=0; mj<temp_RData.size(); mj++){
                    temp_pair_data(mi,mj) = temp_RData[mj][mi];
                }
            }

            // if size of k is 2 then it specifies the maximum and minimum range.
            // if size of k is 1 then it specifies the maximum range
            if(k.size() == 2 || k.size() == 1){

                // get distance matrix
                dist_mat = DistMat(temp_RData);

                // default k
                temp_k_min = 2;
                temp_k_max = 2;

                // determine minimum and maximum k
                if(k.size() == 1){

                    // we will change temp_k_max iff the unique number of points
                    // are larger than 10
                    if( unq > 10 ){
                        if ((int) (unq / 5) < k[0]){
                            temp_k_max = (int) (unq / 5);
                        } else {
                            temp_k_max = k[0];
                        }
                    }

                } else {

                    // we will change temp_k_min and temp_k_max iff the unique number of points
                    // are larger than 10
                    if(unq > 10){
                        if ((int) (unq / 5) < k[0]){
                            temp_k_min = (int) (unq / 5);
                        } else {
                            temp_k_min = k[0];
                        }

                        if ((int) (unq / 5) < k[1]){
                            temp_k_max = (int) (unq / 5);
                        } else {
                            temp_k_max = k[1];
                        }

                    }
                }

                best_k = temp_k_min;
                best_sil_score = (double)(-INT_MAX);

                // find the optimal k using Silhouette score
                for(int m=temp_k_min; m<=temp_k_max; m++){

                    // get cluster assignments
                    temp_classn = GetClusterAssignments(temp_pair_data, m);

                    // get Silhouette score for each observation
                    sil_scores = silhouette_score(temp_classn, dist_mat);

                    // obtain mean Silhouette score
                    mean_sil_score = GetSilMeanScore(sil_scores);

                    // obtain best sil score / best k
                    if(best_sil_score <= mean_sil_score){
                        best_sil_score = mean_sil_score;
                        best_k = m;
                    }
                }

                // finally apply GOC using the best k
                // get cluster assignments
                temp_classn = GetClusterAssignments(temp_pair_data, best_k);

                // get discretized data
                temp_Ddata = ApplyGOC(temp_classn, temp_RData, best_k);

            } else if(k.size() == indices[0].size()){
                // if size of k is the same as the total number of interactions, then it specifies the exact k
                // for each combination

                // get cluster assignments
                temp_classn = GetClusterAssignments(temp_pair_data, k[int_index[i]]);

                // get discretized data
                temp_Ddata = ApplyGOC(temp_classn, temp_RData, k[int_index[i]]);
            }

            // segment discretized data
            temp_cond_Ddata = Segment_expr_matr_by_conds(temp_Ddata, exp_conds, conds);

            // construct tables by experimental conditions
            for(int cnd = 0; cnd < exp_conds; cnd++){
                tables[cnd] = tableCpp(temp_cond_Ddata[cnd][0],
                                       temp_cond_Ddata[cnd][1],
                                       (*std::max_element(temp_Ddata[0].begin(),
                                                          temp_Ddata[0].end()))+1,
                                       (*std::max_element(temp_Ddata[1].begin(),
                                                          temp_Ddata[1].end()))+1);
            }

            // Apply SharmaSong on tables
            stat = SharmaSongTest(tables, pvalue, estimate, df);

            // update results
            fill_stat_collector(sc[int_index[i]], indices[0][int_index[i]], indices[1][int_index[i]],
                                gnames[indices[0][int_index[i]]], gnames[indices[1][int_index[i]]], pvalue,
                                estimate);
        }

    }
}

/*
 * Master -- Parallel differential coexpression network builder that discretizes on the fly.
 * Parameters
 *
 * exp_matr: A continuous sample x gene matrix representing
 * gene expression at given samples.
 *
 * exp_conds: Total number of experimental conditions.
 *
 * conds: A vector with the same size as the number of samples in exp_matr, representing the
 * experimental condition id (starting 0) of each sample.
 *
 * indices: A matrix containing the indices of interaction pairs to be extracted from exp_matr. The first column
 * contains parent indices and the second column contains child indices.
 *
 * gnames: A vector of strings containing gene names for all genes in exp_matr
 *
 * k: a vector of integers that can contain the following:
 * i) k with size two: the minimum and maximum number of discrete levels allowed.
 * ii) k with size one: the maximum number of discrete levels allowed (minimum would be assumed to be 2).
 * iii) k with size n (where n is indices[0].size()): the number of discrete levels to be used for each
 * interaction.
 *
 * nthreads: an integer representing the numbers of thread to be used.
 *
 */
vec<stat_collector> diffcoexpnet(const frame<double> & exp_matr,
                                 int exp_conds,
                                 const vec<int> & conds,
                                 const frame<int> & indices,
                                 const vec<std::string> & gnames,
                                 const vec<int> & k,
                                 int nthreads=-1)
{
    // fix nthreads
    int true_threads = (int)boost::thread::hardware_concurrency()-1;

    if(nthreads > true_threads) {
        nthreads = true_threads;
    }

    // if(nthreads < 2) {
    //     std::cout <<"WARNING: No multithreading enabled!"<<std::endl;
    // }

    // create threads
    boost::thread_group threads;

    // // a vector of stat_collector data structure
    vec<stat_collector> sc((int)indices[0].size());

    // load balancing on number of interactions for multi-threads

    frame<int> int_index;
    vec<int> temp_indices;
    int eq_int = (int)ceil((double)indices[0].size()/nthreads); // total interaction for each thread

    for(int intid=0; intid!=indices[0].size(); intid++){
        if((intid+1) % eq_int != 0){ // prepare interaction indices for each thread
            temp_indices.push_back(intid);
        } else {
            temp_indices.push_back(intid);
            int_index.push_back(temp_indices); // move on to next thread
            temp_indices.clear(); //vec<int>().swap(temp_indices); // clear temp_indices
        }
    }
    if(temp_indices.size()!=0){
        int_index.push_back(temp_indices); // residual interaction indices
    }

    // Create threads for each set of interactions
    for(int i=0; i<int_index.size(); i++){

        if(nthreads >= 2){ // multiple threads
            threads.create_thread(boost::bind(&diffcoexpnet_thread,
                                              boost::cref(exp_matr),
                                              exp_conds,
                                              boost::cref(conds),
                                              boost::cref(indices),
                                              boost::cref(gnames),
                                              boost::cref(k),
                                              boost::cref(int_index[i]),
                                              boost::ref(sc)));
        } else { // single thread
            diffcoexpnet_thread(exp_matr, exp_conds, conds, indices, gnames,
                                k, int_index[i], sc);
        }

    }

    // join threads
    if(nthreads >= 2){
        threads.join_all();
    }

    return sc;
}
