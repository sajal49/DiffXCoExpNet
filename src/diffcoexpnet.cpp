#include "methods.h"
#include <boost/thread.hpp>

void discrete_diffcoexpnet_thread(const frame<int> & exp1, const frame<int> & exp2, const vec<int> & parents,
                                  const vec<int> & children, const vec<int> & plevels, const vec<int> & clevels,
                                  const vec<int> & int_index, frame<mydouble> & stat_collect, std::string method)
{
  // to store temporary tables
  std::vector<frame<int > > tables(2);

  // temporary variables for stat collection
  mydouble estimate=0, stat=0, pvalue=1;

  for(int i=0; i<int_index.size(); i++){

    // make two contingency tables

    tables[0] = tableCpp(exp1[parents[int_index[i]]], exp1[children[int_index[i]]], plevels[int_index[i]],
                         clevels[int_index[i]]);
    tables[1] = tableCpp(exp2[parents[int_index[i]]], exp2[children[int_index[i]]], plevels[int_index[i]],
                         clevels[int_index[i]]);

    if(method == "conserved"){
      stat = ConservedTest(tables, pvalue, estimate);
    } else {
      stat = SharmaSongTest(tables, pvalue, estimate);
    }

    stat_collect[0][int_index[i]] = children[int_index[i]];
    stat_collect[1][int_index[i]] = parents[int_index[i]];
    stat_collect[2][int_index[i]] = pvalue;
    stat_collect[3][int_index[i]] = estimate;

  }
}

frame<mydouble> discrete_diffcoexpnet(const frame<int> & exp1, const frame<int> & exp2, const vec<int> & parents,
                                      const vec<int> & children, const vec<int> & plevels, const vec<int> & clevels,
                                      int nthreads, std::string method)
{
  // fix nthreads
  int true_threads = boost::thread::hardware_concurrency()-1;

  if(nthreads > true_threads) {
    nthreads = true_threads;
  }

  if(nthreads < 2) {
    std::cout <<"WARNING: No multithreading enabled!"<<std::endl;
  }

  // create threads
  boost::thread_group threads;

  // stat collector data structure (4 vectors, each of nparents-nchildren elements)
  frame<mydouble> stat_collect(4, vec<mydouble>(parents.size(), 0));

  // load balancing on number of interactions for multi-threads

  frame<int> int_index;
  vec<int> temp_indices;
  int eq_int = (int)ceil((double)parents.size()/nthreads); // total interaction for each thread

  for(int intid=0; intid!=parents.size(); intid++){
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
                                        boost::cref(exp1),
                                        boost::cref(exp2),
                                        boost::cref(parents),
                                        boost::cref(children),
                                        boost::cref(plevels),
                                        boost::cref(clevels),
                                        boost::cref(int_index[i]),
                                        boost::ref(stat_collect),
                                        method));
    } else { // single thread
      discrete_diffcoexpnet_thread(exp1, exp2, parents, children, plevels, clevels, int_index[i],
                                   stat_collect, method);
    }

  }

  // join threads
  if(nthreads >= 2){
    threads.join_all();
  }

  return stat_collect;
}
