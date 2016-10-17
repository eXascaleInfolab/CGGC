//============================================================================
// Name        : RGModularityClusterer.h
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : implements the randomized greedy algorithm and the CGGC scheme
//============================================================================


#include <vector>
#include <list>

#include <boost/unordered_map.hpp>


#ifndef MODOPTIMIZER_H_
#define MODOPTIMIZER_H_

using namespace std;

class Partition;
class Graph;
class ActiveRowSet;
class SparseClusteringMatrix;

class ModOptimizer {
public:
    ModOptimizer(Graph* graph);
    virtual ~ModOptimizer();

    Partition* get_clusters();

    double ClusterRG(int sample_size, int runs);
    double ClusterCGGC(int ensemble_size, int sample_size_restart,
        bool iterative);
    double GetModularityFromClustering(Graph* graph, Partition* clusters);

private:
    Graph* graph_;
    ActiveRowSet* active_rows_;
    SparseClusteringMatrix* cluster_matrix_;
    Partition* clusters_;

    Partition* CompareClusters(Graph* graph, Partition* partition1,
        Partition* partition2);
    double PerformJoins(int sample_size);
    Partition* PerformJoinsRestart(Graph* graph, Partition* partition,
        int sample_size_restart);
    Partition* RefineCluster(Graph* graph, Partition* clusters);
    Partition* GetPartitionFromJoins(vector<pair<int, int> > joins,
        const int &best_step,  Partition* partition);
    vector<int>* GetMembershipFromPartition(Partition* partition,
        int vertex_count);
};

#endif /* MODOPTIMIZER_H_ */
