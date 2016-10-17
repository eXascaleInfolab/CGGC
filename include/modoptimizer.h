//============================================================================
// Name        : RGModularityClusterer.h
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : implements the randomized greedy algorithm and the CGGC scheme
//============================================================================


#ifndef MODOPTIMIZER_H_
#define MODOPTIMIZER_H_

#include "basetypes.h"


class Partition;
class Graph;
class ActiveRowSet;
class SparseClusteringMatrix;

class ModOptimizer {
public:
    ModOptimizer(Graph* graph);
    virtual ~ModOptimizer();

    ModOptimizer(const ModOptimizer&)=delete;
    ModOptimizer& operator =(const ModOptimizer&)=delete;

    Partition* get_clusters();

    double ClusterRG(size_t sample_size, uint16_t runs);
    void ClusterCGGC(size_t ensemble_size, size_t sample_size_restart, bool iterative);
    double GetModularityFromClustering(Graph* graph, Partition* clusters);

private:
    Graph* graph_;
    ActiveRowSet* active_rows_;
    SparseClusteringMatrix* cluster_matrix_;
    Partition* clusters_;

    Partition* CompareClusters(Graph* graph, Partition* partition1, Partition* partition2);
    double PerformJoins(size_t sample_size);
    Partition* PerformJoinsRestart(Graph* graph, Partition* partition
        , size_t sample_size_restart);
    Partition* RefineCluster(Graph* graph, Partition* clusters);
    Partition* GetPartitionFromJoins(t_idpair_vector joins,
        const size_t best_step, Partition* partition);
    t_id_vector* GetMembershipFromPartition(Partition* partition, size_t vertex_count);
};

#endif /* MODOPTIMIZER_H_ */
