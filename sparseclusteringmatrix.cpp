//============================================================================
// Name        : SparseClusteringMatrix.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : Storing the sparse matrix e, which stores the fractions of 
//               edges connecting a pair of vertices (row and column)
//               matrix implemented only for undirected graphs
//============================================================================


#include "sparseclusteringmatrix.h"

#include "graph.h"
#include "partition.h"

SparseClusteringMatrix::SparseClusteringMatrix(Graph* graph) {
    rows_ = NULL;
    row_sums_ = NULL;

    init(graph);
}

SparseClusteringMatrix::SparseClusteringMatrix(Graph* graph, Partition* clusters) {
    dimension_ = clusters->get_partition_vector()->size();

    int* clustermap = new int[graph->get_vertex_count()]; // maps vertex_id -> cluster_id
    for (int i = 0; i < dimension_; i++) {
        t_id_list* cluster = clusters->get_partition_vector()->at(i);
        int cluster_row = *(cluster->begin()); // cluster will be stored in row of first vertex

        for (t_id_list::iterator iter = cluster->begin(); iter != cluster->end(); ++iter) {
            int vertexid = *iter;
            clustermap[vertexid] = cluster_row;
        }
    }

    rows_ = new t_row_value_map[graph->get_vertex_count()];
    row_sums_ = new double[graph->get_vertex_count()];

    double initvalue = 1.0 / (2 * graph->get_edge_count()); // initial value
                                                          // 1 / (2*|E|)


    // for every neighbor fill field in sparse matrix (== insert hash table )
    for (int i = 0; i < graph->get_vertex_count(); i++) {
        int cluster1 = clustermap[i];

        vector<int>* neighbors = graph->GetNeighbors(i);
        int neighborCount = neighbors->size();
        for (int j = 0; j < neighborCount; j++) {
            int cluster2 = clustermap[neighbors->at(j)];

            if (rows_[cluster1].find(cluster2) != rows_[cluster1].end())
                rows_[cluster1][cluster2] += initvalue;
            else
                rows_[cluster1][cluster2] = initvalue;
        }
    }

    for (int i = 0; i < graph->get_vertex_count(); i++) {
        double sum = 0.0;
        if (rows_[i].size() > 0) {
            for (t_row_value_map::iterator j = rows_[i].begin(); j != rows_[i].end(); ++j) {
                sum += rows_[i][j->first];
            }
        }
        row_sums_[i] = sum;
    }

    delete [] clustermap;
}

SparseClusteringMatrix::~SparseClusteringMatrix() {
    if (rows_ != NULL) delete [] rows_;
    if (row_sums_ != NULL) delete [] row_sums_;
}

void SparseClusteringMatrix::init(Graph* graph) {
    dimension_ = graph->get_vertex_count();

    rows_ = new t_row_value_map[dimension_];
    row_sums_ = new double[dimension_];

    double initvalue = 1.0 / (2 * graph->get_edge_count()); // initial value
                                                          // 1 / (2*|E|)

    // for every neighbor fill field in sparse matrix (== insert hash table )
    for (int i = 0; i < dimension_; i++) {
        vector<int>* neighbors = graph->GetNeighbors(i);
        int neighbor_count = neighbors->size();
        rows_[i].rehash(neighbor_count * 1.1);

        for (int j = 0; j < neighbor_count; j++)
            rows_[i][neighbors->at(j)] = initvalue;
    }

    for (int i = 0; i < dimension_; i++) {
        row_sums_[i] = initvalue * rows_[i].size();
    }
}

t_row_value_map* SparseClusteringMatrix::GetRow(int &rowIndex) {
    return &rows_[rowIndex];
}

double& SparseClusteringMatrix::GetRowSum(int &rowIndex) {
    return row_sums_[rowIndex];
}

int SparseClusteringMatrix::GetRowEntries(int &rowIndex) {
    return rows_[rowIndex].size();
}

double& SparseClusteringMatrix::Get(int &rowIndex, int &columnIndex) {
    t_row_value_map::iterator iter = rows_[rowIndex].find(columnIndex);
    return iter->second;
}

/*
 *  Joins two clusters by adding row b to row a. For better performance, row b should have
 *  less entries than row a.
 */
void SparseClusteringMatrix::JoinCluster(int &a, int &b) {
    // Adjust matrix E
    for (t_row_value_map::iterator iter = rows_[b].begin(); iter != rows_[b].end(); ++iter) {
        int column = iter->first;
        double value = iter->second;

        double new_value = value;
        if (rows_[a].find(column) != rows_[a].end())
            new_value += rows_[a].find(column)->second;

        rows_[a][column] = new_value;
        rows_[column][a] = new_value;

        if (column != b) {
            t_row_value_map::iterator iter = rows_[column].find(b);
            rows_[column].erase(iter);
        }
    }

    rows_[a][a] = rows_[a][a] + rows_[a][b];
    rows_[a].erase(b);

    // Adjust vector A
    row_sums_[a] += row_sums_[b];
    row_sums_[b] = 0;
}
