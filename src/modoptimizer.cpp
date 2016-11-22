//============================================================================
// Name        : RGModularityClusterer.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : implements the randomized greedy algorithm and the CGGC scheme
//============================================================================


#include <cstdlib>  // atoi, strtoul
#ifdef DEBUG
#include <cassert>
#endif // DEBUG

#include "sparseclusteringmatrix.h"
#include "activerowset.h"
#include "graph.h"
#include "partition.h"
#include "modoptimizer.h"

using std::make_pair;


ModOptimizer::ModOptimizer(Graph* graph): graph_(graph), active_rows_(nullptr)
    , cluster_matrix_(nullptr), clusters_(nullptr)  {}

ModOptimizer::~ModOptimizer() {
    delete clusters_;
}

Partition* ModOptimizer::get_clusters() {
    return clusters_;
}

double ModOptimizer::ClusterRG(size_t k, uint16_t runs) {
    Partition* best_partition = nullptr;
    double best_q = -1;

    for (uint16_t i = 0; i < runs; i++) {
        double Q = PerformJoins(k);
        if (Q > best_q) {
            if (best_q != -1)
                delete best_partition;

            best_q= Q;
            best_partition = clusters_;
        }
        else
            delete clusters_;
    }

    clusters_ = RefineCluster(graph_, best_partition);
    delete best_partition;
    return best_q;
}

void ModOptimizer::ClusterCGGC(size_t initclusters, size_t restartk, bool iterative) {
    Partition* clusterings[initclusters];
    Partition* unjoined_clusters[initclusters];

    ClusterRG(1, 1);
    unjoined_clusters[0] = clusters_;
    for (size_t i = 1; i < initclusters; i++) {
        ClusterRG(1, 1);
        clusterings[i] = clusters_;
        unjoined_clusters[i] = CompareClusters(graph_, unjoined_clusters[i - 1],
                                               clusterings[i]);

        delete clusterings[i];
        delete unjoined_clusters[i - 1];
    }

    Partition* tmp_clustering = unjoined_clusters[initclusters - 1];

    if (iterative) {
        double cur_q = GetModularityFromClustering(graph_, tmp_clustering);
        double last_q = 0;

        while ((cur_q - last_q) > 0.0001) {
            unjoined_clusters[0] = PerformJoinsRestart(graph_, tmp_clustering, 1);
            for (size_t i = 1; i < initclusters; i++) {
                clusterings[i] = PerformJoinsRestart(graph_, tmp_clustering, 1);
                unjoined_clusters[i] = CompareClusters(graph_,
                        unjoined_clusters[i - 1], clusterings[i]);

                delete clusterings[i];
                delete unjoined_clusters[i - 1];
            }
            last_q = cur_q;
            cur_q = GetModularityFromClustering(graph_,
                    unjoined_clusters[initclusters - 1]);


            if (cur_q > last_q) {
                delete tmp_clustering;
                tmp_clustering = unjoined_clusters[initclusters - 1];
            } else
                delete unjoined_clusters[initclusters - 1];
        }
    }

    Partition* joinrestartclusters = PerformJoinsRestart(graph_, tmp_clustering,
                                                         restartk);
    delete tmp_clustering;
    Partition* result = RefineCluster(graph_, joinrestartclusters);
    delete joinrestartclusters;
    clusters_ = result;
    //return GetModularityFromClustering(graph_, clusters_);
}

t_id_vector* ModOptimizer::GetMembershipFromPartition(Partition* partition
, size_t vertex_count) {
    t_id_vector* membership = new t_id_vector(vertex_count);
    for (size_t i = 0; i < partition->get_partition_vector()->size(); i++) {
        t_id_list* cluster = partition->get_partition_vector()->at(i);
        for (t_id vertex_id: *cluster) {
            membership->at(vertex_id) = i;
        }
    }
    return membership;
}

Partition* ModOptimizer::CompareClusters(Graph* graph, Partition* partition1
, Partition* partition2) {
    vector<t_id_vector*> maps;
    maps.push_back(GetMembershipFromPartition(partition1,
                                              graph->get_vertex_count()));
    maps.push_back(GetMembershipFromPartition(partition2,
                                              graph->get_vertex_count()));

    Partition* result_clustering = new Partition();
    vector<bool> assigned(graph->get_vertex_count(), false);

    t_id_list* newcluster = nullptr;
    for (size_t i = 0; i < partition1->get_partition_vector()->size(); i++) {
        t_id_list* cluster = partition1->get_partition_vector()->at(i);
        for (t_id vertex1: *cluster) {
            if (!assigned[vertex1]) {
                newcluster = new t_id_list();
                result_clustering->get_partition_vector()->push_back(newcluster);
                newcluster->push_back(vertex1);
                assigned[vertex1] = true;
            } else
                continue;

            for (t_id vertex2: *cluster) {
                if (!assigned[vertex2]) {
                    if (maps[1]->at(vertex1) == maps[1]->at(vertex2)) {
                        newcluster->push_back(vertex2);
                        assigned[vertex2] = true;
                    }
                }
            }
        }
    }

    return result_clustering;
}

double ModOptimizer::PerformJoins(size_t sample_size) {
    ActiveRowSet active_rows(graph_->get_vertex_count());
    SparseClusteringMatrix cluster_matrix(graph_);

    size_t dimension = graph_->get_vertex_count();

    // Ensure that dimension is positive
    if(!dimension)
        return -1;
    t_idpair_vector joins(dimension - 1);
    size_t best_step = -1;  // Max value that means not initialized
    double best_step_q = -1;  // Note: evaluated Q always > -0.5 > best_step_q

    //**********
    // calc initial Q
    //**********
    double Q = 0;
    for (size_t i = 0; i < dimension; i++) {
        double a_i = cluster_matrix.GetRowSum(i);
        Q -= a_i * a_i;
    }

    //**********
    // perform joins
    //**********

    for (size_t step = 0; step < graph_->get_vertex_count() - 1; step++) {

        size_t max_sample = sample_size;
        if (sample_size < graph_->get_vertex_count() / 2) {
            max_sample = 1;
        } else if (sample_size < (graph_->get_vertex_count() - 1 - step)) {
            max_sample = sample_size;
        } else {
            max_sample = graph_->get_vertex_count() - 1 - step;
        }


        // *******
        // find join
        // *******

        double max_delta_q = -1;
        t_index join_a = -1; //  the two clusters to join
        t_index join_b = -1;

        for (size_t sample_num = 0; sample_num < max_sample; sample_num++) {

            t_index row_num = -1;
            if (max_sample == graph_->get_vertex_count() - 1 - step)
                row_num = active_rows.Get(sample_num);
            else
                row_num = active_rows.GetRandomElement();

            t_row_value_map* sample_row = cluster_matrix.GetRow(row_num);

            for (t_row_value_map::iterator entry = sample_row->begin(); entry != sample_row->end(); ++entry) {
                t_index column_num = entry->first;
                t_value value = entry->second;

                if (column_num == row_num) continue;

                double delta_q = 2 * (value - cluster_matrix.GetRowSum(row_num)*
                        cluster_matrix.GetRowSum(column_num));

                if (delta_q > max_delta_q) {
                    max_delta_q = delta_q;
                    if (cluster_matrix.GetRowEntries(row_num) >=
                            cluster_matrix.GetRowEntries(column_num)) {
                        join_a = row_num;
                        join_b = column_num;
                    } else {
                        join_a = column_num;
                        join_b = row_num;
                    }
                }
            }

            if (sample_num == max_sample - 1 && max_delta_q < 0
            && max_sample < graph_->get_vertex_count() - 1 - step) {
                if(max_sample < graph_->get_vertex_count()/2) {
                    max_sample++;
                } else {
                    max_sample = graph_->get_vertex_count() - 1 - step;
                    sample_num = 1;
                }
            }
        }

        // if there is no valid merge, stop merge process
        // (can only occur for unconnected graph)
        if (join_a == static_cast<t_index>(-1)) break;

        // *******
        // execute join
        // *******
        cluster_matrix.JoinCluster(join_a, join_b);
        active_rows.Remove(join_b);
        joins[step].first = join_a;
        joins[step].second = join_b;
        Q += max_delta_q;

        if (Q > best_step_q) {
            best_step_q = Q;
            best_step = step;
        }
    }

    clusters_ = GetPartitionFromJoins(joins, best_step, nullptr);
#ifdef DEBUG
    assert(best_step_q > -0.5 && "PerformJoins(), modularity should E (-0.5, 1]");
#endif // DEBUG
    return best_step_q;
}

Partition* ModOptimizer::PerformJoinsRestart(Graph* graph, Partition* clusters
, size_t k_restart_) {
    SparseClusteringMatrix cluster_matrix(graph, clusters);
    ActiveRowSet active_rows(clusters);
    size_t dimension = clusters->get_partition_vector()->size();

    // Ensure that dimension is positive
    if(!dimension)
        return nullptr;
    t_idpair_vector joins(dimension - 1);

    size_t best_step = -1;  // Max value that means not initialized
    double best_step_q = -1;  // Note: evaluated Q always > -0.5 > best_step_q

    double modularity = 0; // not the actual start value of Q,

    //**********
    // perform joins
    //**********
    for (size_t step = 0; step < clusters->get_partition_vector()->size() - 1; step++) {

        size_t max_sample = k_restart_;
        if (k_restart_ < (clusters->get_partition_vector()->size() - 1 - step)) {
            max_sample = k_restart_;
        } else {
            max_sample = clusters->get_partition_vector()->size() - 1 - step;
        }

        // *******
        // find join
        // *******
        double max_delta_q = -1;
        t_index join_a = -1; //  the two clusters to join
        t_index join_b = -1;

        for (size_t sample_num = 0; sample_num < max_sample; sample_num++) {
            t_index row_num;
            if (max_sample == clusters->get_partition_vector()->size() - 1 - step) {
                row_num = active_rows.Get(sample_num);
            } else
                row_num = active_rows.GetRandomElement();

            t_row_value_map* sample_row = cluster_matrix.GetRow(row_num);

            for (t_row_value_map::iterator entry = sample_row->begin();
                    entry != sample_row->end(); ++entry) {
                t_index column_num = entry->first;
                t_value value = entry->second;

                if (column_num == row_num) continue;

                double delta_q = 2 * (value - cluster_matrix.GetRowSum(row_num)
                        * cluster_matrix.GetRowSum(column_num));
                if (delta_q > max_delta_q) {
                    max_delta_q = delta_q;
                    if (cluster_matrix.GetRowEntries(row_num) >=
                            cluster_matrix.GetRowEntries(column_num)) {
                        join_a = row_num;
                        join_b = column_num;
                    } else {
                        join_a = column_num;
                        join_b = row_num;
                    }
                }
            }
            if (sample_num == max_sample - 1 && max_delta_q < 0
            && max_sample < dimension - 1 - step)
                max_sample++;
        }

        // if there is no valid merge, stop merge process
        // (can only occur for unconnected graph)
        if (join_a == static_cast<t_index>(-1)) break;

        // *******
        // execute join
        // *******
        cluster_matrix.JoinCluster(join_a, join_b);
        active_rows.Remove(join_b);
        joins[step].first = join_a;
        joins[step].second = join_b;
        modularity += max_delta_q;

        if (modularity > best_step_q) {
            best_step_q = modularity;
            best_step = step;
        }
    }
#ifdef DEBUG
    assert(best_step_q > -0.5 && "PerformJoinsRestart(), modularity should E (-0.5, 1]");
#endif // DEBUG
    Partition* new_clusters = GetPartitionFromJoins(joins, best_step, clusters);
    return new_clusters;
}

Partition* ModOptimizer::GetPartitionFromJoins(t_idpair_vector joins
, const size_t bestStep, Partition* partial_partition) {
    Partition* result_partition = new Partition();
    if (partial_partition == nullptr) { // create new singleton partition
        // Initialize clusters
        for (size_t i = 0; i < graph_->get_vertex_count(); i++) {
            t_id_list* vlist = new t_id_list();
            vlist->push_back(i);
            result_partition->get_partition_vector()->push_back(vlist);
        }
    } else { // rearrange input partition
        // we need to create the complete list
        for (size_t i = 0; i < graph_->get_vertex_count(); i++) {
            t_id_list* vlist = new t_id_list();
            result_partition->get_partition_vector()->push_back(vlist);
        }

        for (size_t i = 0; i < partial_partition->get_partition_vector()->size(); i++) {
            // the first element of the list determines where to put the list
            t_id pos = *(partial_partition->get_partition_vector()->at(i)->begin());
            for (t_id vertex:
                    *(partial_partition->get_partition_vector()->at(i))) {
                result_partition->get_partition_vector()->at(pos)->
                        push_back(vertex);
            }
        }
    }

    //join clusters according to join list
    for (size_t step = 0; step <= bestStep; step++) {
        t_id_list* list1 =
                result_partition->get_partition_vector()->at(joins[step].first);
        t_id_list* list2 =
                result_partition->get_partition_vector()->at(joins[step].second);

        list1->splice(list1->end(), *list2);
        delete result_partition->get_partition_vector()->at(joins[step].second);
        result_partition->get_partition_vector()->at(joins[step].second) = nullptr;
    }

    result_partition->RemoveEmptyEntries();
    return result_partition;
}

Partition* ModOptimizer::RefineCluster(Graph* graph, Partition* clusters) {
    clusters->RemoveEmptyEntries();

    size_t cluster_count = clusters->get_partition_vector()->size();
    t_id_vector clusterdegree(cluster_count); // sum of degrees of all vertices of a cluster
    t_id_vector clustermap(graph->get_vertex_count()); // maps vertex_id -> cluster_id

    vector<unordered_map<t_id, long>> links(graph->get_vertex_count());  // Note: long is used to correctly evaluate difference

    /*
     *   Create and fill data structure
     */
    for (size_t i = 0; i < cluster_count; i++) {
        t_id_list* cluster = clusters->get_partition_vector()->at(i);

        size_t cdegree = 0;
        for (t_id vertexid: *cluster) {
            cdegree += graph->GetNeighbors(vertexid)->size();
            clustermap[vertexid] = i;
        }
        clusterdegree[i] = cdegree;
    }


    double edgeCount = 0;

    for (size_t i = 0; i < graph->get_vertex_count(); i++) {
        t_id_vector* neighbors = graph->GetNeighbors(i);
        for (size_t j = 0; j < neighbors->size(); j++) {
            t_id neighbor_id = neighbors->at(j);
            if (static_cast<t_id>(i) == neighbor_id) continue;

            t_id neighborcluster = clustermap[neighbor_id];

            size_t newvalue = 1;
            if (links[i].find(neighborcluster) != links[i].end()) {
                newvalue += links[i][neighborcluster];
            }

            links[i][neighborcluster] = newvalue;
            edgeCount++;
        }
    }
    edgeCount /= 2; // we counted all edges twice

    /*
     *   Calculate and execute vertex moves
     */
    bool improvement_found = true;
    size_t movecount = 0;
    double sum_delta_q = 0.0;
    while (improvement_found) {
        improvement_found = false;
        for (size_t vertex_id = 0; vertex_id < graph->get_vertex_count(); vertex_id++) {

            t_id best_move_cluster = -1;  // Max value, means not initialized
            double bestDeltaQ = 0;

            t_id current_cluster_id = clustermap[vertex_id];

            // for all adjacent clusters of the cluster of vertexid
            for (auto iter = links[vertex_id].begin();
                    iter != links[vertex_id].end(); ++iter) {
                t_id cluster_id = iter->first;

                if (current_cluster_id == cluster_id) continue;

                double term1 = static_cast<double>(links[vertex_id][cluster_id] -
                        links[vertex_id][current_cluster_id]) / edgeCount;
                double term2 = clusterdegree[cluster_id] -
                        clusterdegree[current_cluster_id];
                term2 += graph->GetNeighbors(vertex_id)->size();
                term2 *= graph->GetNeighbors(vertex_id)->size();
                term2 /= 2.0;
                term2 /= edgeCount;
                term2 /= edgeCount;

                double deltaQ = term1 - term2;

                if (deltaQ > bestDeltaQ) {
#ifdef DEBUG
                    assert(deltaQ <= 1.5 && "RefineCluster(), Q E (-0.5, 1] => dQ <= 1.5");
#endif // DEBUG
                    bestDeltaQ = deltaQ;
                    best_move_cluster = cluster_id;
                }
            }

            // move vertex
            if (bestDeltaQ > 0) {
#ifdef DEBUG
                assert(best_move_cluster != -1
                    && "RefineCluster(), best_move_cluster should be initialized");
#endif // DEBUG
                sum_delta_q += bestDeltaQ;
                clusterdegree[current_cluster_id] -=
                        graph->GetNeighbors(vertex_id)->size();
                clusterdegree[best_move_cluster] +=
                        graph->GetNeighbors(vertex_id)->size();

                for (size_t i = 0; i < graph->GetNeighbors(vertex_id)->size(); i++) {
                    t_id neighborid = graph->GetNeighbors(vertex_id)->at(i);

                    links[neighborid][current_cluster_id]--;
                    if (links[neighborid].find(best_move_cluster) !=
                            links[neighborid].end())
                        links[neighborid][best_move_cluster]++;
                    else
                        links[neighborid][best_move_cluster] = 1;
                }

                clustermap[vertex_id] = best_move_cluster;
                improvement_found = true;
                movecount++;
            }
        }
    }

    Partition* resultclusters = new Partition(cluster_count);
    for (size_t i = 0; i < graph->get_vertex_count(); i++) {
        t_id x = clustermap[i];
        t_id_list* cluster = resultclusters->get_partition_vector()->at(x);
        cluster->push_back(i);
    }

    return resultclusters;
}

double ModOptimizer::GetModularityFromClustering(Graph* graph,
        Partition* clusters) {
    size_t cluster_count = clusters->get_partition_vector()->size();

    t_id_vector clustermap (graph->get_vertex_count()); // maps vertex_id -> cluster_id
    for (size_t i = 0; i < cluster_count; i++) {
        t_id_list* cluster = clusters->get_partition_vector()->at(i);
        size_t csize = 0;
        for (t_id vertex_id: *cluster) {
            clustermap[vertex_id] = i;
            csize++;
        }
    }

    typedef vector<t_row_value_map*> t_sparse_matrix;

    t_sparse_matrix e;
    for (size_t i = 0; i < cluster_count; i++)
        e.push_back(new t_row_value_map());

    size_t edge_count = 0; // will be 2*|E|
    for (size_t i = 0; i < graph->get_vertex_count(); i++) {
        t_id_vector* neighbors = graph->GetNeighbors(i);
        for (size_t j = 0; j < neighbors->size(); j++) {
            if (static_cast<t_id>(i) == neighbors->at(j)) continue; // disregard loops

            t_id from = clustermap[i];
            t_id to = clustermap[neighbors->at(j)];
            if (e[from]->find(to) != e[from]->end())
                e[from]->at(to) += 1.0;
            else e[from]->insert(make_pair(to,1.0));

            edge_count++;
        }
    }

    vector<t_value> a(cluster_count);
    for (size_t i = 0; i < cluster_count; i++) {
        a[i] = 0.0;

    	for (t_row_value_map::iterator iter = e[i]->begin();
                iter != e[i]->end(); ++iter) {
            t_index column = iter->first;
            e[i]->at(column) /= static_cast<t_value>(edge_count);
            a[i] += e[i]->at(column);
        }
    }

    double Q = 0.0;
    for (size_t i = 0; i < cluster_count; i++)
	if (e[i]->find(i) != e[i]->end())
        	Q += e[i]->at(i) - a[i] * a[i];
	else
		Q += 0 - a[i] * a[i];

    for (size_t i = 0; i < e.size(); i++)
        delete e[i];

    return Q;
}
