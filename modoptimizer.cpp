//============================================================================
// Name        : RGModularityClusterer.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : implements the randomized greedy algorithm and the CGGC scheme
//============================================================================


#include "modoptimizer.h"

#include <boost/foreach.hpp>

#include "sparseclusteringmatrix.h"
#include "activerowset.h"
#include "graph.h"
#include "partition.h"

using namespace std;

ModOptimizer::ModOptimizer(Graph* graph) {
    graph_ = graph;
    clusters_ = NULL;
}

ModOptimizer::~ModOptimizer() {
    delete clusters_;
}

Partition* ModOptimizer::get_clusters() {
    return clusters_;
}

double ModOptimizer::ClusterRG(int k, int runs) {
    Partition* best_partition = NULL;
    double best_q = -1;

    for (int i = 0; i < runs; i++) {
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

double ModOptimizer::ClusterCGGC(int initclusters, int restartk,
        bool iterative) {
    Partition* clusterings[initclusters];
    Partition* unjoined_clusters[initclusters];

    ClusterRG(1, 1);
    unjoined_clusters[0] = clusters_;
    for (int i = 1; i < initclusters; i++) {
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
            for (int i = 1; i < initclusters; i++) {
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
}

vector<int>* ModOptimizer::GetMembershipFromPartition(Partition* partition,
                                                     int vertex_count) {
    vector<int>* membership = new vector<int>(vertex_count);
    for (int i = 0; i < partition->get_partition_vector()->size(); i++) {
        list<int>* cluster = partition->get_partition_vector()->at(i);
        BOOST_FOREACH(int vertex_id, *cluster) {
            membership->at(vertex_id) = i;
        }
    }
    return membership;
}

Partition* ModOptimizer::CompareClusters(Graph* graph,
        Partition* partition1, Partition* partition2) {
    
    vector<vector<int>* > maps;
    maps.push_back(GetMembershipFromPartition(partition1, 
                                              graph->get_vertex_count()));
    maps.push_back(GetMembershipFromPartition(partition2,
                                              graph->get_vertex_count()));

    Partition* result_clustering = new Partition();
    std::vector<bool> assigned(graph->get_vertex_count(), false);

    list<int>* newcluster = NULL;
    for (int i = 0; i < partition1->get_partition_vector()->size(); i++) {
        list<int>* cluster = partition1->get_partition_vector()->at(i);
        BOOST_FOREACH(int vertex1, *cluster) {
            if (!assigned[vertex1]) {
                newcluster = new list<int>();
                result_clustering->get_partition_vector()->push_back(newcluster);
                newcluster->push_back(vertex1);
                assigned[vertex1] = true;
            } else
                continue;

            BOOST_FOREACH(int vertex2, *cluster) {
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

double ModOptimizer::PerformJoins(int sample_size) {
    ActiveRowSet active_rows(graph_->get_vertex_count());
    SparseClusteringMatrix cluster_matrix(graph_);

    int dimension = graph_->get_vertex_count();
    vector<pair<int, int> > joins(dimension - 1);
    int best_step = -1;
    double best_step_q = -1;

    //**********
    // calc initial Q
    //**********
    double Q = 0;
    for (int i = 0; i < dimension; i++) {
        double a_i = cluster_matrix.GetRowSum(i);
        Q -= a_i * a_i;
    }

    //**********
    // perform joins
    //**********

    for (int step = 0; step < graph_->get_vertex_count() - 1; step++) {

        int max_sample = sample_size;
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
        int join_a = -1; //  the two clusters to join
        int join_b = -1;

        max_delta_q = -1;
        for (int sample_num = 0; sample_num < max_sample; sample_num++) {

            int row_num;
            if (max_sample == graph_->get_vertex_count() - 1 - step)
                row_num = active_rows.Get(sample_num);
            else
                row_num = active_rows.GetRandomElement();

            t_row_value_map* sample_row = cluster_matrix.GetRow(row_num);

            for (t_row_value_map::iterator entry = sample_row->begin(); entry != sample_row->end(); ++entry) {
                int column_num = entry->first;
                double value = entry->second;

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

            if (sample_num == max_sample - 1 && max_delta_q < 0 &&
                    max_sample < graph_->get_vertex_count() - 1 - step)
                
                if(max_sample < graph_->get_vertex_count()/2)
                    max_sample++;
                else {
                    max_sample = graph_->get_vertex_count() - 1 - step;
                    sample_num = 1;
                }    
        }
        
        // if there is no valid merge, stop merge process
        // (can only occur for unconnected graph)
        if (join_a == -1) break; 
                
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

    clusters_ = GetPartitionFromJoins(joins, best_step, NULL);
    return best_step_q;
}

Partition* ModOptimizer::PerformJoinsRestart(Graph* graph, Partition* clusters,
                                             int k_restart_) {
    SparseClusteringMatrix cluster_matrix(graph, clusters);
    ActiveRowSet active_rows(clusters);

    int dimension = clusters->get_partition_vector()->size();
    vector<pair<int, int> > joins(dimension - 1);

    int best_step = -1;
    double best_step_q = -1;

    double modularity = 0; // not the actual start value of Q,

    //**********
    // perform joins
    //**********
    for (int step = 0; step < clusters->get_partition_vector()->size() - 1; step++) {

        int max_sample = k_restart_;
        if (k_restart_ < (clusters->get_partition_vector()->size() - 1 - step)) {
            max_sample = k_restart_;
        } else {
            max_sample = clusters->get_partition_vector()->size() - 1 - step;
        }

        // *******
        // find join
        // *******
        double max_delta_q = -1;
        int join_a = -1; //  the two clusters to join
        int join_b = -1;

        max_delta_q = -1;
        for (int sample_num = 0; sample_num < max_sample; sample_num++) {
            int row_num;
            if (max_sample == clusters->get_partition_vector()->size() - 1 - step)
                row_num = active_rows.Get(sample_num);
            else
                row_num = active_rows.GetRandomElement();

            t_row_value_map* sample_row = cluster_matrix.GetRow(row_num);

            for (t_row_value_map::iterator entry = sample_row->begin();
                    entry != sample_row->end(); ++entry) {
                int column_num = entry->first;
                double value = entry->second;

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
            if (sample_num == max_sample - 1 && max_delta_q < 0 &&
                    max_sample < dimension - 1 - step)
                max_sample++;
        }

        // if there is no valid merge, stop merge process
        // (can only occur for unconnected graph)
        if (join_a == -1) break;
        
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
    Partition* new_clusters = GetPartitionFromJoins(joins, best_step, clusters);
    return new_clusters;
}

Partition* ModOptimizer::GetPartitionFromJoins(
        vector<pair<int, int> > joins,
        const int &bestStep,
        Partition* partial_partition) {
    
    Partition* result_partition = new Partition();
    if (partial_partition == NULL) { // create new singleton partition
        // Initialize clusters
        for (int i = 0; i < graph_->get_vertex_count(); i++) {
            list<int>* vlist = new list<int>();
            vlist->push_back(i);
            result_partition->get_partition_vector()->push_back(vlist);
        }
    } else { // rearrange input partition
        // we need to create the complete list
        for (int i = 0; i < graph_->get_vertex_count(); i++) {
            list<int>* vlist = new list<int>();
            result_partition->get_partition_vector()->push_back(vlist);
        }

        for (int i = 0; i < partial_partition->get_partition_vector()->size(); i++) {
            // the first element of the list determines where to put the list
            int pos = *(partial_partition->get_partition_vector()->at(i)->begin());
            BOOST_FOREACH(int vertex,
                    *(partial_partition->get_partition_vector()->at(i))) {
                result_partition->get_partition_vector()->at(pos)->
                        push_back(vertex);
            }
        }
    }

    //join clusters according to join list
    for (int step = 0; step <= bestStep; step++) {
        list<int>* list1 =
                result_partition->get_partition_vector()->at(joins[step].first);
        list<int>* list2 =
                result_partition->get_partition_vector()->at(joins[step].second);

        list1->splice(list1->end(), *list2);
        delete result_partition->get_partition_vector()->at(joins[step].second);
        result_partition->get_partition_vector()->at(joins[step].second) = NULL;
    }

    result_partition->RemoveEmptyEntries();
    return result_partition;
}

Partition* ModOptimizer::RefineCluster(Graph* graph, Partition* clusters) {
    typedef boost::unordered_map<int, int> t_id_id_mapping;

    clusters->RemoveEmptyEntries();

    int cluster_count = clusters->get_partition_vector()->size();
    vector<int> clusterdegree(cluster_count); // sum of degrees of all vertices of a cluster
    vector<int> clustermap(graph->get_vertex_count()); // maps vertex_id -> cluster_id

    vector<t_id_id_mapping> links(graph->get_vertex_count());
    //for (int i=0; i<)

    /*
     *   Create and fill data structure
     */
    for (int i = 0; i < cluster_count; i++) {
        list<int>* cluster = clusters->get_partition_vector()->at(i);

        int cdegree = 0;
        BOOST_FOREACH(int vertexid, *cluster) {
            cdegree += graph->GetNeighbors(vertexid)->size();
            clustermap[vertexid] = i;
        }
        clusterdegree[i] = cdegree;
    }


    double edgeCount = 0;

    for (int i = 0; i < graph->get_vertex_count(); i++) {
        vector<int>* neighbors = graph->GetNeighbors(i);
        for (int j = 0; j < neighbors->size(); j++) {
            int neighbor_id = neighbors->at(j);
            if (i == neighbor_id) continue;

            int neighborcluster = clustermap[neighbor_id];

            int newvalue = 1;
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
    int movecount = 0;
    double sum_delta_q = 0.0;
    while (improvement_found) {
        improvement_found = false;
        for (int vertex_id = 0; vertex_id < graph->get_vertex_count(); vertex_id++) {

            int best_move_cluster = -1;
            double bestDeltaQ = 0;

            int current_cluster_id = clustermap[vertex_id];

            // for all adjacent clusters of the cluster of vertexid
            for (t_id_id_mapping::iterator iter = links[vertex_id].begin();
                    iter != links[vertex_id].end(); ++iter) {
                int cluster_id = iter->first;

                if (current_cluster_id == cluster_id) continue;

                double term1 = (double) (links[vertex_id][cluster_id] -
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
                    bestDeltaQ = deltaQ;
                    best_move_cluster = cluster_id;
                }
            }

            // move vertex
            if (bestDeltaQ > 0) {
                sum_delta_q += bestDeltaQ;
                clusterdegree[current_cluster_id] -=
                        graph->GetNeighbors(vertex_id)->size();
                clusterdegree[best_move_cluster] +=
                        graph->GetNeighbors(vertex_id)->size();

                for (int i = 0; i < graph->GetNeighbors(vertex_id)->size(); i++) {
                    int neighborid = graph->GetNeighbors(vertex_id)->at(i);

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
    for (int i = 0; i < graph->get_vertex_count(); i++) {
        int x = clustermap[i];
        list<int>* cluster = resultclusters->get_partition_vector()->at(x);
        cluster->push_back(i);
    }

    return resultclusters;
}

double ModOptimizer::GetModularityFromClustering(Graph* graph,
        Partition* clusters) {
    int cluster_count = clusters->get_partition_vector()->size();

    vector<int> clustermap (graph->get_vertex_count()); // maps vertex_id -> cluster_id
    for (int i = 0; i < cluster_count; i++) {
        list<int>* cluster = clusters->get_partition_vector()->at(i);
        int csize = 0;
        BOOST_FOREACH (int vertex_id, *cluster) {    
            clustermap[vertex_id] = i;
            csize++;
        }
    }

    typedef boost::unordered_map<int,double> t_sparse_row_vector;
    typedef vector<t_sparse_row_vector*> t_sparse_matrix;

    t_sparse_matrix e;
    for (int i = 0; i < cluster_count; i++)
        e.push_back(new t_sparse_row_vector());

    int edge_count = 0; // will be 2*|E|
    for (int i = 0; i < graph->get_vertex_count(); i++) {
        vector<int>* neighbors = graph->GetNeighbors(i);
        for (int j = 0; j < neighbors->size(); j++) {
            if (i == neighbors->at(j)) continue; // disregard loops

            int from = clustermap[i];
            int to = clustermap[neighbors->at(j)];
	    if (e[from]->find(to) != e[from]->end())	            
		e[from]->at(to) += 1.0;
	    else
		e[from]->insert(make_pair(to,1.0));	

            edge_count++;
        }
    }

    vector<double> a(cluster_count);
    for (int i = 0; i < cluster_count; i++) {
        a[i] = 0.0;

    	for (t_sparse_row_vector::iterator iter = e[i]->begin();
                iter != e[i]->end(); ++iter) {    
            int column = iter->first;
            e[i]->at(column) /= (double) edge_count;
            a[i] += e[i]->at(column);
        }
    }

    double Q = 0.0;
    for (int i = 0; i < cluster_count; i++)
	if (e[i]->find(i) != e[i]->end())
        	Q += e[i]->at(i) - a[i] * a[i];
	else
		Q += 0 - a[i] * a[i];

    for (int i = 0; i < e.size(); i++)
        delete e[i];

    return Q;
}
