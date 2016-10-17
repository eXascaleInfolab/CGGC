//============================================================================
// Name        : Graph.h
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : reading graphs from file and storing graphs in memory
//============================================================================


#ifndef GRAPH_H_
#define GRAPH_H_

#include <string>
#include <vector>
#include <list>

#include <boost/unordered_map.hpp>

#include "partition.h"

using namespace std;

class Graph {
public:
    Graph(std::string filename);
    Graph(Graph* ingraph, list<int>* vertexlist);
    Graph(int vertexcount, list<pair<int, int> >* elist);
    ~Graph();

    int get_vertex_count();
    int get_edge_count();
    boost::unordered_map<int, int>* get_id_mapper();
    
    vector<int>* GetNeighbors(int &vertex_id);
    Partition* GetConnectedComponents();

private:
    int vertex_count_;
    int edge_count_;
    vector<vector<int>* > neighbors_;
    boost::unordered_map<int, int>* id_mapper_;
    
    void LoadFromFile(std::string filename);
    void LoadSubgraph(Graph* ingraph, list<int>* vertexlist);
    void LoadFromEdgelist(int vertexcount, list<pair<int, int> >* elist);
};

#endif /* GRAPH_H_ */
