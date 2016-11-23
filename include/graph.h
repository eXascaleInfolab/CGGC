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
#include <unordered_map>
#include <fstream>

#include "basetypes.h"
#include "partition.h"

using std::string;
using std::unordered_map;
using std::ifstream;


typedef list<t_idpair>  t_idpair_list;
typedef unordered_map<t_id, t_id>  t_id_id_map;

class Graph {
public:
    static char inputFormat(const string& filename);

    Graph(const string& filename, char fmt);
    Graph(Graph* ingraph, t_id_list* vertexlist);
    Graph(size_t vertexcount, t_idpair_list* elist);
    ~Graph();

    Graph(Graph&&)=default;
    Graph(const Graph&)=delete;

    Graph& operator =(Graph&&);
    Graph& operator =(const Graph&)=delete;

    size_t get_vertex_count();
    size_t get_edge_count();
    t_id_id_map* get_id_mapper();

    t_id_vector* GetNeighbors(t_id vertex_id);
    Partition* GetConnectedComponents();

    unordered_map<t_id, t_id>  ieids;  // Map from internal to external ids of nodes
private:
    size_t vertex_count_;
    size_t edge_count_;  // The number of edges, not arcs!
    vector<t_id_vector*> neighbors_;
    t_id_id_map* id_mapper_;
protected:
    void LoadFromFile(const string& filename, char fmt);
    void LoadSubgraph(Graph* ingraph, t_id_list* vertexlist);
    void LoadFromEdgelist(size_t vertexcount, t_idpair_list* elist);

    void loadNSL(ifstream& inpfile, bool directed);
    void loadMetis(ifstream& inpfile);
    void loadPajek(ifstream& inpfile);
};

#endif /* GRAPH_H_ */
