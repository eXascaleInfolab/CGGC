//============================================================================
// Name        : Graph.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : reading graphs from file and storing graphs in memory
//============================================================================


#include "graph.h"

#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#define MAX_LINE_LENGTH 100

using namespace std;
using namespace boost;



Graph::Graph(std::string filename) {
    id_mapper_ = NULL;
    LoadFromFile(filename);
}

Graph::Graph(Graph* ingraph, t_id_list* vertexlist) {
    id_mapper_ = NULL;
    LoadSubgraph(ingraph, vertexlist);
}

Graph::Graph(int vertexcount, list<pair<int, int> >* elist) {
    id_mapper_ = NULL;
    LoadFromEdgelist(vertexcount, elist);
}

int Graph::get_vertex_count() {
    return vertex_count_;
}

int Graph::get_edge_count() {
    return edge_count_;
}

vector<int>* Graph::GetNeighbors(int &vertex_id) {
    return neighbors_.at(vertex_id);
}

unordered_map<int, int>* Graph::get_id_mapper() {
    return id_mapper_;
}

/*
 * loads undirected graph from file
 * file must be in Pajek .net or METIS .graph file format with appropriate
 * file extension, reads only undirected, unweighted files
 */
void Graph::LoadFromFile(std::string filename) {
    vertex_count_ = 0;
    edge_count_ = 0;

    char line[255];
    std::ifstream infile(filename.data());

    if (!infile) {
        std::cout << "Could not open file." << std::endl;
    }
    else if (filename.rfind(".graph") != std::string::npos) { // read METIS .graph file
        char line[MAX_LINE_LENGTH];
        infile.getline(line, MAX_LINE_LENGTH);

        char *tok;
        tok = strtok(line, " ");
        vertex_count_ = atoi(tok);

        for (int i = 0; i < vertex_count_; i++) {
            std::string line;
            int from = i;
            neighbors_.push_back(new vector<int>());
            getline(infile, line);

            char_separator<char> sep(" ");
            tokenizer<char_separator<char> > tokens(line, sep);

            BOOST_FOREACH(std::string tok, tokens) {
                int to = atoi(tok.data()) - 1;
                if (from != to) {
                    neighbors_.at(from)->push_back(to);
                    edge_count_++;
                }
            }
        }
        edge_count_ = edge_count_ / 2;
    }
    else if (filename.rfind(".net") != std::string::npos) { // read Pajek .net file
        infile.getline(line, 255);
        char *tok;
        tok = strtok(line, " "); // Skip *Vertices
        tok = strtok(NULL, " ");

        vertex_count_ = atoi(tok); // read vertex count

        // read vertex degrees
        for (int i = 0; i < vertex_count_; i++) { // read vertex degree
            infile.getline(line, 255);
            neighbors_.push_back(new vector<int>());
        }

        infile.getline(line, 255); // skip "*Edges"

        // read edges
        while (infile.getline(line, 255)) {
            char *tok;

            tok = strtok(line, " ");
            int from = atol(tok) - 1;

            tok = strtok(NULL, " ");
            int to = atol(tok) - 1;

            tok = strtok(NULL, " ");

            if (from != to) { // do not add loops
                neighbors_.at(from)->push_back(to);
                neighbors_.at(to)->push_back(from);

                edge_count_++;
            }
        }
    }
    else {
        std::cerr << "Unsupported file format. Quitting." << std::endl;
        exit(1);
    }
}

void Graph::LoadSubgraph(Graph* ingraph, t_id_list* vertexlist) {
    vertex_count_ = vertexlist->size();
    edge_count_ = 0;

    typedef unordered_map<int, int> t_id_id_map;
    t_id_id_map* reverse_mapping = new t_id_id_map(); // map original_id from source graph -> new id in this graph
    id_mapper_ = new t_id_id_map(); // maps new id in this graph -> original_id from source graph

    // read vertex degrees
    int i = 0;
    BOOST_FOREACH(int original_id, *vertexlist) {
        neighbors_.push_back(new vector<int>());

        t_id_id_map::value_type* rentry = new t_id_id_map::value_type(original_id, i);
        reverse_mapping->insert(*rentry);
        delete rentry;

        t_id_id_map::value_type* nentry = new t_id_id_map::value_type(i, original_id);
        id_mapper_->insert(*nentry);
        delete nentry;
        i++;
    }

    // read edges
    BOOST_FOREACH(int vertex_id, *vertexlist) {
        vector<int>* t_neighbors = ingraph->GetNeighbors(vertex_id);
        int from = reverse_mapping->at(vertex_id);

        for (int j = 0; j < t_neighbors->size(); j++) {
            // if edge goes to vertex outside of sub-group of vertices (vertexlist) ignore this edge
            if (reverse_mapping->find(t_neighbors->at(j)) == reverse_mapping->end())
                continue;

            int to = reverse_mapping->at(t_neighbors->at(j));

            if (from != to) { // do not add loops
                neighbors_.at(from)->push_back(to);
                neighbors_.at(to)->push_back(from);
                edge_count_++;
            }
        }
    }
    
    delete reverse_mapping;
}

void Graph::LoadFromEdgelist(int vertexcount, list<pair<int, int> >* elist) {
    this->vertex_count_ = vertexcount;
    for (int i = 0; i < vertex_count_; i++)
        neighbors_.push_back(new vector<int>());

    edge_count_ = 0;
    typedef pair<int,int> t_intpair;
    BOOST_FOREACH(t_intpair edge, *elist) {
        this->neighbors_.at(edge.first)->push_back(edge.second);
        this->neighbors_.at(edge.second)->push_back(edge.first);
        edge_count_++;
    }
}

void recursive_visit(Graph* graph, t_id_list* cluster, int i, std::vector<bool>* visited) {
    if (visited->at(i))
        return;

    cluster->push_back(i);
    visited->at(i) = true;
    for (int n = 0; n < graph->GetNeighbors(i)->size(); n++)
        if (!visited->at(graph->GetNeighbors(i)->at(n)))
            recursive_visit(graph, cluster, graph->GetNeighbors(i)->at(n), visited);
}

Partition* Graph::GetConnectedComponents() {
    std::vector<bool>* visited = new std::vector<bool>(this->get_vertex_count(), false);

    int cc_counter = 0;
    Partition* sccs = new Partition();

    for (int i = 0; i < this->get_vertex_count(); i++) {
        if (!visited->at(i)) {
            t_id_list* new_cluster = new t_id_list();
            recursive_visit(this, new_cluster, i, visited);
            cc_counter++;
            sccs->get_partition_vector()->push_back(new_cluster);
        }
    }

    delete visited;
    return sccs;
}

Graph::~Graph() {
    for (int i = 0; i < vertex_count_; i++)
        delete neighbors_[i];

    delete id_mapper_;
}

