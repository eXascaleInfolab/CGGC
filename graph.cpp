//============================================================================
// Name        : Graph.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : reading graphs from file and storing graphs in memory
//============================================================================


#include <fstream>
#include <cstring>
#include "graph.h"

using std::ifstream;
using std::endl;


Graph::Graph(string filename)
: vertex_count_(0), edge_count_(0), neighbors_(), id_mapper_(nullptr) {
    LoadFromFile(filename);
}

Graph::Graph(Graph* ingraph, t_id_list* vertexlist)
: vertex_count_(0), edge_count_(0), neighbors_(), id_mapper_(nullptr) {

    LoadSubgraph(ingraph, vertexlist);
}

Graph::Graph(size_t vertexcount, t_idpair_list* elist)
: vertex_count_(0), edge_count_(0), neighbors_(), id_mapper_(nullptr) {
    LoadFromEdgelist(vertexcount, elist);
}

size_t Graph::get_vertex_count() {
    return vertex_count_;
}

size_t Graph::get_edge_count() {
    return edge_count_;
}

t_id_vector* Graph::GetNeighbors(t_id vertex_id) {
    return neighbors_.at(vertex_id);
}

t_id_id_map* Graph::get_id_mapper() {
    return id_mapper_;
}

/*
 * loads undirected graph from file
 * file must be in Pajek .net or METIS .graph file format with appropriate
 * file extension, reads only undirected, unweighted files
 */
void Graph::LoadFromFile(string filename) {
    vertex_count_ = 0;
    edge_count_ = 0;

    ifstream infile(filename);

    if (!infile) {
        std::cout << "Could not open file." << endl;
    }
    else if (filename.rfind(".graph") != string::npos) { // read METIS .graph file
        string line;
        getline(infile, line);

        char *tok;
        tok = strtok(const_cast<char*>(line.data()), " ");
        vertex_count_ = atoi(tok);

        for (int i = 0; i < vertex_count_; i++) {
            int from = i;
            neighbors_.push_back(new t_id_vector());
            getline(infile, line);

            for(char* tok = strtok(const_cast<char*>(line.data()), " ")
            ; tok != nullptr; tok = strtok(nullptr, " ")) {
                int to = atoi(tok) - 1;
                if (from != to) {
                    neighbors_.at(from)->push_back(to);
                    edge_count_++;
                }
            }
        }
        edge_count_ /= 2;
    }
    else if (filename.rfind(".net") != string::npos) { // read Pajek .net file
        string line;
        getline(infile, line);
        char *tok;
        tok = strtok(const_cast<char*>(line.data()), " "); // Skip *Vertices
        tok = strtok(nullptr, " ");

        vertex_count_ = atoi(tok); // read vertex count

        // read vertex degrees
        for (int i = 0; i < vertex_count_; i++) { // read vertex degree
            getline(infile, line);
            neighbors_.push_back(new t_id_vector());
        }

        getline(infile, line); // skip "*Edges"

        // read edges
        while (getline(infile, line)) {
            char *tok;

            tok = strtok(const_cast<char*>(line.data()), " ");
            int from = atol(tok) - 1;

            tok = strtok(nullptr, " ");
            int to = atol(tok) - 1;

            tok = strtok(nullptr, " ");

            if (from != to) { // do not add loops
                neighbors_.at(from)->push_back(to);
                neighbors_.at(to)->push_back(from);

                edge_count_++;
            }
        }
    }
    else {
        std::cerr << "Unsupported file format. Quitting." << endl;
        exit(1);
    }
}

void Graph::LoadSubgraph(Graph* ingraph, t_id_list* vertexlist) {
    vertex_count_ = vertexlist->size();
    edge_count_ = 0;

    t_id_id_map* reverse_mapping = new t_id_id_map(); // map original_id from source graph -> new id in this graph
    id_mapper_ = new t_id_id_map(); // maps new id in this graph -> original_id from source graph

    // read vertex degrees
    size_t i = 0;
    for (t_id original_id: *vertexlist) {
        neighbors_.push_back(new t_id_vector());

        t_id_id_map::value_type* rentry = new t_id_id_map::value_type(original_id, i);
        reverse_mapping->insert(*rentry);
        delete rentry;

        t_id_id_map::value_type* nentry = new t_id_id_map::value_type(i, original_id);
        id_mapper_->insert(*nentry);
        delete nentry;
        i++;
    }

    // read edges
    for (int vertex_id: *vertexlist) {
        t_id_vector* t_neighbors = ingraph->GetNeighbors(vertex_id);
        int from = reverse_mapping->at(vertex_id);

        for (size_t j = 0; j < t_neighbors->size(); j++) {
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

void Graph::LoadFromEdgelist(size_t vertexcount, t_idpair_list* elist) {
    this->vertex_count_ = vertexcount;
    for (size_t i = 0; i < vertex_count_; i++)
        neighbors_.push_back(new t_id_vector());

    edge_count_ = 0;
    typedef t_idpair t_intpair;
    for (const t_intpair& edge: *elist) {
        this->neighbors_.at(edge.first)->push_back(edge.second);
        this->neighbors_.at(edge.second)->push_back(edge.first);
        edge_count_++;
    }
}


typedef vector<bool>  t_bool_vector;

void recursive_visit(Graph* graph, t_id_list* cluster, int i, t_bool_vector* visited) {
    if (visited->at(i))
        return;

    cluster->push_back(i);
    visited->at(i) = true;
    for (size_t n = 0; n < graph->GetNeighbors(i)->size(); n++)
        if (!visited->at(graph->GetNeighbors(i)->at(n)))
            recursive_visit(graph, cluster, graph->GetNeighbors(i)->at(n), visited);
}

Partition* Graph::GetConnectedComponents() {
    t_bool_vector* visited = new t_bool_vector(this->get_vertex_count(), false);

    size_t cc_counter = 0;
    Partition* sccs = new Partition();

    for (size_t i = 0; i < this->get_vertex_count(); i++) {
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

    if(id_mapper_)
        delete id_mapper_;
    id_mapper_ = nullptr;
}

auto Graph::operator =(Graph&& gr) -> Graph&
{
    this->~Graph();
    vertex_count_ = gr.vertex_count_;
    edge_count_ = gr.edge_count_;
    neighbors_ = move(gr.neighbors_);
    id_mapper_ = gr.id_mapper_;
    gr.id_mapper_ = nullptr;

    return *this;
}
