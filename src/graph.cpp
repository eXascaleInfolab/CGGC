//============================================================================
// Name        : Graph.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : reading graphs from file and storing graphs in memory
//============================================================================


#include <fstream>
#include <cstring>
#include <cstdlib>  // atoi, strtoul
#include <cstdio>  // perror
#include <cerrno>  // errno
#include <system_error>
#include <stdexcept>
#include "graph.h"

using std::ifstream;
using std::endl;
using std::system_error;
using std::domain_error;


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


//! \brief Print the error and throws exception
//!
//! \param msg=nullptr const char*  - operation-related message to be shown
//! \param origin=nullptr const char*  - origin, i.e. the value that caused the error
void raiseError(const char* msg=nullptr, const char* origin=nullptr) {
    string comment = "ERROR";
    if(msg || origin) {
        comment += ", ";
        if(msg) {
            comment += msg;
            if(origin)
                comment.append(" ").append(origin);
        } else comment += origin;
    }
    //comment.append(": ") += strerror(errno);
    perror(comment.c_str());
    throw system_error(errno, std::system_category());
}

/*
 * loads undirected graph from file
 * file must be in Pajek .net or METIS .graph file format with appropriate
 * file extension, reads only undirected, unweighted files
 */
void Graph::LoadFromFile(string filename) {
#ifdef DEBUG
    std::cout << ">> LoadFromFile(), loading the graph from " << filename << endl;
#endif // DEBUG
    vertex_count_ = 0;
    edge_count_ = 0;
    constexpr char MSG_CONVERSION_FAILED[] = "Value conversion failed of the parsed text";

    ifstream infile(filename);

    if (!infile) {
        std::cout << "Could not open file." << endl;
    }
    else if (filename.rfind(".graph") != string::npos) { // read METIS .graph file
        // Metis format of the input graph (.graph):
        // % Comments are marked with '%' symbol
        // % Header:
        // <vertices_num> <endges_num> [<format_bin> [vwnum]]
        // % Body, vertices_num lines without the comments:
        // [vsize] [vweight_1 .. vweight_vwnum] vid1 [eweight1]  vid2 [eweight2] ...
        // ...
        //
        // where:
        //  vertices_num  - the number of vertices in the network (graph)
        //  endges_num  - the number of edges (not directed, A-B and B-A counted
        //      as a single edge)
        //  format_bin - up to 3 digits {0, 1}: <vsized><vweighted><eweighted>
        //      vsized  - whether the size of each vertex is specified (vsize)
        //      vweighted  - whether the vertex weights are specified (vweight_1 .. vweight_vmnum)
        //      eweighted  - whether edges weights are specified eweight<i>
        //  vm_num  - the number of weights in each vertex (not the same as the number of edges)
        //
        // vsize  - size of the vertex, integer >= 0. NOTE: doe not used normally
        // vweight  - weight the vertex, integer >= 0
        // vid  - vertex id, integer >= 1. ATTENTION: can't be equal to 0
        // eweight  - edge weight, integer >= 1

        string line;

        // Read the header skipping the comments '%' and empty lines
        char *tok;
        bool eweighted = false;
        size_t vweights = 0;  // vmnum
        bool vsized = false;
        while(getline(infile, line)) {
            // Skip comments starting with '%' and empty strings
            tok = strtok(const_cast<char*>(line.data()), " ");  // The number of nodes
            if(!tok || tok[0] == '%')
                continue;

            // Parse the header
            errno = 0;
            vertex_count_ = strtoul(tok, nullptr, 10);
            if(errno)
                raiseError(MSG_CONVERSION_FAILED, tok);

            tok = strtok(nullptr, " ");  // The number of edges (not directed)
            if(tok)
                tok = strtok(nullptr, " ");  // Weights of nodes and links in binary format: 011
            errno = 0;
            const unsigned format = tok ? strtoul(tok, nullptr, 2) : 0;
            if(errno)
                raiseError(MSG_CONVERSION_FAILED, tok);
            eweighted = format & 0b1;
            vweights = format & 0b10;
            vsized = format & 0b100;

            // Read vwnum if requried
            if(tok && vweights)
                tok = strtok(nullptr, " ");
            vweights = strtoul(tok, nullptr, 10);
            if(errno)
                raiseError(MSG_CONVERSION_FAILED, tok);
            break;
        }

        // Parse body
        t_id  from = 0;  // Currently loading vertex id
        while(getline(infile, line)) {
            // Skip comments starting with '%' and empty strings
            tok = strtok(const_cast<char*>(line.data()), " ");  // The number of nodes
            if(!tok || tok[0] == '%')
                continue;

            // Skip vertex size
            if(vsized)
                strtok(nullptr, " ");
            // Skip vertex weights
            if(vweights)
                for(size_t i = 0; tok && i < vweights; tok = strtok(nullptr, " "));

            // Parse vertices and their weights
            neighbors_.push_back(new t_id_vector());
            for(; tok != nullptr; tok = strtok(nullptr, " ")) {
                errno = 0;
                t_id to = strtoul(tok, nullptr, 10);
                if(errno)
                    raiseError(MSG_CONVERSION_FAILED, tok);
                // Id of the vertices >= 1 by specification, but internally we store from if = 0
                if(to < 1)
                    throw domain_error("Id validation failed (Metis specifies vertex id >= 1)\n");
                to -= 1;
                if (from != to) {
                    neighbors_.at(from)->push_back(to);
                    edge_count_++;
                }

                // Read and skip weight if required, because this algorithm supports
                // only unweighted input networks
                if(eweighted && !strtok(nullptr, " "))
                    break;
            }
            ++from;
        }
        edge_count_ /= 2;
    }
    else if (filename.rfind(".net") != string::npos) { // read Pajek .net file
        fputs("ERROR: Pajek format is not implemented\n", stderr);
        throw domain_error("Not implemented\n");
        // Note: Below only some specially formated Pajek files are parsed well,
        // the implementation below is error-prone for the general Pajek format.

        string line;
        getline(infile, line);
        char *tok;
        tok = strtok(const_cast<char*>(line.data()), " "); // Skip *Vertices
        tok = strtok(nullptr, " ");

        errno = 0;
        vertex_count_ = strtoul(tok, nullptr, 10); // read vertex count
        if(errno)
            raiseError(MSG_CONVERSION_FAILED, tok);

        // read vertex degrees
        for (size_t i = 0; i < vertex_count_; i++) { // read vertex degree
            getline(infile, line);
            neighbors_.push_back(new t_id_vector());
        }

        getline(infile, line); // skip "*Edges"

        // read edges
        while (getline(infile, line)) {
            char *tok;

            tok = strtok(const_cast<char*>(line.data()), " ");
            errno = 0;
            t_id from = atol(tok) - 1;
            if(errno)
                raiseError(MSG_CONVERSION_FAILED, tok);

            tok = strtok(nullptr, " ");
            errno = 0;
            t_id to = atol(tok) - 1;
            if(errno)
                raiseError(MSG_CONVERSION_FAILED, tok);

            tok = strtok(nullptr, " ");

            if (from != to) { // do not add loops
                neighbors_.at(from)->push_back(to);
                neighbors_.at(to)->push_back(from);

                edge_count_++;
            }
        }
    }
    else throw domain_error("Unsupported file format\n");
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
    for (t_id vertex_id: *vertexlist) {
        t_id_vector* t_neighbors = ingraph->GetNeighbors(vertex_id);
        t_id from = reverse_mapping->at(vertex_id);

        for (size_t j = 0; j < t_neighbors->size(); j++) {
            // if edge goes to vertex outside of sub-group of vertices (vertexlist) ignore this edge
            if (reverse_mapping->find(t_neighbors->at(j)) == reverse_mapping->end())
                continue;

            t_id to = reverse_mapping->at(t_neighbors->at(j));

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

void recursive_visit(Graph* graph, t_id_list* cluster, t_id i, t_bool_vector* visited) {
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
    for (size_t i = 0; i < vertex_count_; i++)
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
