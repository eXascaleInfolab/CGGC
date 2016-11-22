//============================================================================
// Name        : Graph.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : reading graphs from file and storing graphs in memory
//============================================================================


#include <map>
#include <cstring>  // strerror
#include <cstdlib>  // atoi, strtoul
#include <cstdio>  // perror
#include <cerrno>  // errno
#include <system_error>
#include <stdexcept>
#include <cassert>
#include "graph.h"

using std::cerr;
using std::map;
using std::endl;
using std::system_error;
using std::domain_error;


// Accessory functions ---------------------------------------------------------
void tolower(string& text)
{
    for(auto& c: text)
        c = tolower(c);
}

void tolower(char* text)
{
	if(!text)
		return;
	while(*text)
		*text++ = tolower(*text);
}

// Graph definition ------------------------------------------------------------
Graph::Graph(const string& filename, char fmt)
: vertex_count_(0), edge_count_(0), neighbors_(), id_mapper_(nullptr), ieids() {
    LoadFromFile(filename, fmt);
}

Graph::Graph(Graph* ingraph, t_id_list* vertexlist)
: vertex_count_(0), edge_count_(0), neighbors_(), id_mapper_(nullptr), ieids() {

    LoadSubgraph(ingraph, vertexlist);
}

Graph::Graph(size_t vertexcount, t_idpair_list* elist)
: vertex_count_(0), edge_count_(0), neighbors_(), id_mapper_(nullptr), ieids() {
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

/**
 * loads undirected graph from file
 * file must be in Pajek .net or METIS .graph file format with appropriate
 * file extension, reads only undirected, unweighted files
 */
void Graph::LoadFromFile(const string& filename, char fmt) {
#ifdef DEBUG
    std::cout << ">> LoadFromFile(), loading the graph from " << filename << endl;
#endif // DEBUG
    vertex_count_ = 0;
    edge_count_ = 0;

    // Identify input format
    if(!fmt)
        fmt = inputFormat(filename);
    ifstream inpfile(filename);

    if(!inpfile) {
        cerr << "Could not open the file " << filename << endl;
		throw ifstream::failure(strerror(errno));
    }

    switch(fmt) {
    case 'e':
    case 'a':
        loadNSL(inpfile, fmt == 'a');
        break;
    case 'm':
        loadMetis(inpfile);
        break;
    case 'p':
        loadPajek(inpfile);
        break;
    default:
        throw domain_error("Unknown format of the input file\n");
    }
}

char Graph::inputFormat(const string& filename)
{
    constexpr char  fmtdefault = 'e';
    map<string, char>  fileExt = {
        // NSL formats
        {"nse", 'e'},
        {"nsa", 'a'},
        // Pajek network format
        {"pjk", 'p'},
        {"pajek", 'p'},
        {"net", 'p'},
        // Metis graph format
        {"mts", 'm'},
        {"graph", 'm'}
    };

    auto iext = filename.rfind('.');
    if(iext == string::npos)
        return fmtdefault;

    auto ires = fileExt.find(filename.c_str() + iext + 1);
    if(ires == fileExt.end())
        return fmtdefault;
    return ires->second;
}


void Graph::loadNSL(ifstream& finp, bool directed)
{
	// Intermediate data structures to fill internal data structures
	unordered_map<t_id, t_id>  eiids;  // Map from external to internal ids

	// Parse NSE/A file
	string line;

	// Parse the header
	// [Nodes: <nodes_num>[,]	<Links>: <links_num>[,] [Weighted: {0, 1}]]
	// Note: the comma is either always present as a delimiter or always absent
	while(getline(finp, line)) {
		// Skip empty lines
		if(line.empty())
			continue;
		// Consider only subsequent comments
		if(line[0] != '#')
			break;

		// 1. Replace the staring comment mark '#' with space to allow "#clusters:"
		// 2. Replace ':' with space to allow "Clusters:<clsnum>"
		for(size_t pos = 0; pos != string::npos; pos = line.find(':', pos + 1))
			line[pos] = ' ';

		// Parse nodes num
        char *tok = strtok(const_cast<char*>(line.data()), " \t");
        if(!tok)
			continue;
		tolower(tok);
		if(strcmp(tok, "nodes"))
			continue;
		// Read nodes num
		tok = strtok(nullptr, " \t");
		if(tok) {
			// Note: optional trailing ',' is allowed here
			vertex_count_ = strtoul(tok, nullptr, 10);
			// Read the number of links
			tok = strtok(nullptr, " \t");
			if(tok) {
				tolower(tok);
				if(directed ? !strcmp(tok, "edges") : !strcmp(tok, "arcs"))
                    throw domain_error(string(tok).insert(0, "Unexpected internal file format: ") += "\n");
				tok = strtok(nullptr, " \t");
				if(tok) {
					// Note: optional trailing ',' is allowed here
					edge_count_ = strtoul(tok, nullptr, 10);
					// Read Weighted flag
					tok = strtok(nullptr, " \t");
					if(tok && (tolower(tok), !strcmp(tok, "weighted"))
					&& (tok = strtok(nullptr, " \t")) && strtoul(tok, nullptr, 10) != 0)
						fputs("WARNING, the network is weighted and this algorithm does not support weights"
							", so the weights are omitted.\n", stderr);
				}
			}
		}
		// Get the following line to unify the payload processing
		getline(finp, line);
	}

	// Preallocate containers if possible
	if(vertex_count_)
		neighbors_.reserve(vertex_count_);

	// Parse the body
	// Note: the processing is started from the read line
    size_t iline = 0;  // Payload line index, links counter
    do {
		// Skip empty lines and comments
		if(line.empty() || line[0] == '#')
			continue;

        char *tok = strtok(const_cast<char*>(line.data()), " \t");
        if(!tok)
			continue;
		t_id sid = strtoul(tok, nullptr, 10);  // External source id
        tok = strtok(nullptr, " \t");
        if(!tok)
			throw domain_error(string(line).insert(0
				, "Destination link id is not specified in this line: ").c_str());
		t_id did = strtoul(tok, nullptr, 10);  // External destination id
		// Make the mappings and fill the nodes
		auto ies = eiids.find(sid);  // Index of the external src id
		if(ies == eiids.end()) {
			ies = eiids.emplace(sid, eiids.size()).first;
			ieids.emplace(ieids.size(), sid);
			neighbors_.push_back(new t_id_vector());
		}
		auto ied = eiids.find(did);  // Index of the external dst id
		if(ied == eiids.end()) {
			ied = eiids.emplace(did, eiids.size()).first;
			ieids.emplace(ieids.size(), did);
			neighbors_.push_back(new t_id_vector());
		}
		neighbors_[ies->second]->push_back(ied->second);
		//fprintf(stderr, "+ arc: %u %u  [%u %u]\n", ies->second, ied->second, sid, did);

		// Insert back arc in case of edges
        if(!directed)
			neighbors_[ied->second]->push_back(ies->second);
		++iline;
    } while(getline(finp, line));
    assert(eiids.size() == ieids.size() && neighbors_.size() == ieids.size()
		&& "Node mappings are not synchronized");

	// Initialize internal data structures using nodes
	if(!vertex_count_)
		vertex_count_ = neighbors_.size();
	if(!edge_count_)
		edge_count_ = iline * (1 + !directed);
}

void Graph::loadMetis(ifstream& inpfile)
{
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
    //      vsized  - whether the size of each vertex is specified (vsize), usually 0
    //      vweighted  - whether the vertex weights are specified (vweight_1 .. vweight_vmnum)
    //      eweighted  - whether edges weights are specified eweight<i>
    // 	  ATTENTION: when the fmt parameter is not provided, it is assumed that the
    // 		    vertex sizes, vertex WEIGHTS, and edge weights are all equal to 1 and NOT present in the file
    //  vwnum  - the number of weights in each vertex (not the same as the number of edges), >= 0
    //      Note: if vwnum > 0 then format_bin.vweighted should be 1
    //
    //  vsize  - size of the vertex, integer >= 0. NOTE: does not used normally
    //  vweight  - weight the vertex, integer >= 0
    //  vid  - vertex id, integer >= 1. ATTENTION: can't be equal to 0
    //  eweight  - edge weight, integer >= 1

    constexpr char MSG_CONVERSION_FAILED[] = "Value conversion failed of the parsed text";
    string line;

    // Read the header skipping the comments '%' and empty lines
    char *tok;
    bool eweighted = false;
    size_t vweights = 0;  // vmnum
    bool vsized = false;
    //bool frcsln = false;  // Force the self link to have node weight equal to 1 (see the specification)
    while(getline(inpfile, line)) {
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
        vsized = (format & 0b100) || vweights;
        //frcsln = !tok;

        // Read vwnum if requried
        if(tok && vweights) {
            tok = strtok(nullptr, " ");
            vweights = strtoul(tok, nullptr, 10);
            if(errno)
                raiseError(MSG_CONVERSION_FAILED, tok);
        }  // Otherwise vweights = 1 if format does not reset them
        break;
    }

    // Parse body
    neighbors_.reserve(vertex_count_);
    t_id  from = 0;  // Currently loading vertex id
    while(getline(inpfile, line)) {
        // Skip comments starting with '%' and empty strings
        tok = strtok(const_cast<char*>(line.data()), " ");  // The number of nodes
        if(!tok || tok[0] == '%')
            continue;

        // Skip vertex size
        if(vsized)
            strtok(nullptr, " ");
        // Skip vertex weights (this algorithm supports only unweighted networks)
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

            // Read and skip the weight if exist (this algorithm supports
            // only unweighted networks)
            if(eweighted && !strtok(nullptr, " "))
                break;
        }
        // NOTE: self weights can't be considered, because this algorithm supports
        // only unweighted networks without selflinks
        //// Consider self weight forcing
        //if(frcsln) {
        //    neighbors_.at(from)->push_back(from);
        //    edge_count_++;
        //}
        ++from;
    }
    edge_count_ /= 2;
    assert(vertex_count_ == neighbors_.size() && "Vertex number validation failed");
}

void Graph::loadPajek(ifstream& inpfile)
{
    fputs("ERROR: Pajek format is not implemented\n", stderr);
    throw domain_error("Not implemented\n");
    // Note: Below only some specially formated Pajek files are parsed well,
    // the implementation below is error-prone for the general Pajek format.

    constexpr char MSG_CONVERSION_FAILED[] = "Value conversion failed of the parsed text";
    string line;

    while(getline(inpfile, line)) {
        // Skip empty lines and comments
        if(line.empty() || line[0] == '%')
            continue;
        if(line[0] == '*') {
            tolower(line);
            // Read the number of vertices
            if(line == "*vertices") {
                auto vnumstr = line.c_str() + strlen("*vertices");
                errno = 0;
                vertex_count_ = strtoul(vnumstr, nullptr, 10); // read vertex count
                if(errno)
                    raiseError(MSG_CONVERSION_FAILED, vnumstr);
                break;
            }

            // Check for the edges
            if(line == "*edges")
                break;
        }

    }
    if(line.empty())
        throw domain_error("The pajek input file does not contain an Edges section\n");

    // Load the edges
    errno = 0;
    while (getline(inpfile, line)) {
        char *tok = strtok(const_cast<char*>(line.data()), " ");
        // Id of the vertices >= 1 by specification, but internally we store from if = 0
        t_id from = atol(tok) - 1;
        if(errno)
            raiseError(MSG_CONVERSION_FAILED, tok);

        tok = strtok(nullptr, " ");
        t_id to = atol(tok) - 1;
        if(errno)
            raiseError(MSG_CONVERSION_FAILED, tok);

        //tok = strtok(nullptr, " ");  // Weights are not supported the the clustering algorithm
        if (from != to) { // do not add loops
            neighbors_.at(from)->push_back(to);
            neighbors_.at(to)->push_back(from);

            edge_count_++;
        }
    }
    assert((!vertex_count_ || vertex_count_ == neighbors_.size())
        && "Vertex number validation failed");
    if(!vertex_count_)
        vertex_count_ = neighbors_.size();
}

void Graph::LoadSubgraph(Graph* ingraph, t_id_list* vertexlist) {
    vertex_count_ = vertexlist->size();
    edge_count_ = 0;

    t_id_id_map reverse_mapping; // map original_id from source graph -> new id in this graph
    id_mapper_ = new t_id_id_map(); // maps new id in this graph -> original_id from source graph

    // read vertex degrees
    size_t i = 0;
    for (t_id original_id: *vertexlist) {
        neighbors_.push_back(new t_id_vector());

        t_id_id_map::value_type* rentry = new t_id_id_map::value_type(original_id, i);
        reverse_mapping.insert(*rentry);
        delete rentry;

        t_id_id_map::value_type* nentry = new t_id_id_map::value_type(i, original_id);
        id_mapper_->insert(*nentry);
        delete nentry;
        i++;
    }

    // read edges
    for (t_id vertex_id: *vertexlist) {
        t_id_vector* t_neighbors = ingraph->GetNeighbors(vertex_id);
        t_id from = reverse_mapping.at(vertex_id);

        for (size_t j = 0; j < t_neighbors->size(); j++) {
            // if edge goes to vertex outside of sub-group of vertices (vertexlist) ignore this edge
            if (reverse_mapping.find(t_neighbors->at(j)) == reverse_mapping.end())
                continue;

            t_id to = reverse_mapping.at(t_neighbors->at(j));

            if (from != to) { // do not add loops
                neighbors_.at(from)->push_back(to);
                neighbors_.at(to)->push_back(from);
                edge_count_++;
            }
        }
    }
}

void Graph::LoadFromEdgelist(size_t vertexcount, t_idpair_list* elist) {
    vertex_count_ = vertexcount;
    for (size_t i = 0; i < vertex_count_; i++)
        neighbors_.push_back(new t_id_vector());

    edge_count_ = 0;
    typedef t_idpair t_intpair;
    for (const t_intpair& edge: *elist) {
        neighbors_.at(edge.first)->push_back(edge.second);
        neighbors_.at(edge.second)->push_back(edge.first);
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
    t_bool_vector visited(this->get_vertex_count(), false);

    size_t cc_counter = 0;
    Partition* sccs = new Partition();

    for (size_t i = 0; i < this->get_vertex_count(); i++) {
        if (visited.at(i)) {
            t_id_list* new_cluster = new t_id_list();
            recursive_visit(this, new_cluster, i, &visited);
            cc_counter++;
            sccs->get_partition_vector()->push_back(new_cluster);
        }
    }

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
