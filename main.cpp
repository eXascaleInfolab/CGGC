//============================================================================
// Name        : RGMC.cpp
// Author      : Michael Ovelgönne
//               extended by Artem Lutov <artem@exascale.info>, UniFR
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : Test application for randomized greedy modularity clustering
//============================================================================

#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <vector>

#include <boost/program_options.hpp>

#include "modoptimizer.h"
#include "graph.h"

using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;

namespace po = boost::program_options;


//! \brief Save resulting clusters to the specified file
//!
//! \param outfname string  - output file name
//! \param final_clusters Partition*  - clusters
//! \param graph Graph*  - input graph
//! \param fmt char  - output clusters format: members or vertices owner (cluster) id
//! \return void
void StoreClustering(string outfname, Partition* final_clusters, Graph* graph, char fmt) {
#ifdef DEBUG
    cout << ">> Saving resulting clustering to the " << outfname << endl;
#endif // DEBUG
    ofstream out;
    out.exceptions(ofstream::badbit);
    out.open(outfname.c_str());
    if (!out) {
        cerr << "Cannot open the output file.\n";
        return;
    }
    const auto&  cls = *final_clusters->get_partition_vector();
    const auto&  ieids = graph->ieids;
    switch(fmt) {
    // List member vertices of each cluster
    case 'l':  // A special case of the cnl format
    default:
        // Output the CNL header
        out <<  "# Clusters: " << cls.size() << ",  Nodes: " <<  graph->get_vertex_count()
            << ", Fuzzy: 0\n";
    case 'c':  // without the header
        for (size_t i = 0; i < cls.size(); i++) {
            for (t_id vertex_id: *(cls.at(i)))
                out << (!ieids.empty() ? ieids.at(vertex_id) : ++vertex_id) << ' ' ;  // Note: ++ to map back to the original id for Metis and Pajek
            out << endl;
        }
        break;
    // List clusters of each vertex
    case 'v':
        t_id_vector assingments(graph->get_vertex_count(), -1);
        for (size_t i = 0; i < cls.size(); i++)
            for (t_id vertex_id: *(cls.at(i)))
                assingments[vertex_id] = i + 1;  // Note: +1 to output cluster ids starting from 1
        // Input vertices in the CNL format might form non-solid range and, hence,
        // their ids should be explicitly specified
        if(!ieids.empty())
            out << "# Vertex Cluster";
        for (size_t i = 0; i < assingments.size(); i++) {
            if(!ieids.empty())
                out << ieids.at(i) << ' ';
            out << assingments[i] << endl;
        }
        break;
    }
}

int main(int argc, char* argv[]) {
    char inpfmt;
    string filename;
    char outfmt;
    string outfname;
    int k;
    int finalk;
    int runs;
    int ensemblesize;
    bool iterative;
    bool adv;
    int alg;
    int seed;

    if(argc <= 0) {
        cerr << "Error. Input arguments are expected.\n";
        return 1;
    }

    po::options_description desc(string("Performs clustering of the unweighed undirected"
        " network (graph) using RG, CGGC_RG or CGGCi_RG algorithms.\n"
        "\nUsage: ").append(argv[0]).append("[options] inpfile\n"
        "\nSupported Arguments"));

    po::positional_options_description pdesc;
    pdesc.add("inpfile", -1);

    desc.add_options()
            ("help,h", "Display this message")
            ("inpfmt,i", po::value<char> (&inpfmt)->default_value(0),
             "input network format (inferred from the file extension if not specified explicitly):"
                "\n\te - nse format: header in comments + each line specifies a single edge: <sid> <did>,"
                "\n\ta - nsa format: header in comments + each line specifies a single arc: <sid> <did>,"
                "\n\tm - Metis format of the unweighted network,"
                "\n\tp - Pajek format of the undirected unweighted network"
            )
            ("inpfile", po::value<string>(&filename), "input network (graph) file")
            ("startk,s", po::value<int>(&k)->default_value(2), "sample size of RG")
            ("finalk,f", po::value<int>(&finalk)->default_value(2000), "sample size for final RG step")
            ("runs,r", po::value<int>(&runs)->default_value(1), "number of runs from which to pick the best result")
            ("ensemblesize,e", po::value<int>(&ensemblesize)->default_value(-1), "size of ensemble for ensemble algorithms (-1 = ln(#vertices))")
            ("algorithm,a", po::value<int>(&alg)->default_value(2), "algorithm: 1: RG, 2: CGGC_RG, 3: CGGCi_RG")  // Note: 3 is the most accurate, but also the heaviest; 2 and 3 are ensemble clustering
            ("outfmt,o", po::value<char> (&outfmt)->default_value('l'), "output clusters format:"
                "\n\tl - cnl format: header in comments + each line corresponds to the cluster and contains ids of the member vertices,"
                "\n\tc - each line corresponds to the cluster and contains ids of the member vertices,"
                "\n\tv - each line corresponds to the vertex and contains id of the owner cluster"
            )
            ("outfile,c", po::value<string> (&outfname), "file to store the detected communities if required")
            ("seed,d", po::value<int> (&seed), "seed value to initialize random number generator")
            ;
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
        .options(desc).positional(pdesc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 0;
    }

    if (!vm.count("inpfile")) {
        cerr << "Error. No filename given.\n";
        cout << desc << "\n";
        return 1;
    }

    if (!vm.count("seed")) {
        time_t t;
        time(&t);
        srand((unsigned int) t);
    }
    else srand((unsigned int) seed);

    Graph graph(filename, inpfmt);

    if (ensemblesize == -1) ensemblesize = log(graph.get_vertex_count());

    switch (alg) {
        case 1:
            adv = 0;
            iterative = 0;
            break;
        case 2:
            adv = 1;
            iterative = 0;
            break;
        case 3:
            adv = 1;
            iterative = 1;
            break;
        default:
            cerr << "Error. Invalid parameter for '--algorithm'.\n";
            return 1;
    }


    clock_t start, end;
    double time;
#ifdef DEBUG
    cerr << "> Starting ModOptimizer\n";
#endif // DEBUG
    ModOptimizer gclusterer(&graph);
    start = clock();
    if (adv) {
#ifdef DEBUG
        cout << "> Starting CGGC clustering\n";
#endif // DEBUG
        gclusterer.ClusterCGGC(ensemblesize, finalk, iterative);
    } else {
#ifdef DEBUG
        cout << "> Starting RG clustering\n";
#endif // DEBUG
        gclusterer.ClusterRG(k, runs);
    }
#ifdef DEBUG
        cout << "> The clustering is completed\n";
#endif // DEBUG
    Partition* final_clusters = gclusterer.get_clusters();

    end = clock();
    time = (double(end) - double(start)) / CLOCKS_PER_SEC;

    double Q = gclusterer.GetModularityFromClustering(&graph, final_clusters);
    cout << "Q: " << Q  << ", clusters: " << final_clusters->get_partition_vector()->size()
        << ", time [sec]: "<< time << endl;

    if (!outfname.empty())
        StoreClustering(outfname, final_clusters, &graph, outfmt);
}

