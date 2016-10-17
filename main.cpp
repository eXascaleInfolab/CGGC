//============================================================================
// Name        : RGMC.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : Test application for randomized greedy modularity clustering
//============================================================================

#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <vector>

//#include <boost/version.hpp>
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
//! \param out_filename string  - output file name
//! \param final_clusters Partition*  - clusters
//! \param graph Graph*  - input graph
//! \param clmbs bool  - output cluster members or vertices owner (cluster) id
//! \return void
void StoreClustering(string out_filename, Partition* final_clusters, Graph* graph, bool clmbs) {
#ifdef DEBUG
    cout << ">> Saving resulting clustering to the " << out_filename << endl;
#endif // DEBUG
    ofstream out(out_filename.c_str());
    if (!out) {
        cerr << "Cannot open output file.\n";
        return;
    }
    if(clmbs) {
        for (size_t i = 0; i < final_clusters->get_partition_vector()->size(); i++) {
            for (t_id vertex_id: *(final_clusters->get_partition_vector()->at(i)))
                out << vertex_id + 1 << ' ' ;
            out << endl;
        }
    } else {
        t_id_vector assingments(graph->get_vertex_count(), -1);
        for (size_t i = 0; i < final_clusters->get_partition_vector()->size(); i++)
            for (t_id vertex_id: *(final_clusters->get_partition_vector()->at(i)))
                assingments[vertex_id] = i + 1;
        for (size_t i = 0; i < graph->get_vertex_count(); i++)
            out << assingments[i] << "\n";
    }
    out.close();
}

int main(int argc, char* argv[]) {
    string filename;
    string out_filename;
    char out_fmt;
    int k;
    int finalk;
    int runs;
    int ensemblesize;
    bool iterative;
    bool adv;
    int alg;
    int seed;

    po::options_description desc("Supported Arguments");
    desc.add_options()
            ("help", "Display this message")
            ("file", po::value<string > (&filename), "input graph file")
            ("k", po::value<int>(&k)->default_value(2), "sample size of RG")
            ("finalk", po::value<int>(&finalk)->default_value(2000), "sample size for final RG step")
            ("runs", po::value<int>(&runs)->default_value(1), "number of runs from which to pick the best result")
            ("ensemblesize", po::value<int>(&ensemblesize)->default_value(-1), "size of ensemble for ensemble algorithms (-1 = ln(#vertices))")
            ("algorithm", po::value<int>(&alg)->default_value(1), "algorithm: 1: RG, 2: CGGC_RG, 3: CGGCi_RG")
            ("outfmt", po::value<char> (&out_fmt)->default_value('c'), "output clusters format:"
                "\n\tc - each line corresponds to the cluster and contains ids of the member vertices,"
                "\n\tv - each line corresponds to the vertex and contains id of the owner cluster")
            ("outfile", po::value<string> (&out_filename), "file to store the detected communities")
            ("seed", po::value<int> (&seed), "seed value to initialize random number generator")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 0;
    }

    if (!vm.count("file")) {
        cout << "No filename given. Exit." << endl;
        return 1;
    }

    if (vm.count("outfmt")) {
        if (out_fmt != 'v' && out_fmt != 'c') {
            cout << "Invalid file output format specified. Exit." << endl;
            return 1;
        }
        if (out_filename.empty()) {
            cout << "Output filename is not specified. Exit." << endl;
            return 1;
        }
    }

    if (!vm.count("seed")) {
        time_t t;
        time(&t);
        srand((unsigned int) t);
    }
    else {
        srand((unsigned int) seed);
    }

    Graph graph(filename);

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
            cout << "Invalid parameter for '--algorithm'." << endl;
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
    cout << "Q: " << Q  << " clusters: " << final_clusters->get_partition_vector()->size()
        << "  time [sec]: "<< time << endl;

    if (!out_filename.empty())
        StoreClustering(out_filename, final_clusters, &graph, out_fmt == 'c');
}

