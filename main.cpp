//============================================================================
// Name        : RGMC.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : Test application for randomized greedy modularity clustering
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <string>
#include <fstream>
#include <vector>

#include <boost/version.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include "modoptimizer.h"
#include "graph.h"

namespace po = boost::program_options;

void StoreClustering(std::string out_filename, Partition* final_clusters,
        Graph* graph) {
    std::ofstream out(out_filename.data());
    if (!out) {
        std::cerr << "Cannot open output file.\n";
        return;
    }
    std::vector<int> assingments(graph->get_vertex_count(), -1);
    for (int i = 0; i < final_clusters->get_partition_vector()->size(); i++) {
        BOOST_FOREACH(int vertex_id,
                      *(final_clusters->get_partition_vector()->at(i)) ) {
            assingments[vertex_id] = i + 1;
        }    
    }
    for (int i = 0; i < graph->get_vertex_count(); i++)
        out << assingments[i] << "\n";
    out.close();
}

int main(int argc, char* argv[]) {
    std::string filename;
    std::string out_filename;
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
            ("file", po::value<std::string > (&filename), "input graph file")
            ("k", po::value<int>(&k)->default_value(2), "sample size of RG")
            ("finalk", po::value<int>(&finalk)->default_value(2000), "sample size for final RG step")
            ("runs", po::value<int>(&runs)->default_value(1), "number of runs from which to pick the best result")
            ("ensemblesize", po::value<int>(&ensemblesize)->default_value(-1), "size of ensemble for ensemble algorithms (-1 = ln(#vertices))")
            ("algorithm", po::value<int>(&alg)->default_value(1), "algorithm: 1: RG, 2: CGGC_RG, 3: CGGCi_RG")
            ("outfile", po::value<std::string> (&out_filename), "file to store the detected communities")
            ("seed", po::value<int> (&seed), "seed value to initialize random number generator")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    if (!vm.count("file")) {
        std::cout << "No filename given. Exit." << std::endl;
        exit(0);
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
            std::cout << "Invalid parameter for '--algorithm'." << std::endl;
            exit(1);
    }


    clock_t start, end;
    double time;
    ModOptimizer gclusterer(&graph);
    start = clock();
    if (adv) 
        gclusterer.ClusterCGGC(ensemblesize, finalk, iterative);
    else
        gclusterer.ClusterRG(k, runs);
    Partition* final_clusters = gclusterer.get_clusters();

    end = clock();
    time = (double(end) - double(start)) / CLOCKS_PER_SEC;

    double Q = gclusterer.GetModularityFromClustering(&graph, final_clusters);
    std::cout << "Q: " << Q  << "  time [sec]: "<< time << std::endl;

    if (vm.count("clusterfile")) {
        StoreClustering(out_filename, final_clusters, &graph);
    }
}

