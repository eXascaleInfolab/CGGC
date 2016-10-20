# Randomized Greedy (RG) modularity clustering algorithm and the use of RG in the CGGC/CGGCi scheme.

The Paper: "An Ensemble Learning Strategy for Graph Clustering" by
Michael Ovelgönne and Andreas Geyer-Schulz at Graph Partitioning and
Graph Clustering, Contemporary Mathematics, American Mathematical
Society, 2013
Winner of modularity maximization and modularity pareto (time + quality)
challenges, 10th DIMACS Implementation Challenge

http://www.umiacs.umd.edu/~mov/

> This is a modified version of the original rgmc clustering algorithm
(algorithms: RG, CGGC_RG, CGGCi_RG) with the fixed bugs to outpt the
resulting clusters, with additional I/O (full support of the official
Metis format of the input graph) and significantly refactored. The
core algorithm itself is untouched (and the implementation is not
optimal).

## Disclamer
Copyright 2009-2012 Michael Ovelgönne and Karlsruhe Institute of Technology.
Updated by Artem V L, University of Fribourg.

Licensed under the [LGPL v2.1](License.md).

## Reference
This is a slightly altered and cleaned-up version of the code used for the
10th DIMACS Implementation Challenge. 

Reference:
An Ensemble Learning Strategy for Graph Clustering,
Michael Ovelgönne and Andreas Geyer-Schulz,
10th DIMACS Implementation Challenge - Graph Partitioning and Graph Clustering, 2012
http://www.cc.gatech.edu/dimacs10/papers/%5B18%5D-dimacs10_ovelgoennegeyerschulz.pdf

## Contact
Michael Ovelgönne mov -a-t- umiacs.umd.edu  - for the algorithm-related questions.
Please use github Issues for the implementation-related issues.

## Prerequisites ---------------------------------------------
1. Boost program options library (v. >= 1.42)
"sudo apt-get install libboost-program-options-dev"

2. Make
Only needed, if you want to use the Makefile

## Build
The source code comes with a makefile. Run make to build the program.
Compiler/Linker errors are most likely due to incorrect 
settings of include and lib paths.  


## Run
Run rgmc with the following parameters:
  --file arg               input graph file in Metis (,graph) or Pajek (.net) format (vid >= 1), UNWEIGHTED
  --k arg (=2)             sample size of RG
  --finalk arg (=2000)     sample size for final RG step
  --runs arg (=1)          number of RG runs from which to pick the best result
  --ensemblesize arg (=-1) size of ensemble for ensemble algorithms (-1 = 
                           ln(#vertices))
  --algorithm arg (=1)     algorithm: 1: RG, 2: CGGC_RG, 3: CGGCi_RG
  --outfile arg            file to store the detected communities
  --outfmt {c, v}          output file format, default: c:
    c  - each line contains unwrapped clusters, i.e. list of the member vertices (ids)
    v  - each lie corresponds to the vertex and contains it's owner cluster (id)
  --seed arg               seed value to initialize random number generator


### Example

$ rgmc --file=test.graph --algorithm=3 --outfile=test.out
runs the CGGCi_RG algorithm on the graph test.graph and writes the results to
test.out
$ rgmc --file=email.graph
