# Randomized Greedy (RG) modularity clustering algorithm with CGGC/CGGCi (Core Groups Graph ensemble Clustering Iterative) scheme.

The Paper: "An Ensemble Learning Strategy for Graph Clustering" by
Michael Ovelgönne and Andreas Geyer-Schulz at Graph Partitioning and
Graph Clustering, Contemporary Mathematics, American Mathematical
Society, 2013
Winner of modularity maximization and modularity pareto (time + quality)
challenges, 10th DIMACS Implementation Challenge  
Original sources: http://www.umiacs.umd.edu/~mov/

> This is a modified version of the original rgmc (RG, CGGC_RG, CGGCi_RG clustering
algorithms) extended with additional I/O formats and some minor fixes (mainly in
the I/O): full support of the official Metis format of the input graph, refined
parsing of the Pajek, added NSL input format and CNL output format to be applicable
in the [PyCABeM clustering benchmark](https://github.com/eXascaleInfolab/PyCABeM).
The core algorithm itself is untouched (though, the implementation is not
optimal).

## Disclamer
Copyright 2009-2012 Michael Ovelgönne and Karlsruhe Institute of Technology.  
Modified by Artem Lutov <artem@exascale.info>, University of Fribourg.

Licensed under the [LGPL v2.1](License.md).

## Reference
This is a slightly altered and cleaned-up version of the code used for the
10th DIMACS Implementation Challenge. The algorithm performs a modularity-based
greedy (RG) and optionally ensemble (CGGC) clustering of the unweighted undirected
input network (graph).

Reference:  
An Ensemble Learning Strategy for Graph Clustering,
Michael Ovelgönne and Andreas Geyer-Schulz,
10th DIMACS Implementation Challenge - Graph Partitioning and Graph Clustering, 2012
http://www.cc.gatech.edu/dimacs10/papers/%5B18%5D-dimacs10_ovelgoennegeyerschulz.pdf

## Contact
Michael Ovelgönne <mov -a-t- umiacs.umd.edu>  - for the algorithm-related questions.
Please use github Issues for the implementation-related issues.

## Prerequisites
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
```
$ ./rgmc -h
Performs clustering of the unweighed undirected network (graph) using RG, CGGC_RG or CGGCi_RG algorithms.

Supported Arguments:
  -h [ --help ]                   Display this message
  -i [ --inpfmt ] arg (=)        input network format (inferred from the file
                                  extension if not specified explicitly):
                                  e - nse format: header in comments + each
                                  line specifies a single edge: <sid> <did>,
                                  a - nse format: header in comments + each
                                  line specifies a single arc: <sid> <did>,
                                  m - Metis format of the unweighted network,
                                  p - Pajek format of the undirected unweighted
                                  network
  -s [ --startk ] arg (=2)        sample size of RG
  -f [ --finalk ] arg (=2000)     sample size for final RG step
  -r [ --runs ] arg (=1)          number of runs from which to pick the best
                                  result
  -e [ --ensemblesize ] arg (=-1) size of ensemble for ensemble algorithms (-1
                                  = ln(#vertices))
  -a [ --algorithm ] arg (=3)     algorithm: 1: RG, 2: CGGC_RG, 3: CGGCi_RG
  -o [ --outfmt ] arg (=l)        output clusters format:
                                  l - cnl format: header in comments + each
                                  line corresponds to the cluster and contains
                                  ids of the member vertices,
                                  c - each line corresponds to the cluster and
                                  contains ids of the member vertices,
                                  v - each line corresponds to the vertex and
                                  contains id of the owner cluster
  -f [ --outfile ] arg            file to store the detected communities
  -d [ --seed ] arg               seed value to initialize random number
                                  generator
```

### Example
```
$ rgmc --algorithm=2 --outfile=test.out test.graph
runs the CGGCi_RG algorithm on the graph test.graph and writes the results to
test.out
$ rgmc email.graph
```
