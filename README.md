# Randomized Greedy (RG) modularity clustering algorithm with Core Groups Graph ensemble Clustering Iterative (CGGC/CGGCi) scheme

The Paper: [An Ensemble Learning Strategy for Graph Clustering](http://www.cc.gatech.edu/dimacs10/papers/%5B18%5D-dimacs10_ovelgoennegeyerschulz.pdf)
by Michael Ovelgönne and Andreas Geyer-Schulz at Graph Partitioning and
Graph Clustering, Contemporary Mathematics, American Mathematical
Society, 2013  
Winner of modularity maximization and modularity pareto (time + quality)
challenges, 10th DIMACS Implementation Challenge  
Original sources: http://www.umiacs.umd.edu/~mov/

Copyright 2009-2012 Michael Ovelgönne <mov -a-t- umiacs.umd.edu>, Karlsruhe Institute of Technology.  
Modified by Artem Lutov <artem@exascale.info>, University of Fribourg.

Licensed under the [LGPL v2.1](License.md).

> This is a modified version of the original rgmc (RG, CGGC_RG, CGGCi_RG clustering
algorithms) extended with additional I/O formats and some minor fixes (mainly in
the I/O): full support of the official Metis format of the input graph, refined
parsing of the Pajek, added NSL input format and CNL output format to be applicable
in the [PyCABeM clustering benchmark](https://github.com/eXascaleInfolab/PyCABeM).
The core algorithm itself is untouched (though, the implementation is not
optimal).

## Content
- [Deployment](#deployment)
	- [Dependencies](#dependencies)
	- [Compilation](#compilation)
- [Usage](#usage)
  - [Input](#input)
  - [Output](#output)
- [Related Projects](#related-projects)

# Deployment
## Dependencies
Boost program options library (v. >= 1.42)
"sudo apt-get install libboost-program-options-dev"

## Compilation
`$ make [release|debug]`  
Compiler/Linker errors are most likely due to incorrect settings of include and lib paths.  

# Usage
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

For example:
```
$ rgmc --algorithm=2 --outfile=test.out test.graph
runs the CGGCi_RG algorithm on the graph test.graph and writes the results to
test.out
```
or just `$ rgmc email.graph`.

## Input
The undirected unweighted input network to be clustered can be specified in the NSL (nsa/nse), [Metis graph](http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html) or [Pajek network](https://gephi.org/users/supported-graph-formats/pajek-net-format/) (using only `*Edges`) formats.

NSL format (nsa - arcs, directed network; nse - edges, undirected network) specifies network links in each line of the file as node ids separated by the space delimiter with optional `#` line comments and an optional header:
```
# Example Network .nse (edges - undirected)
# Nodes: 3  Edges: 3   Weighted: 0
# Note that the number of links corresponds to the number of payload lines in the file
0 1
# Empty lines and comments are allowed
0 2
2 1
```

## Output
The CNL (clusters nodes list) output is the default and standard output format, generalization of the Stanford SNAP ground-truth communities. It is an input format for various NMI-evaluation algorithms. Each line of the file corresponds to the single resulting cluster, where member nodes are specified separated by the space/tab with optional share. For example:
```
# Clusters: 2, Nodes: 5, Fuzzy: 0
0
1 3 2 4
```
The clusters output for each vertex is also supported.

# Related Projects
- [GenConvNMI](https://github.com/eXascaleInfolab/GenConvNMI) - Overlapping NMI evaluation that is compatible with the original NMI (unlike the `onmi`).
- [OvpNMI](https://github.com/eXascaleInfolab/OvpNMI) - Another method of the NMI evaluation for the overlapping clusters (communities) that is not compatible with the standard NMI value unlike GenConvNMI, but it is much faster and yields exact results unlike probabilistic results with some variance in GenConvNMI.
- [ExecTime](https://bitbucket.org/lumais/exectime/)  - A lightweight resource consumption (RSS RAM, CPU, etc.) profiler.
- [PyCABeM](https://github.com/eXascaleInfolab/PyCABeM) - Python Benchmarking Framework for the Clustering Algorithms Evaluation. Uses extrinsic (NMIs) and intrinsic (Q) measures for the clusters quality evaluation considering overlaps (nodes membership by multiple clusters).
