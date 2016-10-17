##############################################################################
#                                                                            #
#    CGGC Runner                                                             #
#                                                                            #    
#    This script is intended to be used to test whether existing algorithms  #  
#    for clustering/partitioning can benefit from being embedded in the      # 
#    CGGC scheme.                                                            #
#                                                                            #
#    Reference: An Ensembling Learning Strategy for Graph Clustering         #
#               Michael Ovelgonne, Andreas Geyer-Schulz, 2012                #
#                                                                            #
##############################################################################
import sys
import subprocess
from optparse import OptionParser
import igraph
import numpy
import tempfile

def read_partition(filename, size):
    f = open(filename)
    linecount = size
    partition = numpy.empty(linecount, dtype=int)
    i = 0
    for line in f:
        id = int(line.strip().split()[1])
        partition[i] = id
        i = i + 1
    assert(i == size)   
    return partition    


#########################################################
#  Calculate the maximal overlap of a list of partitions
#########################################################
def maximal_overlap(partitions):
    overlap =  numpy.empty(len(partitions[0]), dtype=int)
    distinct_ids = {}
    dist_id = 0
    for i in range(len(partitions[0])):
        concat = 0
        # concatenate the bits of all cluster ids for one vertex
        for p in range(len(partitions)):
            concat = (concat << 32) + partitions[p][i]
        new_id =  concat
        # assign new cluster id,
        # we map a new_id to a value in the range [0,#overlap_clusters]
        if not new_id in distinct_ids:  
            distinct_ids[new_id] = dist_id
            dist_id = dist_id + 1
        overlap[i] = distinct_ids[new_id]
    return overlap         
            

#########################################################
#  Create the graph induced by a partition of a graph,
#  i.e. the induced graph has a vertex for every cluster
#  in the partition
#########################################################
def create_induced_graph(graph, core_groups_partition):
    cluster_count = max(core_groups_partition)+1
    induced_graph = igraph.Graph(cluster_count)
    edge_list = set([])
    edge_weights = {}
    for edge in graph.es:
        new_from = core_groups_partition[edge.source]
        new_to = core_groups_partition[edge.target]
        
        if new_from > new_to:
            tmp = new_to
            new_to = new_from
            new_from = tmp
    
        if (new_from, new_to) in edge_list:
            edge_weights[(new_from, new_to)] = edge_weights[(new_from, new_to)] + edge["weight"]
        else:
            edge_list.add((new_from, new_to))
            edge_weights[(new_from, new_to)] = edge["weight"]
    
    edge_list = [(int(x), int(y)) for (x,y) in edge_list]
    induced_graph.add_edges(edge_list)
    induced_graph.es["weight"] = [edge_weights[(new_from, new_to)] for (new_from, new_to) in edge_list]
    
    assert(sum(graph.es["weight"]) == sum(induced_graph.es["weight"]))
    return induced_graph    

def recursive_resolve(id, mappings, level):
    if level < len(mappings)-1:
        return recursive_resolve(mappings[level][id], mappings, level+1)
    else:
        return mappings[level][id]

#########################################################
#  From a list of iterative vertex-cluster mappings 
#  (partitions), this function maps the combined 
#  mapping back to the initial partition
#########################################################
def resolves_mapping(mappings):
    result_partition = numpy.empty(len(mappings[0]), dtype=int)
    for id in range(len(mappings[0])):
        cluster_id = recursive_resolve(id, mappings, 0)
        result_partition[id] = cluster_id    
    return result_partition

def store_partition(partition, filename):
    f = open(filename,"w")
    for i in range(len(partition)):
        f.write(str(i)+" "+str(partition[i])+"\n")
    f.close()     

def run_cggc_induced_graph(graph_filename, num_iterations, ensemble_size):
    # iteratively create core groups     
    mappings = []
    tmp_graph_filename = graph_filename
    for i in range(num_iterations):
        graph = igraph.Graph.Read_Ncol(tmp_graph_filename, weights=True)
        assert(sum(graph.es["weight"])>0)
        
        tmp_partition_filename = tempfile.NamedTemporaryFile(prefix="/tmp/").name
        partitions = []
        for j in range(ensemble_size):
            ########################################################################
            # Edit the command string according to the parameters of your program  #
            ########################################################################
            command_string = [options.solver_filename, graph_filename, tmp_partition_filename]
            #print "run command: ", command_string
            subprocess.call(command_string)
            #print "returned form subprocess"
            partition = read_partition(tmp_partition_filename, graph.vcount())
            partitions.append(partition)
            
        core_groups_partition = maximal_overlap(partitions)
        mappings.append(core_groups_partition)

        induced_graph = create_induced_graph(graph, core_groups_partition)
        tmp_graph_filename = tempfile.NamedTemporaryFile(prefix="/tmp/", suffix=".ncol").name
        induced_graph.write_ncol(tmp_graph_filename, weights="weight", names=None)

    # Final Optimization
    graph = igraph.Graph.Read_Ncol(tmp_graph_filename, weights=True)
    tmp_partition_filename = tempfile.NamedTemporaryFile(prefix="/tmp/final_partition_").name
    command_string = [options.solver_filename, tmp_graph_filename, tmp_partition_filename]
    subprocess.call(command_string)
    partition = read_partition(tmp_partition_filename, graph.vcount())
    mappings.append(partition)
    
    # Resolve the final partition from the multi-level mapping
    final_partition = resolves_mapping(mappings)      
    return final_partition


def run_cggc_partition(graph_filename, num_iterations, ensemble_size):
    # iteratively create core groups     
    mappings = []
    graph = igraph.Graph.Read_Ncol(graph_filename, weights=True)
    tmp_new_partition_filename = tempfile.NamedTemporaryFile(prefix="/tmp/initial_partition_").name
    store_partition(range(graph.vcount()), tmp_new_partition_filename)
    for i in range(num_iterations):
        partitions = []
        tmp_partition_filename = tmp_new_partition_filename
        for j in range(ensemble_size):
            ########################################################################
            # Edit the command string according to the parameters of your program  #
            ########################################################################
            tmp_weak_partition_filename = tempfile.NamedTemporaryFile(prefix="/tmp/iteration"+str(j)).name
            command_string = [options.solver_filename, graph_filename, tmp_partition_filename, tmp_weak_partition_filename]
            subprocess.call(command_string)
            partition = read_partition(tmp_weak_partition_filename, graph.vcount())
            partitions.append(partition)
            
        core_groups_partition = maximal_overlap(partitions)
        mappings.append(core_groups_partition)

        tmp_new_partition_filename = tempfile.NamedTemporaryFile(prefix="/tmp/", suffix=".partition").name
        store_partition(core_groups_partition, tmp_new_partition_filename)

    # Final Optimization
    tmp_partition_filename = tmp_new_partition_filename
    tmp_new_partition_filename = tempfile.NamedTemporaryFile(prefix="/tmp/").name
    command_string = [options.solver_filename, graph_filename, tmp_partition_filename, tmp_new_partition_filename]
    subprocess.call(command_string)
    final_partition = read_partition(tmp_new_partition_filename, graph.vcount())
    
    return final_partition

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--solver", action="store", dest="solver_filename")
    parser.add_option("-g", "--graph", action="store", dest="graph_filename")
    parser.add_option("-i", "--iterations", action="store", dest="num_iterations", type=int, default=1)    
    parser.add_option("-e", "--ensemblesize", action="store", dest="ensemble_size", type=int, default=5)
    parser.add_option("-p", "--acceptspartition", action="store_true", dest="accepts_partition", default=False)
    (options, args) = parser.parse_args()

    if options.solver_filename == None or options.graph_filename == None:
        print "Filenames for solver and graph required."
        sys.exit(1)
  
    if options.accepts_partition:
        final_partition = run_cggc_partition(options.graph_filename, options.num_iterations, options.ensemble_size) 
    else:    
        final_partition = run_cggc_induced_graph(options.graph_filename, options.num_iterations, options.ensemble_size) 

    for i in range(len(final_partition)):
        print i, final_partition[i]
        
