//============================================================================
// Name        : SparseClusteringMatrix.h
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : Storing the sparse matrix e, which stores the fractions of 
//               edges connecting a pair of vertices (row and column)
//               matrix implemented only for undirected graphs
//============================================================================


#ifndef SPARSECLUSTERINGMATRIX_H_
#define SPARSECLUSTERINGMATRIX_H_

#include <iostream>
#include <vector>
#include <list>

#include <boost/unordered_map.hpp>

typedef boost::unordered_map<int, double> t_row_value_map;
typedef boost::unordered_map<int, double>::value_type t_row_value_map_entry;


class Graph;
class Partition;

class SparseClusteringMatrix {
public:
	SparseClusteringMatrix(Graph* graph);
	SparseClusteringMatrix(Graph* graph, Partition* clusters);
	virtual ~SparseClusteringMatrix();

	void JoinCluster(int &a, int &b);
	double& Get(int &rowIndex, int &columnIndex);
	t_row_value_map* GetRow(int &rowIndex);
	double& GetRowSum(int &rowIndex);
	int GetRowEntries(int &rowIndex);

private:
	t_row_value_map* rows_; // matrix E
	double* row_sums_;   // vector A
	int dimension_;	   // number of rows/columns of E

	void init(Graph* graph);
};

#endif /* SPARSECLUSTERINGMATRIX_H_ */
