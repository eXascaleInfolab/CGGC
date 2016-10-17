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

#include <unordered_map>
#include "basetypes.h"

using std::unordered_map;


typedef unordered_map<t_index, t_value> t_row_value_map;
typedef t_row_value_map::value_type t_row_value_map_entry;

class Graph;
class Partition;

class SparseClusteringMatrix {
public:
	SparseClusteringMatrix(Graph* graph);
	SparseClusteringMatrix(Graph* graph, Partition* clusters);
	virtual ~SparseClusteringMatrix();

	SparseClusteringMatrix(const SparseClusteringMatrix&)=delete;
	SparseClusteringMatrix& operator =(const SparseClusteringMatrix&)=delete;

	void JoinCluster(int a, int b);
	t_value& Get(t_index rowIndex, t_index columnIndex);
	t_row_value_map* GetRow(t_index rowIndex);
	t_value& GetRowSum(t_index rowIndex);
	int GetRowEntries(t_index rowIndex);

private:
	t_row_value_map* rows_; // matrix E
	t_value* row_sums_;   // vector A
	size_t dimension_;	   // number of rows/columns of E

	void init(Graph* graph);
};

#endif /* SPARSECLUSTERINGMATRIX_H_ */
