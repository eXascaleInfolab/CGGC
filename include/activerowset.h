//============================================================================
// Name        : ActiveRowSet.h
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : stores which rows of an SparseClusteringMatrix are still in use,
//               rows and columns that are not in use have been merged with other
//               rows/columns
//============================================================================


#ifndef ACTIVEROWSET_H_
#define ACTIVEROWSET_H_

#include <vector>
#include <unordered_map>
#include "basetypes.h"

using std::vector;
using std::unordered_map;


typedef vector<t_index>  t_index_vector;
typedef unordered_map<t_index, t_index>  t_index_index_map;

class Partition;

class ActiveRowSet {
public:
    ActiveRowSet(size_t size);
    ActiveRowSet(Partition* clusters);
    virtual ~ActiveRowSet();

    void Remove(t_index element);
    t_index GetRandomElement();
    t_index Get(t_index index);
    size_t GetActiveRowCount();

private:
    t_index_vector elements_;
    t_index_index_map element_lookup_;
    size_t num_elements_;
};

#endif /* ACTIVEROWSET_H_ */
