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

#include <boost/unordered_map.hpp>

class Partition;

class ActiveRowSet {
public:
    ActiveRowSet(int size);
    ActiveRowSet(Partition* clusters);
    virtual ~ActiveRowSet();

    void Remove(int &element);
    int GetRandomElement();
    int Get(int &index);
    int GetActiveRowCount();

private:
    std::vector<int> elements_;
    boost::unordered_map<int, int> element_lookup_;
    int num_elements_;
};

#endif /* ACTIVEROWSET_H_ */
