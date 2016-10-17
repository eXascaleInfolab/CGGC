//============================================================================
// Name        : ActiveRowSet.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : stores which rows of an SparseClusteringMatrix are still in use,
//               rows and columns that are not in use have been merged with other
//               rows/columns
//============================================================================


#include "activerowset.h"

#include "partition.h"

ActiveRowSet::ActiveRowSet(int size) {
    num_elements_ = size;
    elements_.resize(size);

    for (int i = 0; i < size; i++) {
        elements_[i] = i;
        element_lookup_[i] = i;
    }
}

ActiveRowSet::ActiveRowSet(Partition* clusters) {
    elements_.resize(clusters->get_partition_vector()->size());

    for (int i = 0; i < clusters->get_partition_vector()->size(); i++) {
        elements_[i] = *(clusters->get_partition_vector()->at(i)->begin());
        element_lookup_[*(clusters->get_partition_vector()->at(i)->begin())] = i;
    }

    num_elements_ = clusters->get_partition_vector()->size();
}

ActiveRowSet::~ActiveRowSet() {
}

int ActiveRowSet::GetRandomElement() {
    int randnumber = rand() % num_elements_;
    return elements_[randnumber];
}

int ActiveRowSet::Get(int &index) {
    return elements_[index];
}

int ActiveRowSet::GetActiveRowCount() {
    return num_elements_;
}

void ActiveRowSet::Remove(int &element) {
    // copy id from row in last active bucket to bucket of deleted row
    int bucket_id = element_lookup_[element];
    element_lookup_[elements_[num_elements_ - 1]] = bucket_id;
    elements_[bucket_id] = elements_[num_elements_ - 1];
    num_elements_--;
}
