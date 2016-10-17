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

ActiveRowSet::ActiveRowSet(size_t size)
: elements_(), element_lookup_(), num_elements_(size)  {
    elements_.resize(num_elements_);

    for (t_index i = 0; i < num_elements_; i++) {
        elements_[i] = i;
        element_lookup_[i] = i;
    }
}

ActiveRowSet::ActiveRowSet(Partition* clusters)
: elements_(), element_lookup_(), num_elements_(0) {
    elements_.resize(clusters->get_partition_vector()->size());

    for (size_t i = 0; i < clusters->get_partition_vector()->size(); i++) {
        elements_[i] = *(clusters->get_partition_vector()->at(i)->begin());
        element_lookup_[*(clusters->get_partition_vector()->at(i)->begin())] = i;
    }

    num_elements_ = clusters->get_partition_vector()->size();
}

ActiveRowSet::~ActiveRowSet() {
}

t_index ActiveRowSet::GetRandomElement() {
    t_index randnumber = rand() % num_elements_;
    return elements_[randnumber];
}

t_index ActiveRowSet::Get(t_index index) {
    return elements_[index];
}

size_t ActiveRowSet::GetActiveRowCount() {
    return num_elements_;
}

void ActiveRowSet::Remove(size_t element) {
    // copy id from row in last active bucket to bucket of deleted row
    size_t bucket_id = element_lookup_[element];
    element_lookup_[elements_[num_elements_ - 1]] = bucket_id;
    elements_[bucket_id] = elements_[num_elements_ - 1];
    num_elements_--;
}
