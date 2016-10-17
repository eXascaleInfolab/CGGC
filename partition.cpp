//============================================================================
// Name        : Partition.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : stores a partition
//============================================================================

#include "partition.h"

#include <algorithm>

#include <boost/foreach.hpp>


Partition::Partition(int size) {
    partition_vector_.reserve(size);
    for (int i = 0; i < size; i++)
        partition_vector_.push_back(new list<int>());
}

Partition::~Partition() {
    for (int i = 0; i<this->partition_vector_.size(); i++)
        delete partition_vector_[i];
}

void Partition::print(std::ostream file) {
    for (int i = 0; i<partition_vector_.size(); i++) {
        list<int>* cluster = partition_vector_[i];
        BOOST_FOREACH(int vertexid, *cluster) {   
            file << vertexid << " ";
        }
        file << std::endl;
    }
}

t_partition* Partition::get_partition_vector() {
    return &partition_vector_;
}

bool IsEmpty(list<int>* list) {
    bool result = (list == NULL || list->begin() == list->end());
    if (result)
        delete list;
    return result;
}

/*
 Removes from the internal partition vector those idlist entries
 that are NULL or contain no entries (empty clusters)
 */
void Partition::RemoveEmptyEntries() {
    this->partition_vector_.erase(
            std::remove_if(this->partition_vector_.begin(),
            this->partition_vector_.end(), IsEmpty),
            this->partition_vector_.end());
}