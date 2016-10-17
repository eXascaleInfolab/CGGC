//============================================================================
// Name        : Partition.cpp
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : stores a partition
//============================================================================

#include <algorithm>
#include "partition.h"

using std::endl;


Partition::Partition(size_t size): partition_vector_() {
    partition_vector_.reserve(size);
    for (size_t i = 0; i < size; i++)
        partition_vector_.push_back(new t_id_list());
}

Partition::~Partition() {
    for (size_t i = 0; i<this->partition_vector_.size(); i++)
        delete partition_vector_[i];
}

void Partition::print(ostream& file) {
    for (size_t i = 0; i<partition_vector_.size(); i++) {
        t_id_list* cluster = partition_vector_[i];
        for (t_id vertexid: *cluster) {
            file << vertexid << " ";
        }
        file << endl;
    }
}

t_partition* Partition::get_partition_vector() {
    return &partition_vector_;
}

bool IsEmpty(t_id_list* list) {
    bool result = (list == nullptr || list->begin() == list->end());
    if (result)
        delete list;
    return result;
}

/*
 Removes from the internal partition vector those idlist entries
 that are nullptr or contain no entries (empty clusters)
 */
void Partition::RemoveEmptyEntries() {
    this->partition_vector_.erase(
            remove_if(this->partition_vector_.begin(),
            this->partition_vector_.end(), IsEmpty),
            this->partition_vector_.end());
}
