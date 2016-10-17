//============================================================================
// Name        : Partition.h
// Author      : Michael Ovelgönne
// Version     :
// Copyright   : 2009-2012 Karlsruhe Institute of Technology
// Description : stores a partition
//============================================================================


#ifndef PARTITION_H
#define	PARTITION_H

#include <list>
#include <iostream>
#include "basetypes.h"

using std::list;
using std::ostream;

typedef list<t_id>  t_id_list;
typedef vector<t_id_list*>  t_partition;


class Partition {
public:
    Partition(size_t size = 0);
    virtual ~Partition();

    void RemoveEmptyEntries();
    void print(ostream& file = std::cout);

    t_partition* get_partition_vector();

private:
    t_partition partition_vector_;
};

#endif	/* PARTITION_H */

