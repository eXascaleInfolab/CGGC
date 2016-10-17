//============================================================================
// Name        : basetypes.h
// Author      : Artem Lutov
// Version     :
// Copyright   : 2016, Artem lutov <artem@exascale.info>
// Description : Definition of the base common types
//============================================================================

#ifndef BASETYPES_H
#define BASETYPES_H

#include <vector>

using std::vector;
using std::pair;


typedef int  t_id;  //! Vertex t_id type
typedef pair<t_id, t_id>  t_idpair;
typedef vector<t_id>  t_id_vector;
typedef vector<t_idpair>  t_idpair_vector;

typedef size_t  t_index;  //! Index type
typedef double  t_value;  //! Value type

#endif // BASETYPES_H
