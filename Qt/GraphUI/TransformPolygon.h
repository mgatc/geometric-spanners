#ifndef GSNUNF_TRANSFORMPOLYGON_H
#define GSNUNF_TRANSFORMPOLYGON_H

#include <unordered_map>

#include "SpanningGraph.h"


namespace gsnunf {

//typedef unordered_map<Vertex_handle, short int> VisitsAllowedTable;

template< typename T >
typename T::VertexMap TransformPolygon( const T& G ) {
    typename T::VertexMap out;
    int visits = 0;
    for( auto A : G._E ) {
        visits = A.second.size();
        assert( visits <= 3 ); // SpanningGraph should already be guaranteed to have degree <= 3

        if( A.first->info().on_outer_face )
            --visits;

        out.insert( make_pair( A.first, visits ) );
    }
    return out;
}
} // namespace gsnunf

#endif // GSNUNF_TRANSFORMPOLYGON_H

