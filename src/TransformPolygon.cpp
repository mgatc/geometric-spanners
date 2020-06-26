#include "TransformPolygon.h"

#include "SpanningGraph.h"
#include "Graph.h"



namespace gsnunf {

VisitsAllowedTable TransformPolygon( const SpanningGraph &SG ) {
    VisitsAllowedTable out;
    int visits = 0;
    for( auto A : SG._E ) {
        visits = A.second.size();
        assert( visits <= 3 ); // SpanningGraph should already be guaranteed to have degree <= 3

        if( A.first->info().on_outer_face )
            --visits;

        out.insert( make_pair( A.first, visits ) );
    }
    return out;
}

} // namespace gsnunf
