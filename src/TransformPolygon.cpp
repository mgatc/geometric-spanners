#include "TransformPolygon.h"

#include "SpanningGraph.h"
#include "Graph.h"



namespace gsnunf {

VisitsAllowedTable TransformPolygon( const SpanningGraph &SG ) {
    VisitsAllowedTable out;
    int visits = 0, degree = 0;
    for( auto A : SG._E ) {
        degree = A.second.size();
        assert( degree <= 3 );

        if( !A.first->info().on_outer_face )
            visits = degree;
        else if( degree == 3 )
            visits = 2;
        else
            visits = 1;

        out.insert( make_pair( A.first, visits ) );
    }
    return out;
}

} // namespace gsnunf
