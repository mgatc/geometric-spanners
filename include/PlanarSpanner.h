#ifndef GSNUNF_PLANARSPANNER_H
#define GSNUNF_PLANARSPANNER_H

#include <iostream>
#include <list>

#include "Timer.h"
#include "DelaunayGraph.h"
#include "SpanningGraph.h"
#include "TransformPolygon.h"
#include "PolygonSpanner.h"
#include "GeometricSpannerPrinter.h"

namespace gsnunf {

template< typename RandomAccessIterator, typename OutputIterator >
void PlanarSpanner( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result ) {
    //GeometricSpannerPrinter printer( .25f );
    Timer t(",");

    DelaunayGraph G( pointsBegin, pointsEnd );
    SpanningGraph( G );

    SplitVertexSet V;
    SplitVertexEdgeMap P;
    {
        Timer timer(",");
        TransformPolygon( G, V, P );
    }
    {
        Timer timer(",");
        PolygonSpanner( G, V, P );
    }
    //printer.draw( G, "PolygonSpanner" );

    // send resulting edge list to output iterator
    for( auto const& adj : G._E ) {
        Vertex_handle v_1 = adj.first;
        for( auto const& v_2 : adj.second ) {
            *result = make_pair( v_1->point(), v_2->point() );
            ++result;
        }
    }
}

} // namespace gsnunf

#endif // GSNUNF_PLANARSPANNER_H
