#ifndef GSNUNF_BGS2002_H
#define GSNUNF_BGS2002_H

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
void BGS2002( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result ) {
    //GeometricSpannerPrinter printer;
    //Timer t(",");
    cout<<"DelaunayGraph\n";
    DelaunayGraph G( pointsBegin, pointsEnd );

//    //printer.drawEdges( G._DT, { {"color", "gray"} } );
    cout<<"SpanningGraph\n";
//
    SpanningGraph( G );
//    //printer.drawEdges( G, { {"", "thick"} } );
//
    size_t split_size_estimate = G._DT.number_of_vertices();
    SplitVertexSet V( split_size_estimate );
    SplitVertexEdgeMap P;
//    {
//        //Timer timer(",");
            cout<<"TransformPolygon\n";
//
        TransformPolygon( G, V, P );
//    }
//    {
//        //Timer timer(",");
            cout<<"PolygonSpanner\n";
//
        PolygonSpanner( G, V, P );
//    }
//    cout<<"\n";
//    printer.drawEdges( G, {{"color", "blue"}} );
//    printer.drawVertices( G._DT, {{"color","red"}} );
//    printer.print( "PolygonSpanner" );
//    cout<<StretchFactor(G)<<",";
//    cout<<G.degree()<<",";

    // send resulting edge list to output iterator
//    for( auto const& adj : G._E ) {
//        Vertex_handle v_1 = adj.first;
//        for( auto const& v_2 : adj.second ) {
//            *result = make_pair( v_1->point(), v_2->point() );
//            ++result;
//        }
//    }
}

} // namespace gsnunf

#endif // GSNUNF_BGS2002_H
