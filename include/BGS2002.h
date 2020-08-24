#ifndef GSNUNF_BGS2002_H
#define GSNUNF_BGS2002_H

#include <algorithm> //sort
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
    GeometricSpannerPrinter printer;
    //Timer t(",");
    //cout<<"DelaunayGraph\n";
//    vector<Point> S( pointsBegin, pointsEnd );
//    sort( S.begin(), S.end() );
    DelaunayGraph G( pointsBegin, pointsEnd );


    printer.drawEdges( G._DT, {
        { "line width", "2pt" },
        { "color", "gray" }
    });
    //cout<<"SpanningGraph\n";
//
    SpanningGraph( G );
    printer.drawEdges( G, {
        { "line width", "13pt" },
        { "color", "cyan" }
    });
//
    //size_t split_size_estimate = G._DT.number_of_vertices();
    SplitVertexSet V;
    SplitVertexEdgeMap P;
//    {
//        //Timer timer(",");
           // cout<<"TransformPolygon\n";
//
        TransformPolygon( G, V, P );
//    }
    //cout<<sizeof(*P.begin())<<"\n";
////    {
////        //Timer timer(",");
            //cout<<"PolygonSpanner\n";

        PolygonSpanner( G, V, P );

//        cout<<"V:"<< sizeof(*V.index.begin())*V.V.size()<<" "
//            <<"E:"<< sizeof(*P.begin())*P.size()<<" \n";
//    }
//    cout<<"\n";
    printer.drawEdges( G, {
        { "line width", "5pt" },
        { "color", "cyan" }
    });
    printer.drawVertices( G._DT, 5, {
        { "fill", "blue" },
        { "draw", "white" }
    });
    printer.print( "PolygonSpanner" );
    cout<<StretchFactor(G)<<",";
    cout<<G.degree()<<",";

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
