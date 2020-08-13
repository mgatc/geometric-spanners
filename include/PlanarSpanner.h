#ifndef GSNUNF_PLANARSPANNER_H
#define GSNUNF_PLANARSPANNER_H

#include <list>

#include "DelaunayGraph.h"
#include "SpanningGraph.h"
#include "TransformPolygon.h"
#include "PolygonSpanner.h"
#include "GeometricSpannerPrinter.h"



namespace gsnunf {

template< class T >
void PlanarSpanner( DelaunayGraph<T>& G ) {

    //GeometricSpannerPrinter printer( .25f );
    //printer.draw( G._DT, "Triangulation" );

    SpanningGraph( G );
    //printer.draw( G, "SpanningGraph" );

    SplitVertexSet<T> V;
    SplitVertexEdgeMap<T> P;

    TransformPolygon( G, V, P );
    //print( V, P );

    PolygonSpanner( G, V, P );
    //print_vertices<T>(V);
    //printer.draw( G, "PolygonSpanner" );

}

namespace planar_spanner {

} // namespace planar_spanner

} // namespace gsnunf

#endif // GSNUNF_PLANARSPANNER_H
