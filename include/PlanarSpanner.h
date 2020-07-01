#ifndef GSNUNF_PLANARSPANNER_H
#define GSNUNF_PLANARSPANNER_H

#include <list>

#include "CGALComponents.h"
#include "DelaunayGraph.h"
#include "SpanningGraph.h"
#include "PolygonSpanner.h"
#include "GeometricSpannerPrinter.h"



namespace gsnunf {

template< class T >
DelaunayGraph<T>& PlanarSpanner( const DelaunayGraph<T>& DT, double epsilon ) {

    GeometricSpannerPrinter printer( .25f );
    printer.draw( DT._DT, "Triangulation" );

    DelaunayGraph<T> G = SpanningGraph( DT );
    printer.draw( G, "SpanningGraph" );

//    for( auto it = _DT->finite_vertices_begin(); it!=_DT->finite_vertices_end(); ++it ) {
//        cout<< it->point() << " on_outer_face:"<< (it->info().on_outer_face ? "true" : "false" ) << endl;
//    }

    DelaunayGraph<T> G_P = PolygonSpanner( G );
    printer.draw( G_P, "PolygonSpanner" );
//    GreedySpanner( G_P, epsilon, G );

    return G_P;
}

namespace planar_spanner {

} // namespace planar_spanner

} // namespace gsnunf

#endif // GSNUNF_PLANARSPANNER_H
