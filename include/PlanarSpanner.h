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
void PlanarSpanner( DelaunayGraph<T>& G, double epsilon ) {

    GeometricSpannerPrinter printer( .25f );
    //printer.draw( G._DT, "Triangulation" );
    SpanningGraph( G );
    //std::cout<<"SpanningGraph"<<std::endl;
    //printer.draw( G, "SpanningGraph" );

//    for( auto it = _DT->finite_vertices_begin(); it!=_DT->finite_vertices_end(); ++it ) {
//        cout<< it->point() << " on_outer_face:"<< (it->info().on_outer_face ? "true" : "false" ) << endl;
//    }

    PolygonSpanner( G );
    //std::cout<<"PolygonSpanner"<<std::endl;
    //printer.draw( G, "PolygonSpanner" );
//    GreedySpanner( G_P, epsilon, G );

}

namespace planar_spanner {

} // namespace planar_spanner

} // namespace gsnunf

#endif // GSNUNF_PLANARSPANNER_H
