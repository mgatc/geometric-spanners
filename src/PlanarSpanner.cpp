#include "PlanarSpanner.h"

#include <list>

#include "SpanningGraph.h"
#include "PolygonSpanner.h"
#include "GeometricSpannerPrinter.h"

// TO BE REMOVED

namespace gsnunf {

//PlanarSpanner::PlanarSpanner( shared_ptr<DelaunayTriangulation> DT, double epsilon ) : DelaunayGraph(DT) {
//
//    GeometricSpannerPrinter printer( .25f );
//    printer.draw( *_DT, "Triangulation" );
//
//    SpanningGraph SG( _DT );
//    printer.draw( SG, "SpanningGraph" );
//
////    for( auto it = _DT->finite_vertices_begin(); it!=_DT->finite_vertices_end(); ++it ) {
////        cout<< it->point() << " on_outer_face:"<< (it->info().on_outer_face ? "true" : "false" ) << endl;
////    }
//
//    PolygonSpanner G_P( SG );
//    printer.draw( G_P, "PolygonSpanner" );
////    GreedySpanner( G_P, epsilon, G );
//}



}
