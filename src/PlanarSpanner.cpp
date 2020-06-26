#include "PlanarSpanner.h"

#include <list>

#include "SpanningGraph.h"
#include "TransformPolygon.h"



namespace gsnunf {

PlanarSpanner::PlanarSpanner( std::list<Point> &S, double epsilon, Graph &G )
    : _S(S) {
    // Delaunay Triangulation of S

    SpanningGraph SG( S );
    VisitsAllowedTable visitsAllowed = TransformPolygon( SG );
    Graph G_P = PolygonSpanner( SG, visitsAllowed );
//    GreedySpanner( G_P, epsilon, G );
}



}
