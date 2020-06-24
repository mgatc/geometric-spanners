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
//    Graph G_P = PolygonSpanner( P, DT );
//    GreedySpanner( G_P, epsilon, G );
}


void PlanarSpanner::PolygonSpanner( Graph &P, DelaunayTriangulation &DT, Graph &G_P ) {

}

void PlanarSpanner::GreedySpanner( Graph &G_P, double epsilon, Graph &G ) {

}

}
