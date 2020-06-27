#include "PlanarSpanner.h"

#include <list>

#include "SpanningGraph.h"
#include "PolygonSpanner.h"



namespace gsnunf {

PlanarSpanner::PlanarSpanner( Graph &G, double epsilon ) {
    // Graph G already contains a Delaunay Triangulation of the point set S
    SpanningGraph SG( G );
    Graph G_P = PolygonSpanner( SG );
//    GreedySpanner( G_P, epsilon, G );
}



}
