#ifndef GSNUNF_PLANARSPANNER_H
#define GSNUNF_PLANARSPANNER_H

#include "CGALComponents.h"
#include "Graph.h"
#include "SpanningGraph.h"



namespace gsnunf {

class PlanarSpanner {
  public:
    PlanarSpanner( std::list<Point> &S, double epsilon, Graph &G );

  protected:
    std::list<Point> &_S;

    void TransformPolygon( Graph &SG, Graph &P );
    void PolygonSpanner( Graph &P, DelaunayTriangulation &DT, Graph &G_P );
    void GreedySpanner( Graph &G_P, double epsilon, Graph &G );

}; // class PlanarSpanner

} // namespace gsnunf

#endif // GSNUNF_PLANARSPANNER_H
