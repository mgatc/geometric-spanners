#ifndef GSNUNF_POLYGONSPANNER_H
#define GSNUNF_POLYGONSPANNER_H

#include "CGALComponents.h"
#include "SpanningGraph.h"
#include "Graph.h"

#include <utility>



namespace gsnunf {

class PolygonSpanner : public Graph {
  public:

    PolygonSpanner() {}
    ~PolygonSpanner() {}

    PolygonSpanner( const SpanningGraph& SG );

  protected:

    void add_children( const Graph &G, Vertex_circulator C );
    void add_cross_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r );
    void add_forward_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r );
    void add_polygon_spanner_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r );
    double get_angle( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r );

  private:


}; // class PolygonSpanner

} // namespace gsnunf

#endif // GSNUNF_POLYGONSPANNER_H
