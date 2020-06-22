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

    PolygonSpanner( const SpanningGraph &SG );

  protected:

    void add_children( const Graph &G, Vertex_circulator C );
    void add_cross_delaunay_edges( const Vertex_handle &a, const Vertex_handle &b, const Vertex_handle &c );
    void add_incident_delaunay_edges( const Vertex_handle &a, const Vertex_handle &b, const Vertex_handle &c );
    void add_polygon_spanner_edges( const Vertex_handle &a, const Vertex_handle &b, const Vertex_handle &c );
    double get_angle( const Vertex_handle &a, const Vertex_handle &b, const Vertex_handle &c );
    std::optional<std::pair<Vertex_handle,Vertex_handle> >get_next_valid_pair( Vertex_circulator &C );

  private:


}; // class PolygonSpanner

} // namespace gsnunf

#endif // GSNUNF_POLYGONSPANNER_H
