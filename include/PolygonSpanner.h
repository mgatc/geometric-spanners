#ifndef GSNUNF_POLYGONSPANNER_H
#define GSNUNF_POLYGONSPANNER_H

#include <utility>
#include <queue>

#include "CGALComponents.h"
#include "SpanningGraph.h"
#include "DelaunayGraph.h"



namespace gsnunf {

class PolygonSpanner : public DelaunayGraph {
  public:

//    PolygonSpanner() {}
//    ~PolygonSpanner() {}

    PolygonSpanner( SpanningGraph& SG );

  protected:

    void add_children( const DelaunayGraph &G, Vertex_circulator C );
    void add_cross_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r );
    void add_forward_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r );
    void add_polygon_spanner_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r );
    void find_s_1_in_circulator( Vertex_circulator& C, const Incident_vertices& N );

  private:

    //SpanningGraph SG;
    queue<Vertex_handle> _level;
    unordered_set<Vertex_handle> _known;


}; // class PolygonSpanner

} // namespace gsnunf

#endif // GSNUNF_POLYGONSPANNER_H
