#ifndef GSNUNF_GRAPH_H
#define GSNUNF_GRAPH_H

#include <set>
#include <unordered_map>
#include <list>
#include <memory>

#include "CGALComponents.h"



namespace gsnunf {

class DelaunayGraph {
  public:
    typedef set< Vertex_handle > Incident_vertices;
    typedef unordered_map< Vertex_handle, Incident_vertices > Adjacency_list;

    shared_ptr<DelaunayTriangulation> _DT;
    Adjacency_list _E;

    DelaunayGraph() {}
    ~DelaunayGraph() {}

    DelaunayGraph( DelaunayGraph& G );
    DelaunayGraph( DelaunayTriangulation* DT );
    DelaunayGraph( shared_ptr<DelaunayTriangulation> DT );

    void add_edge( Vertex_handle v1, Vertex_handle v2 );
    void remove_edge( Vertex_handle v1, Vertex_handle v2 );

    int count_valid_neighbors( Vertex_circulator C );
    double get_angle( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r );
    void normalize_circulator( Vertex_circulator &C );
    void normalize_properties();

  protected:

    void add_half_edge( Vertex_handle v1, Vertex_handle v2 );
    void remove_half_edge( Vertex_handle v1, Vertex_handle v2 );
    void update_incident_chords( Vertex_handle v_update, std::unordered_set<Vertex_handle> &eligible_vertices, std::unordered_set<Vertex_handle> &reserved_vertices, int vertices_remaining, bool update_neighbors = false );
    void update_eligible_vertices( std::unordered_set<Vertex_handle> &eligible_vertices, Vertex_handle v, std::unordered_set<Vertex_handle> &reserved );
    void canonical_order( std::list<Vertex_handle> &out );

  private:


}; // class DelaunayGraph

} // namespace gsnunf

#endif // GSNUNF_GRAPH_H