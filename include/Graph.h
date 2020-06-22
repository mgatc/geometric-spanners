#ifndef GRAPH_H
#define GRAPH_H

#include "CGALComponents.h"

#include <unordered_set>
#include <unordered_map>
#include <list>



namespace gsnunf {

class Graph {
public:

    typedef std::unordered_set< Vertex_handle >                    Incident_vertices;
    typedef std::unordered_map< Vertex_handle, Incident_vertices >    Adjacency_list;

    Graph(){}
    ~Graph(){}

    Graph( const Graph &G );
    Graph( const std::list<Point> &S );

    void add_edge( Vertex_handle v1, Vertex_handle v2 );
    void remove_edge( Vertex_handle v1, Vertex_handle v2 );
    int count_valid_neighbors( Vertex_circulator C );
    void normalize_circulator( Vertex_circulator &C );

protected:
    DelaunayTriangulation _DT;
    Adjacency_list _E;

    void add_half_edge( Vertex_handle v1, Vertex_handle v2 );
    void remove_half_edge( Vertex_handle v1, Vertex_handle v2 );
    void update_incident_chords( Vertex_handle v_update, std::unordered_set<Vertex_handle> &eligible_vertices, std::unordered_set<Vertex_handle> &reserved_vertices, int vertices_remaining, bool update_neighbors = false );
    void update_eligible_vertices( std::unordered_set<Vertex_handle> &eligible_vertices, Vertex_handle v, std::unordered_set<Vertex_handle> &reserved );
    void canonical_order( std::list<Vertex_handle> &out );

private:


};

}

#endif
