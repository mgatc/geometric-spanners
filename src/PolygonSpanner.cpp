#include "PolygonSpanner.h"

#include <list>



namespace gsnunf {

PolygonSpanner::PolygonSpanner( const SpanningGraph &SG ) : Graph( SG ) {
    // BFS of SG

    // Declare queue
    // Add any valid vertex to queue

    // Process v_1:
    // Vertex_handle v_i = queue.pop(), s_1, s_m;

    // C = _DT.incident_vertices( v_i ),
    // done(C);

    // add_children( SG, C );

    // if( std::tie( s_1, s_m ) = get_next_valid_pair(C) ) {
        // add_polygon_spanner_edges( s_1, v_i, s_m );
    // }
    //
    // Vertex_circulator C;
    //

    // while( !queue.empty() )
        // v_i = queue.pop();
        // v_i->info().is_removed = true;

        // C = _DT.incident_vertices( v_i ),
        // done(C);

        // queue all in C that are not already removed or queued

        // do {

            // if( std::tie( s_1, s_m ) = get_next_valid_pair(C) ) {

                // find s_j and s_k

                // add_polygon_spanner_edges( s_1, v_i, s_j );
                // add_polygon_spanner_edges( s_k, v_i, s_m );

            // }

        // } while( ++C != done );


} // PolygonSpanner::PolygonSpanner( SpanningGraph &SG )



void PolygonSpanner::add_children( const Graph &G, Vertex_circulator C ) {

}


void PolygonSpanner::add_cross_delaunay_edges( const Vertex_handle &a, const Vertex_handle &b, const Vertex_handle &c ) {

}

void PolygonSpanner::add_incident_delaunay_edges( const Vertex_handle &a, const Vertex_handle &b, const Vertex_handle &c ) {

}

void PolygonSpanner::add_polygon_spanner_edges( const Vertex_handle &a, const Vertex_handle &b, const Vertex_handle &c ) {
    add_incident_delaunay_edges( a, b, c );
    add_cross_delaunay_edges( a, b, c );
}

double PolygonSpanner::get_angle( const Vertex_handle &a, const Vertex_handle &b, const Vertex_handle &c ) {

}

std::optional<std::pair<Vertex_handle,Vertex_handle> > PolygonSpanner::get_next_valid_pair( Vertex_circulator &C ) {

}

} // namespace gsnunf

