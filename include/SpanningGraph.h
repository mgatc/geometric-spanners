#ifndef GSNUNF_SPANNINGGRAPH_H
#define GSNUNF_SPANNINGGRAPH_H

#include <vector>

#include "DelaunayGraph.h"

namespace gsnunf {

using namespace std;

namespace spanning_graph {

void add_first_edge( DelaunayGraph& G, Vertex_handle v, Vertex_circulator C ) {
    Vertex_handle v2 = C->handle();
    G.add_edge( v, v2 );

    //G.addToEventQueue( v, 1 );
    //G.addToEventQueue( v2, 2 );
    //G.addToEventQueue( { v, v2 }, true );
}

void add_second_edge( DelaunayGraph& G, Vertex_handle v, Vertex_circulator C ) {
    while( G._DT.is_infinite(++C) );
    Vertex_handle v2 = C->handle();
    G.add_edge( v, v2 );

//        G.addToEventQueue( v, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v, v2 }, true );
}

void add_last_edge( DelaunayGraph& G, Vertex_handle v, Vertex_circulator C, const VertexHash& is_removed ) {
    --C;
    Vertex_circulator done(C);

    while( ( contains( is_removed, C ) || G._DT.is_infinite(C) ) && --C != done );

    Vertex_handle v2 = C->handle();

    G.add_edge( v, v2 );

//        G.addToEventQueue( v, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v, v2 }, true );
}

void remove_first_edge( DelaunayGraph& G, Vertex_circulator C ) {
    Vertex_handle v1 = C->handle(),
                  v2 = (++C)->handle();
    G.remove_edge( v1, v2 );

//        G.addToEventQueue( v1, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v1, v2 }, false );
}

void remove_second_edge( DelaunayGraph& G, Vertex_circulator C ) {
    Vertex_handle v1 = (++C)->handle(),
                  v2 = (++C)->handle();
    G.remove_edge( v1, v2 );

//        G.addToEventQueue( v1, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v1, v2 }, false );
}

void remove_last_edge( DelaunayGraph& G, Vertex_circulator C, const VertexHash& is_removed ) {
    --C;

    Vertex_circulator done(C);

    while( ( contains( is_removed, C ) || G._DT.is_infinite(C) ) && --C != done );

    Vertex_handle v1 = C->handle(),
                     v2 = (--C)->handle();
    G.remove_edge( v1, v2 );

//        G.addToEventQueue( v1, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v1, v2 }, false );
}

}; // namespace spanning_graph

void SpanningGraph( DelaunayGraph& G ) {
    using namespace spanning_graph;

    vector< Vertex_handle > canonical;

    G.canonical_order( inserter( canonical, canonical.end() ) );
    Timer timer(",");

    Vertex_circulator v_n, done;
    size_t i;

    VertexHash is_removed( canonical.begin(), canonical.end() );

    // Add first three vertices from canonical
    for( i=0; i<3; ++i ) { // Add edges of triangle
        is_removed.erase( canonical.at(i) );
        G.add_edge( canonical.at(i), canonical.at((i+1)%3) );
        //G.addToEventQueue( { canonical.at(i), canonical.at((i+1)%3) }, true );
        //G.addToEventQueue( { canonical.at(i), canonical.at((i+1)%3) }, true );
    }
    // Add the rest of the vertices from canonical
    for( i=i; i<canonical.size(); ++i ) {
        is_removed.erase( canonical.at(i) );
//        G.addToEventQueue( *c_iter, 0 );        // activate c_iter
//        G.addToEventQueue( *c_iter, true );

        v_n = G._DT.incident_vertices( canonical.at(i) );
        done = v_n;

        G.normalize_circulator( v_n, is_removed );
        done = v_n;

        int k = G.count_valid_neighbors( v_n, is_removed );

        if( k == 2 ) {
            // remove edge between first two vertices
            remove_first_edge( G, v_n );
            // add edge between canonical iterator and first vertex
            add_first_edge( G, canonical.at(i), v_n );
            // add edge between canonical iterator and second vertex
            add_second_edge( G, canonical.at(i), v_n );

        } else if( k > 2 ) {
            // remove edge between first two vertices
            remove_first_edge( G, v_n );
            // remove edge between last two vertices
            remove_last_edge( G, v_n, is_removed );
            // add edge between canonical iterator and first vertex
            add_first_edge( G, canonical.at(i), v_n );
            // add edge between canonical iterator and second vertex
            add_second_edge( G, canonical.at(i), v_n );
            // add edge between canonical iterator and last vertex
            add_last_edge( G, canonical.at(i), v_n, is_removed );
        }
    }

    // Test assumption
    for( auto it=G._E.begin(); it!=G._E.end(); ++it ) { // for all v_i, 1<=i<=n
        assert( it->second.size() <= 3 );                 // |v_i| <= 3
    }
}

} // namespace gsnunf

#endif // GSNUNF_SPANNINGGRAPH_H
