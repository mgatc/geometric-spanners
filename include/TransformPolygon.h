#ifndef GSNUNF_TRANSFORMPOLYGON_H
#define GSNUNF_TRANSFORMPOLYGON_H

#include <utility>

#include "DelaunayGraph.h"
#include "SplitVertex.h"

namespace gsnunf {

using namespace std;

namespace transform_polygon {

Vertex_handle find_s_1_handle( const DelaunayGraph& SG, const pair< const Vertex_handle, VertexMap< size_t > >& unsplit, const Vertex_handle& v_n ) {
    Vertex_handle v_i = unsplit.first;
    Vertex_circulator N = SG._DT.incident_vertices(v_i); // get circulator around unsplit.first
    while( (++N)->handle() != v_n );
    while( !contains( unsplit.second, N->handle() ) ) ++N; // orient to a neighbor in unsplit.second

    return N;
}

} // namespace transform_polygon

void TransformPolygon( const DelaunayGraph& SG, SplitVertexSet& V, SplitVertexEdgeMap& P ) {
    using namespace transform_polygon;
    P.clear();
    Vertex_handle v_inf = SG._DT.infinite_vertex();

    // Convex hull circulator
    Vertex_circulator s_1_finder = SG._DT.incident_vertices( v_inf ),
                      v_1_finder = SG._DT.incident_vertices( s_1_finder );

    v_1_finder = SG.orient_circulator( v_1_finder, v_inf );
    ++v_1_finder; // rotate once CCW

    SplitVertex s_first( s_1_finder, 0 ); // v_1 will need an s_1 to point to
    SplitVertex v_1( v_1_finder, s_first ),
                   v_next( v_1 ); // used in the traversal to hold vertex before insertion

    SplitVertex s_1 = s_first;
    SplitVertex v_i = V.at( add_vertex( V, v_1_finder, s_1 ) );
    //cout<<v_i<<"\n";

    Vertex_circulator N;

    do {
        N = SG._DT.incident_vertices( v_i.v );
        while( (--N)->handle() != s_1.v ); // rotate N until it points to s_1

        const VertexSet& N_SG = SG._E.find( v_i.v )->second; // get neighbors of v_i in SG._E
        while( !contains( N_SG, --N ) ); // rotate CW until reaching a neighbor in SG
        s_1 = v_i;
        v_next = SplitVertex( N, s_1 );

        if( v_next == v_1 ) {       // if we've looped all the way around
            SplitVertex& v_ref = V.V.at(v_1.key);
            v_ref.s_1 = s_1.key; // update v_1's s_1
            v_i = v_ref; // update v_i to explicitly point to v_1
        } else { // add new vertex
            v_i = V.at( add_vertex( V, N, s_1 ) );
        }

        assert( !SG._DT.is_infinite(N) );
        add_edge( SG, P, v_i, s_1 );

    } while( v_i != v_1 );

    for( auto e=SG._DT.finite_edges_begin(); e!=SG._DT.finite_edges_end(); ++e ) {
        Vertex_handle u = e->first->vertex( (e->second+1)%3 ),
                         v = e->first->vertex( (e->second+2)%3 );
        if( !contains( SG._E.at(u), v ) ) { // only add edges that are not in SG
            pair< const Vertex_handle, VertexMap< key_type > >&
                unsplit_u = *V.index.find( u ),
                unsplit_v = *V.index.find( v );
            Vertex_handle s_1_u = find_s_1_handle( SG, unsplit_u, v ),
                             s_1_v = find_s_1_handle( SG, unsplit_v, u );
            add_edge( SG, P, V.at( u, s_1_u ), V.at( v, s_1_v ) );
        }
    }
}


} // namespace gsnunf

#endif // GSNUNF_TRANSFORMPOLYGON_H

