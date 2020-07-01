#include "SpanningGraph.h"

#include <list>

// TO BE REMOVED

namespace gsnunf {

//SpanningGraph::SpanningGraph( DelaunayGraph& G ) : DelaunayGraph( G ) {
//    compute_spanning_graph();
//}
//
//SpanningGraph::SpanningGraph( DelaunayTriangulation* DT ) : DelaunayGraph( DT ) {
//    compute_spanning_graph();
//}
//
//SpanningGraph::SpanningGraph( shared_ptr<DelaunayTriangulation> DT ) : DelaunayGraph( DT ) {
//    compute_spanning_graph();
//}
//
//
//
//void SpanningGraph::add_first_edge( Vertex_handle v, Vertex_circulator C ) {
//    Vertex_handle v2 = C->handle();
//    add_edge( v, v2 );
//}
//
//void SpanningGraph::add_second_edge( Vertex_handle v, Vertex_circulator C ) {
//
//    while( _DT->is_infinite(++C) );
//    Vertex_handle v2 = C->handle();
//    add_edge( v, v2 );
//
//}
//
//void SpanningGraph::add_last_edge( Vertex_handle v, Vertex_circulator C ) {
//    --C;
//    Vertex_circulator done(C);
//
//    while( ( C->info().is_removed || _DT->is_infinite(C) ) && --C != done );
//
//    Vertex_handle v2 = C->handle();
//
//    add_edge( v, v2 );
//}
//
//void SpanningGraph::compute_spanning_graph() {
//
//    std::list<Vertex_handle> canonical;
//    std::list<Vertex_handle>::iterator c_iter;
//    Vertex_handle triangle[3];
//    Vertex_circulator v_n, done;
//    int i;
//
//    canonical_order( canonical );
//
//    // Add first three vertices from canonical
//    for( i=0, c_iter=canonical.begin(); i<3&&i<canonical.size(); ++c_iter, ++i ) {
//        (*c_iter)->info().is_removed = false; // add vertex
//        triangle[i] = *c_iter; // save in array for quick addition of edges
//    }
//
//    for( i=0; i<3&&i<canonical.size(); ++i ) { // Add edges of triangle
//        add_edge( triangle[i], triangle[(i+1)%3] );
//    }
//
//    // Add the rest of the vertices from canonical
//    for( i=i; i<canonical.size(); ++c_iter, ++i ) {
//
//        (*c_iter)->info().is_removed = false;
//
//        v_n = _DT->incident_vertices( *c_iter );
//        done = v_n;
//
//        normalize_circulator( v_n );
//
//        done = v_n;
//
//        int k = count_valid_neighbors( v_n );
//
////        cout << "adding " << (*c_iter)->point() << "\n";
////        cout << "  k = " << k << "\n";
//
//        if( k == 2 ) {
//            // remove edge between first two vertices
//            remove_first_edge( v_n );
//            // add edge between canonical iterator and first vertex
//            add_first_edge( *c_iter, v_n );
//            // add edge between canonical iterator and second vertex
//            add_second_edge( *c_iter, v_n );
//
//        } else if( k > 2 ) {
//            // remove edge between first two vertices
//            remove_first_edge( v_n );
//            // remove edge between last two vertices
//            remove_last_edge( v_n );
//            // add edge between canonical iterator and first vertex
//            add_first_edge( *c_iter, v_n );
//            // add edge between canonical iterator and second vertex
//            add_second_edge( *c_iter, v_n );
//            // add edge between canonical iterator and last vertex
//            add_last_edge( *c_iter, v_n );
//        }
//    }
//    normalize_properties();
//
//    // Test assumption
//    for( auto it=_E.begin(); it!=_E.end(); ++it ) // for all v_i, 1<=i<=n
//        assert( it->second.size() <= 3 );         // |v_i| <= 3
//}
//
//void SpanningGraph::remove_first_edge( Vertex_circulator C ) {
//
//    Vertex_handle v1 = C->handle(),
//                  v2 = (++C)->handle();
//
//    remove_edge( v1, v2 );
//}
//
//void SpanningGraph::remove_second_edge( Vertex_circulator C ) {
//
//    Vertex_handle v1 = (++C)->handle(),
//                  v2 = (++C)->handle();
//
//    remove_edge( v1, v2 );
//}
//
//void SpanningGraph::remove_last_edge( Vertex_circulator C ) {
//    --C;
//
//    Vertex_circulator done(C);
//
//    while( ( C->info().is_removed || _DT->is_infinite(C) ) && --C != done );
//
//    Vertex_handle v1 = C->handle(),
//                  v2 = (--C)->handle();
//
//    remove_edge( v1, v2 );
//}

}

