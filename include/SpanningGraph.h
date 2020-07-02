#ifndef GSNUNF_SPANNINGGRAPH_H
#define GSNUNF_SPANNINGGRAPH_H

#include <list>

#include "CGALComponents.h"
#include "DelaunayGraph.h"



namespace gsnunf {

    namespace spanning_graph {

    template< class T >
    void add_first_edge( T& G, typename T::Vertex_handle v, typename T::Vertex_circulator C ) {
        typename T::Vertex_handle v2 = C->handle();
        G.add_edge( v, v2 );
    }

    template< class T >
    void add_second_edge( T& G, typename T::Vertex_handle v, typename T::Vertex_circulator C ) {
        while( G._DT.is_infinite(++C) );
        typename T::Vertex_handle v2 = C->handle();
        G.add_edge( v, v2 );
    }

    template< class T >
    void add_last_edge( T& G, typename T::Vertex_handle v, typename T::Vertex_circulator C, typename T::VertexHash is_removed ) {
        --C;
        typename T::Vertex_circulator done(C);

        //while( ( C->info().is_removed || G._DT->is_infinite(C) ) && --C != done );
        while( ( contains( is_removed, C ) || G._DT.is_infinite(C) ) && --C != done );

        typename T::Vertex_handle v2 = C->handle();

        G.add_edge( v, v2 );
    }

    template< class T >
    void remove_first_edge( T& G, typename T::Vertex_circulator C ) {
        typename T::Vertex_handle v1 = C->handle(),
                                  v2 = (++C)->handle();

        G.remove_edge( v1, v2 );
    }

    template< class T >
    void remove_second_edge( T& G, typename T::Vertex_circulator C ) {
        typename T::Vertex_handle v1 = (++C)->handle(),
                                  v2 = (++C)->handle();

        G.remove_edge( v1, v2 );
    }

    template< class T >
    void remove_last_edge( T& G, typename T::Vertex_circulator C, typename T::VertexHash is_removed ) {
        --C;

        typename T::Vertex_circulator done(C);

        //while( ( C->info().is_removed || G._DT->is_infinite(C) ) && --C != done );
        while( ( contains( is_removed, C ) || G._DT.is_infinite(C) ) && --C != done );

        typename T::Vertex_handle v1 = C->handle(),
                                  v2 = (--C)->handle();

        G.remove_edge( v1, v2 );
    }

}; // namespace spanning_graph


template< class T >
void SpanningGraph( DelaunayGraph<T>& DT ) {

    using Vertex_handle = typename T::Vertex_handle;
    using Vertex_circulator = typename T::Vertex_circulator;
    //using VertexSet = typename DelaunayGraph<T>::VertexSet;
    using VertexHash = typename DelaunayGraph<T>::VertexHash;

    using namespace spanning_graph;

    //DelaunayGraph SG( DT );
    list<Vertex_handle> canonical;
    typename list<Vertex_handle>::iterator c_iter;

    Vertex_handle triangle[3];
    Vertex_circulator v_n, done;
    int i;

    DT.canonical_order( canonical );

    VertexHash is_removed( canonical.begin(), canonical.end() );

    // Add first three vertices from canonical
    for( i=0, c_iter=canonical.begin(); i<3&&i<canonical.size(); ++c_iter, ++i ) {
        //(*c_iter)->info().is_removed = false; // add vertex
        is_removed.erase( *c_iter );
        triangle[i] = *c_iter; // save in array for quick addition of edges
    }

    for( i=0; i<3&&i<canonical.size(); ++i ) { // Add edges of triangle
        DT.add_edge( triangle[i], triangle[(i+1)%3] );
    }

    // Add the rest of the vertices from canonical
    for( i=i; i<canonical.size(); ++c_iter, ++i ) {

        //(*c_iter)->info().is_removed = false;
        is_removed.erase( *c_iter );

        v_n = DT._DT.incident_vertices( *c_iter );
        done = v_n;

        DT.normalize_circulator( v_n, is_removed );

        done = v_n;

        int k = DT.count_valid_neighbors( v_n, is_removed );

//        cout << "adding " << (*c_iter)->point() << "\n";
//        cout << "  k = " << k << "\n";

        if( k == 2 ) {
            // remove edge between first two vertices
            remove_first_edge( DT, v_n );
            // add edge between canonical iterator and first vertex
            add_first_edge( DT, *c_iter, v_n );
            // add edge between canonical iterator and second vertex
            add_second_edge( DT, *c_iter, v_n );

        } else if( k > 2 ) {
            // remove edge between first two vertices
            remove_first_edge( DT, v_n );
            // remove edge between last two vertices
            remove_last_edge( DT, v_n, is_removed );
            // add edge between canonical iterator and first vertex
            add_first_edge( DT, *c_iter, v_n );
            // add edge between canonical iterator and second vertex
            add_second_edge( DT, *c_iter, v_n );
            // add edge between canonical iterator and last vertex
            add_last_edge( DT, *c_iter, v_n, is_removed );
        }
    }

    // Test assumption
    for( auto it=DT._E.begin(); it!=DT._E.end(); ++it ) { // for all v_i, 1<=i<=n
        assert( it->second.size() <= 3 );                 // |v_i| <= 3
    }
}

} // namespace gsnunf

#endif // GSNUNF_SPANNINGGRAPH_H
