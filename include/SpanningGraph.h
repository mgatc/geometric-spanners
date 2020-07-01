#ifndef GSNUNF_SPANNINGGRAPH_H
#define GSNUNF_SPANNINGGRAPH_H

#include <list>

#include "CGALComponents.h"
#include "DelaunayGraph.h"



namespace gsnunf {

template< class T >
DelaunayGraph<T>& SpanningGraph( const DelaunayGraph<T>& DT ) {

    using Vertex_handle = typename T::Vertex_handle;
    using Vertex_circulator = typename T::Vertex_circulator;
    using VertexSet = typename T::VertexSet;


    DelaunayGraph SG( DT );
    vector<Vertex_handle> canonical;
    typename vector<Vertex_handle>::iterator c_iter;

    VertexSet is_removed;

    Vertex_handle triangle[3];
    Vertex_circulator v_n, done;
    int i;

    canonical_order( canonical );

    // Add first three vertices from canonical
    for( i=0, c_iter=canonical.begin(); i<3&&i<canonical.size(); ++c_iter, ++i ) {
        //(*c_iter)->info().is_removed = false; // add vertex
        triangle[i] = *c_iter; // save in array for quick addition of edges
    }

    for( i=0; i<3&&i<canonical.size(); ++i ) { // Add edges of triangle
        add_edge( triangle[i], triangle[(i+1)%3] );
    }

    // Add the rest of the vertices from canonical
    for( i=i; i<canonical.size(); ++c_iter, ++i ) {

        //(*c_iter)->info().is_removed = false;

        v_n = DT._DT.incident_vertices( *c_iter );
        done = v_n;

        normalize_circulator( v_n );

        done = v_n;

        int k = count_valid_neighbors( v_n );

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
//            normalize_properties();

    // Test assumption
    for( auto it=DT._E.begin(); it!=DT._E.end(); ++it ) // for all v_i, 1<=i<=n
        assert( it->second.size() <= 3 );         // |v_i| <= 3

    return SG;
}

namespace spanning_graph {

    template< class T >
    void add_first_edge( const T& G, typename T::Vertex_handle v, typename T::Vertex_circulator C ) {
        typename T::Vertex_handle v2 = C->handle();
        G.add_edge( v, v2 );
    }

    template< class T >
    void add_second_edge( const T& G, typename T::Vertex_handle v, typename T::Vertex_circulator C ) {
        while( G._DT.is_infinite(++C) );
        typename T::Vertex_handle v2 = C->handle();
        G.add_edge( v, v2 );
    }

    template< class T >
    void add_last_edge( const T& G, typename T::Vertex_handle v, typename T::Vertex_circulator C, typename T::VertexSet is_removed ) {
        --C;
        typename T::Vertex_circulator done(C);

        //while( ( C->info().is_removed || G._DT->is_infinite(C) ) && --C != done );
        while( ( contains( is_removed, C ) || G._DT.is_infinite(C) ) && --C != done );

        typename T::Vertex_handle v2 = C->handle();

        G.add_edge( v, v2 );
    }

    template< class T >
    void remove_first_edge( const T& G, typename T::Vertex_circulator C ) {
        typename T::Vertex_handle v1 = C->handle(),
                                  v2 = (++C)->handle();

        G.remove_edge( v1, v2 );
    }

    template< class T >
    void remove_second_edge( const T& G, typename T::Vertex_circulator C ) {
        typename T::Vertex_handle v1 = (++C)->handle(),
                                  v2 = (++C)->handle();

        G.remove_edge( v1, v2 );
    }

    template< class T >
    void remove_last_edge( const T& G, typename T::Vertex_circulator C, typename T::VertexSet is_removed ) {
        --C;

        typename T::Vertex_circulator done(C);

        //while( ( C->info().is_removed || G._DT->is_infinite(C) ) && --C != done );
        while( ( contains( is_removed, C ) || G._DT.is_infinite(C) ) && --C != done );

        typename T::Vertex_handle v1 = C->handle(),
                                  v2 = (--C)->handle();

        G.remove_edge( v1, v2 );
    }

}; // namespace spanning_graph

} // namespace gsnunf

#endif // GSNUNF_SPANNINGGRAPH_H
