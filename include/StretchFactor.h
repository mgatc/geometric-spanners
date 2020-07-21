#ifndef GSNUNF_STRETCHFACTOR_H
#define GSNUNF_STRETCHFACTOR_H

#include "DelaunayGraph.h"
#include "FloydWarshall.h"

namespace gsnunf {

template< class DG >
inline size_t getN( const DG& G ) {
    return G._DT.number_of_vertices();
}

template< class DG >
inline typename DG::FT getDistance( const DG& G, const typename DG::Vertex_handle a, const typename DG::Vertex_handle b ) {
    return a == b ? 0 : CGAL::sqrt( CGAL::squared_distance( a->point(), b->point() ) );
}

template< class DG >
void createVertexToIndexMap( const DG& G, typename DG::template VertexMap<size_t>& index ) {
    index.clear();
    index.reserve( getN(G) );
    size_t i=0;
    for( auto it = G._DT.finite_vertices_begin();
        it != G._DT.finite_vertices_end();
        ++it
    ) index.emplace( it, i++ );
}

template< class DG >
void EuclideanDistanceMatrix( const DG& G, const typename DG::template VertexMap<size_t>& index, vector< vector< optional<typename DG::FT> > >& euclidean ) {
    using Vertex_handle = typename DG::Vertex_handle;
    size_t N = getN(G);

    // Create an NxN table to hold distances.
    vector< vector< optional<typename DG::FT> > > eucl( N, vector< optional<typename DG::FT> >(N, nullopt) );

    for( auto i = G._DT.finite_vertices_begin(); i != G._DT.finite_vertices_end(); ++i )
        for( auto j = G._DT.finite_vertices_begin(); j != G._DT.finite_vertices_end(); ++j )
            eucl.at(index.at(i)).at(index.at(j)) =
                make_optional(
                    getDistance( G, i->handle(), j->handle() )
                );

    // Make sure we added distances for all pairs, none should be nullopt (infinite)
    for( size_t i=0; i<N; ++i )
        for( size_t j=0; j<N; ++j )
            assert( eucl.at(i).at(j) );

    swap( eucl, euclidean );

    return;
}

template< class DG >
void StretchFactorMatrix( const DG& G, vector< vector< optional<typename DG::FT> > >& stretch ) {
    using Vertex_handle = typename DG::Vertex_handle;
    size_t N = getN(G);

    // First, create a vertex-to-index map
    // Add all vertices to a vertex map and assign an index

    typename DG::template VertexMap< size_t > index;
    createVertexToIndexMap( G, index );

    // Next, conduct Floyd-Warshall to determine all paths' cost
    FloydWarshall( G, index, stretch );
    assert( stretch.size() == N );
    // Next, determine Euclidean distance between all vertices
    vector< vector< optional<typename DG::FT> > > euclidean;
    EuclideanDistanceMatrix( G, index, euclidean );
    assert( euclidean.size() == N );

    vector< vector< optional<typename DG::FT> > > quotient( N, vector< optional<typename DG::FT> >(N) );

    for( size_t i=0; i<N; ++i )
        for( size_t j=0; j<N; ++j )
            quotient.at(i).at(j) =
                stretch.at(i).at(j) ?
                    make_optional( i==j ? 0 : *stretch.at(i).at(j) / *euclidean.at(i).at(j) )
                    : nullopt;

    swap( quotient, stretch );

    return;
}

template< class DG >
typename DG::FT StretchFactor( const DG& G ) {
    vector< vector< optional<typename DG::FT> > > stretch;
    StretchFactorMatrix( G, stretch );

    typename DG::FT maxValue = 1.0;
    // Find max in stretch
    for( auto i : stretch )
        for( auto j : i ) {
            maxValue = CGAL::max( *j, maxValue );
            //cout<<*j<<"\n";
        }

    return maxValue;
}

} // namespace gsnunf


#endif // GSNUNF_STRETCHFACTOR_H

