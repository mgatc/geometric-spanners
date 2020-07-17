#ifndef GSNUNF_STRETCHFACTOR_H
#define GSNUNF_STRETCHFACTOR_H

#include "DelaunayGraph.h"
#include "FloydWarshall.h"

namespace gsnunf {

template< class DG >
void EuclideanDistanceMap( const DG& G, typename DG::template EdgeInfoMap<typename DG::FT>& euclidean ) {
    using Vertex_handle = typename DG::Vertex_handle;
    size_t n = G._DT.number_of_vertices();

    typename DG::template VertexMap< optional< typename DG::FT > > allVertices;

    // Add all vertices to a vertex map with nullopts (infinity)
    for( auto i = G._DT.finite_vertices_begin();
        i != G._DT.finite_vertices_end();
        ++i
    ) allVertices.emplace( i, std::nullopt );

    // Then, create an entry to map each vertex to every other vertex
    euclidean.clear();
    for( auto i = G._DT.finite_vertices_begin();
        i != G._DT.finite_vertices_end();
        ++i
    ) euclidean.emplace( i, allVertices );

    // Now, loop through everything and calculate distances
    for( auto& i : euclidean ) {
        for( auto& j : i.second ) {
            j.second = CGAL::sqrt( CGAL::squared_distance( i.first->point(), j.first->point() ) );
        }
    }
    return;
}

template< class DG >
void StretchFactorMap( const DG& G, typename DG::template EdgeInfoMap<typename DG::FT>& stretch ) {
    using Vertex_handle = typename DG::Vertex_handle;
    using EdgeInfoMap = typename DG::template EdgeInfoMap< typename DG::FT >;

    // First, conduct Floyd-Warshall to determine all paths' cost
    FloydWarshall( G, stretch );

    // Next, determine Euclidean distance between all vertices
    EdgeInfoMap euclidean;
    EuclideanDistanceMap( G, euclidean );

    for( auto& i : stretch ) {
        Vertex_handle u = i.first;
        for( auto& j : i.second ) {
            Vertex_handle v = j.first;
            if( u == v )
                j.second = std::make_optional( typename DG::FT(0) );
            else
                j.second = *j.second / *euclidean[u][v];
        }
    }

    return;
}

template< class DG >
typename DG::FT StretchFactor( const DG& G ) {
    using EdgeInfoMap = typename DG::template EdgeInfoMap< typename DG::FT >;

    EdgeInfoMap stretch;
    StretchFactorMap( G, stretch );

    typename DG::FT maxValue;
    // Find max in stretch
    for( auto& i : stretch )
        for( auto& j : i.second )
            maxValue = CGAL::max( *j.second, maxValue );
    return maxValue;
}

} // namespace gsnunf


#endif // GSNUNF_STRETCHFACTOR_H

