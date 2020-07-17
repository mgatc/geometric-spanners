#ifndef GSNUNF_FLOYDWARSHALL_H
#define GSNUNF_FLOYDWARSHALL_H

#include <algorithm>
#include <optional>
#include <vector>

#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/utils.h>

#include "DelaunayGraph.h"

namespace gsnunf {

//using nullopt = std::nullopt;
template< typename T >
using optional = std::optional<T>;
template< typename T >
using vector = std::vector<T>;
template< typename T >
using OptionalRow = vector< optional<T> >;
template< typename T >
using Matrix = vector< OptionalRow<T> >;

template< typename N >
optional<N> min( const optional<N>& ij, const std::pair< optional<N>,optional<N> > ikj ) {
    optional<N> newPathLength = std::nullopt;
    if( ikj.first && ikj.second ) newPathLength = *ikj.first + *ikj.second;
    if( ij && newPathLength ) return CGAL::min( ij, newPathLength );
    if( ij ) return ij;
    return newPathLength;
}

template< class DG >
void FloydWarshall( const DG& G, typename DG::template EdgeInfoMap<typename DG::FT>& dist ) {
    using Vertex_handle = typename DG::Vertex_handle;
    size_t n = G._DT.number_of_vertices();

    dist.clear();
    typename DG::template VertexMap< optional< typename DG::FT > > allVertices;

    // Add all vertices to a vertex map with nullopts (infinity)
    for( auto i = G._DT.finite_vertices_begin();
        i != G._DT.finite_vertices_end();
        ++i
    ) allVertices.emplace( i, std::nullopt );

    // Then, create an entry to map each vertex to every other vertex
    for( auto i = G._DT.finite_vertices_begin();
        i != G._DT.finite_vertices_end();
        ++i
    ) dist.emplace( i, allVertices );

    // Set all i==j to 0
    for( auto& i : dist )
        i.second[i.first] = std::make_optional( typename DG::FT(0) );

    // Add distance of each edge (u,v) in G._E to dist[u][v]
    for( auto& adjacent : G._E ) {
        Vertex_handle u = adjacent.first;
        for( Vertex_handle v : adjacent.second ) {
            dist[u][v] = CGAL::sqrt( CGAL::squared_distance( u->point(), v->point() ) );
        }
    }

    // Check if going through k yields a shorter path from i to j
    for( auto k : dist )
        for( auto& i : dist )
            for( auto j : dist ) {
                i.second[j.first] = gsnunf::min(
                    i.second[j.first],
                  { i.second[k.first], k.second[j.first] }
                );
            }

    return;
}

} // namespace gsnunf


#endif // GSNUNF_FLOYDWARSHALL_H
