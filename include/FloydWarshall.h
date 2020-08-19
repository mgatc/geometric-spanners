#ifndef GSNUNF_FLOYDWARSHALL_H
#define GSNUNF_FLOYDWARSHALL_H

#include <algorithm> // swap
#include <iostream>
#include <optional>
#include <utility> // pair
#include <vector>

#include <CGAL/number_utils.h> // min
#include <CGAL/squared_distance_2.h>
#include <CGAL/utils.h>

#include "DelaunayGraph.h"

namespace gsnunf {

using namespace std;

namespace floyd_warshall {

template< typename N >
optional<N> min( const optional<N>& ij, const pair< optional<N>,optional<N> >& ikj ) {
//    cout<<"min start"<<endl;

    optional<N> newPathLength = nullopt;

    if( ikj.first && ikj.second ) newPathLength = { *ikj.first + *ikj.second };
    if( ij && newPathLength ) return { CGAL::min( *ij, *newPathLength ) };
    if( ij ) return ij;

    return newPathLength;
}

} // namespace floyd_warshall

void FloydWarshall( const DelaunayGraph& G, const DelaunayGraph::template VertexMap<size_t>& index, vector< vector< optional<DelaunayGraph::FT> > >& distances ) {
    using namespace floyd_warshall;
    using Vertex_handle = DelaunayGraph::Vertex_handle;
    size_t N = G.n();

    // Create an NxN table to hold distances.
    vector< vector< optional<DelaunayGraph::FT> > > dist( N, vector< optional<DelaunayGraph::FT> >( N, nullopt ) );
    // container constructor should initialize optionals using default constructor, aka nullopt, aka infinity

    // Set all i==j to 0 (distance to self)
    for( size_t i=0; i<N; ++i )
        dist.at(i).at(i) = make_optional( DelaunayGraph::FT(0) );

    assert( index.size() == N );

    // Add distance of each edge (u,v) in G._E to dist[u][v]
    // using indices of u and v mapped in index
    for( const auto& adjacent : G._E ) {
        Vertex_handle u = adjacent.first; // get vertex handle
        for( const Vertex_handle v : adjacent.second )
            dist.at(index.at(u)).at(index.at(v)) = make_optional( CGAL::sqrt( CGAL::squared_distance( u->point(), v->point() ) ) );

    }

    // Check if going through k yields a shorter path from i to j
    for( size_t k=0; k<N; ++k )
        for( size_t i=0; i<N; ++i )
            for( size_t j=0; j<N; ++j )
                dist.at(i).at(j) = floyd_warshall::min(
                    dist.at(i).at(j),
                  { dist.at(i).at(k), dist.at(k).at(j) }
                );

    // swap the addresses for array we built with the address given in parameters
    swap( dist, distances );

    return;
}

} // namespace gsnunf


#endif // GSNUNF_FLOYDWARSHALL_H
