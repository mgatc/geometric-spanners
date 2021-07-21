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

namespace unf_spanners {

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

void FloydWarshall( const DelaunayGraph& G,
                    const VertexMap<size_t>& index,
                    vector< vector< optional<number_t> > >& distances ) {
    using namespace floyd_warshall;
    size_t N = G.n();

    // Create an NxN table to hold distances.
    vector< vector< optional<number_t> > > dist( N, vector< optional<number_t> >( N, nullopt ) );
    // container constructor should initialize optionals using default constructor, aka nullopt, aka infinity

    // Set all i==j to 0 (getDistance to self)
    for( size_t i=0; i<N; ++i )
        dist.at(i).at(i) = make_optional( number_t(0) );

    assert( index.size() == N );

    // Add getDistance of each edge (u,v) in G._E to dist[u][v]
    // using indices of u and v mapped in index
    for( const auto& adjacent : G.m_E ) {
        VertexHandle u = adjacent.first; // get vertex handle
        for( const VertexHandle &v : adjacent.second )
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
}

} // namespace unf_spanners


#endif // GSNUNF_FLOYDWARSHALL_H
