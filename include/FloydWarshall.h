#ifndef GSNUNF_FLOYDWARSHALL_H
#define GSNUNF_FLOYDWARSHALL_H

#include <algorithm>
#include <array>
#include <optional>

#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/utils.h>

#include "DelaunayGraph.h"

namespace gsnunf {

using namespace std;

template< typename N >
optional<N> min( const optional<N>& ij, const pair< optional<N>,optional<N> >& ikj ) {
//    cout<<"min start"<<endl;

    optional<N> newPathLength = nullopt;

    if( ikj.first && ikj.second ) newPathLength = { *ikj.first + *ikj.second };
    if( ij && newPathLength ) return { CGAL::min( *ij, *newPathLength ) };
    if( ij ) return ij;

    return newPathLength;
}

template< class DG >
void FloydWarshall( const DG& G, const typename DG::template VertexMap<size_t>& index, vector< vector< optional<typename DG::FT> > >& distances ) {
    using Vertex_handle = typename DG::Vertex_handle;

    size_t N = getN(G);

    // Then, create an NxN table to hold distances.
    vector< vector< optional<typename DG::FT> > > dist( N, vector< optional<typename DG::FT> >( N, nullopt ) );
    // container constructor should initialize optionals using default constructor, aka nullopt

    // Set all i==j to 0 (distance to self)
    for( size_t i=0; i<N; ++i )
        dist.at(i).at(i) = make_optional( typename DG::FT(0) );

    // Add distance of each edge (u,v) in G._E to dist[u][v]
    // using indices of u and v mapped in index
    for( auto& adjacent : G._E ) {
        Vertex_handle u = adjacent.first; // get vertex handle
        for( Vertex_handle v : adjacent.second ) {
            dist.at(index.at(u)).at(index.at(v)) = make_optional( CGAL::sqrt( CGAL::squared_distance( u->point(), v->point() ) ) );
        }
    }

    // Check if going through k yields a shorter path from i to j
    for( size_t k=0; k<N; ++k ) {
        for( size_t i=0; i<N; ++i ) {
            for( size_t j=0; j<N; ++j ) {
                //cout<<"1:"<<*dist.at(i).at(j)<<" 2:"<<*dist.at(i).at(k)<<" "<< *dist.at(k).at(j)<<"\n";
                dist.at(i).at(j) = gsnunf::min(
                    dist.at(i).at(j),
                  { dist.at(i).at(k), dist.at(k).at(j) }
                );
            }
        }
    }
    // swap the addresses for array we built with the address given in parameters
    swap( dist, distances );

    return;
}

} // namespace gsnunf


#endif // GSNUNF_FLOYDWARSHALL_H
