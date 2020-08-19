#ifndef GSNUNF_STRETCHFACTOR_H
#define GSNUNF_STRETCHFACTOR_H

#include <algorithm> // swap
#include <optional>
#include <vector>

#include "DelaunayGraph.h"
#include "FloydWarshall.h"

namespace gsnunf {

using namespace std;

namespace stretch_factor {

inline DelaunayGraph::FT getDistance( const DelaunayGraph::Vertex_handle a, const DelaunayGraph::Vertex_handle b ) {
    return a == b ? 0 : CGAL::sqrt( CGAL::squared_distance( a->point(), b->point() ) );
}

void createVertexToIndexMap( const DelaunayGraph& G, DelaunayGraph::template VertexMap<size_t>& index ) {
    index.clear();
    index.reserve( G.n() );
    size_t i=0;
    for( auto it = G._DT.finite_vertices_begin();
        it != G._DT.finite_vertices_end();
        ++it
    ) index.emplace( it, i++ );
}

} // namespace stretch_factor

void EuclideanDistanceMatrix( const DelaunayGraph& G, const DelaunayGraph::template VertexMap<size_t>& index, vector< vector< optional<DelaunayGraph::FT> > >& euclidean ) {
    using namespace stretch_factor;
    size_t N = G.n();

    // Create an NxN table to hold distances.
    vector< vector< optional<DelaunayGraph::FT> > > eucl( N, vector< optional<DelaunayGraph::FT> >(N, nullopt) );

    for( auto i = G._DT.finite_vertices_begin(); i != G._DT.finite_vertices_end(); ++i )
        for( auto j = G._DT.finite_vertices_begin(); j != G._DT.finite_vertices_end(); ++j )
            eucl.at(index.at(i)).at(index.at(j)) =
                make_optional(
                    getDistance( i->handle(), j->handle() )
                );

    // Make sure we added distances for all pairs, none should be nullopt (infinite)
    for( size_t i=0; i<N; ++i )
        for( size_t j=0; j<N; ++j )
            assert( eucl.at(i).at(j) );

    swap( eucl, euclidean );

    return;
}

void StretchFactorMatrix( const DelaunayGraph& G, vector< vector< optional<DelaunayGraph::FT> > >& stretch ) {
    using namespace stretch_factor;

    size_t N = G.n();

    // First, create a vertex-to-index map
    // Add all vertices to a vertex map and assign an index

    DelaunayGraph::template VertexMap< size_t > index;
    createVertexToIndexMap( G, index );

    // Next, conduct Floyd-Warshall to determine all paths' cost
    FloydWarshall( G, index, stretch );
    // Next, determine Euclidean distance between all vertices
    vector< vector< optional<DelaunayGraph::FT> > > euclidean;
    EuclideanDistanceMatrix( G, index, euclidean );

    vector< vector< optional<DelaunayGraph::FT> > > quotient( N, vector< optional<DelaunayGraph::FT> >(N) );

    for( size_t i=0; i<N; ++i )
        for( size_t j=0; j<N; ++j )
            quotient.at(i).at(j) =
                stretch.at(i).at(j) ?
                    make_optional( i==j ? 0 : *stretch.at(i).at(j) / *euclidean.at(i).at(j) )
                    : nullopt;

    swap( quotient, stretch );

    return;
}

DelaunayGraph::FT StretchFactor( const DelaunayGraph& G ) {
    using namespace stretch_factor;
    vector< vector< optional<DelaunayGraph::FT> > > stretch;
    StretchFactorMatrix( G, stretch );

    DelaunayGraph::FT maxValue = 1.0;
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

