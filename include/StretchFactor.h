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

void createVertexToIndexMaps( const DelaunayGraph& G, DelaunayGraph::template VertexMap<size_t>& handleToIndex, vector<Vertex_handle>& indexToHandle ) {
    handleToIndex.clear();
    indexToHandle.clear();
    handleToIndex.reserve( G.n() );
    indexToHandle.reserve( G.n() );
    size_t i=0;
    for( auto it = G._DT.finite_vertices_begin();
        it != G._DT.finite_vertices_end();
        ++it ) {
        handleToIndex.emplace( it, i++ );
        indexToHandle.emplace_back( it );
    }
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

//void StretchFactorMatrix( const DelaunayGraph& G, vector< vector< optional<DelaunayGraph::FT> > >& stretch ) {
//    using namespace stretch_factor;
//
//    size_t N = G.n();
//
//    // First, create a vertex-to-index map
//    // Add all vertices to a vertex map and assign an index
//
//    DelaunayGraph::template VertexMap< size_t > index;
//    createVertexToIndexMap( G, index );
//
//    // Next, conduct Floyd-Warshall to determine all paths' cost
//    FloydWarshall( G, index, stretch );
//    // Next, determine Euclidean distance between all vertices
//    vector< vector< optional<DelaunayGraph::FT> > > euclidean;
//    EuclideanDistanceMatrix( G, index, euclidean );
//
//    vector< vector< optional<DelaunayGraph::FT> > > quotient( N, vector< optional<DelaunayGraph::FT> >(N) );
//
//    for( size_t i=0; i<N; ++i )
//        for( size_t j=0; j<N; ++j )
//            quotient.at(i).at(j) =
//                stretch.at(i).at(j) ?
//                    make_optional( i==j ? 0 : *stretch.at(i).at(j) / *euclidean.at(i).at(j) )
//                    : nullopt;
//
//    swap( quotient, stretch );
//
//    return;
//}

using StretchFactorIndexEntry = pair<pair<       size_t,       size_t>, DelaunayGraph::FT>;
using StretchFactorVertexHandleEntry = pair<pair<Vertex_handle,Vertex_handle>, DelaunayGraph::FT>;

StretchFactorVertexHandleEntry StretchFactor( const DelaunayGraph& G ) {
    using namespace stretch_factor;
    vector< vector< optional<DelaunayGraph::FT> > > stretch;
    //StretchFactorMatrix( G, stretch );
    size_t N = G.n();

    // First, create a vertex-to-index map
    // Add all vertices to a vertex map and assign an index

    DelaunayGraph::template VertexMap< size_t > handleToIndex;
    vector<Vertex_handle> indexToHandle;
    createVertexToIndexMaps( G, handleToIndex, indexToHandle );

    // Next, conduct Floyd-Warshall to determine all paths' cost
    FloydWarshall( G, handleToIndex, stretch );
    // Next, determine Euclidean distance between all vertices
    vector< vector< optional<DelaunayGraph::FT> > > euclidean;
    EuclideanDistanceMatrix( G, handleToIndex, euclidean );

    vector< vector< optional<DelaunayGraph::FT> > > quotient( N, vector< optional<DelaunayGraph::FT> >(N) );

    for( size_t i=0; i<N; ++i )
        for( size_t j=0; j<N; ++j )
            quotient.at(i).at(j) =
                stretch.at(i).at(j) ?
                    make_optional( i==j ? 0 : *stretch.at(i).at(j) / *euclidean.at(i).at(j) )
                    : nullopt;

    swap( quotient, stretch );

    StretchFactorIndexEntry maxVal = make_pair( make_pair(0, 0), 1.0 );
    // Find max in stretch
    for( size_t i=0; i<stretch.size(); ++i ) {
        for( size_t j=0; j<stretch.at(i).size(); ++j ) {
            if( stretch.at(i).at(j) > maxVal.second ) {
                maxVal = make_pair( make_pair( i, j ), *stretch.at(i).at(j) );
            }
        }
    }
    return make_pair(
        make_pair(
            indexToHandle.at( maxVal.first.first  ),
            indexToHandle.at( maxVal.first.second )
        ), maxVal.second
    );
}

template< typename RandomAccessIterator >
StretchFactorVertexHandleEntry StretchFactor( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
    DelaunayGraph G;
    G.buildFromEdgeList( edgesBegin, edgesEnd );
    return StretchFactor(G);
}

} // namespace gsnunf


#endif // GSNUNF_STRETCHFACTOR_H

