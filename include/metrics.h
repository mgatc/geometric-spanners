#ifndef GSNUNF_METRICS_H
#define GSNUNF_METRICS_H

#include <algorithm> // swap
#include <optional>
#include <vector>

#include <boost/functional/hash.hpp>

#include "DelaunayGraph.h"
#include "FloydWarshall.h"

namespace gsnunf {

using namespace std;

template< typename T1, typename T2, typename F >
void forBoth( const std::pair<T1,T2>& p, F func ) {
    func( p.first, p.second );
    func( p.second, p.first );
}

namespace metrics {

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
        indexToHandle.emplace_back( it );
        handleToIndex.emplace( it, i );
        ++i;
    }
}

} // namespace metrics

void EuclideanDistanceMatrix( const DelaunayGraph& G, const DelaunayGraph::template VertexMap<size_t>& index, vector< vector< optional<DelaunayGraph::FT> > >& euclidean ) {
    using namespace metrics;
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

using StretchFactorIndexEntry = pair<pair<size_t,size_t>, DelaunayGraph::FT>;
using StretchFactorVertexHandleEntry = pair<pair<Vertex_handle,Vertex_handle>, DelaunayGraph::FT>;

StretchFactorVertexHandleEntry StretchFactor( const DelaunayGraph& G ) {
    using namespace metrics;
    vector< vector< optional<DelaunayGraph::FT> > > stretch;
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

struct PointHasher {
    std::size_t operator()(const Point& p) const noexcept {
        size_t seed = 31;
        boost::hash_combine( seed, p.x() );
        boost::hash_combine( seed, p.y() );
        return seed;
    }
};

inline size_t countIncident( std::unordered_map< Point,size_t,PointHasher >& count, const Point& p ) {
    if( !contains( count, p ) )
        count.emplace( p, 1 );
    else
        count[p]++;

    return count.at(p);
}

template< typename RandomAccessIterator >
size_t degree( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
    std::vector<pair<Point,Point>> edges( edgesBegin, edgesEnd );
    std::unordered_map< Point,size_t,PointHasher > count( edges.size() );
    size_t max = 0;
    // for each edge
    for( auto e : edges ) {
        max = CGAL::max( max, countIncident( count, e.first  ) );
        max = CGAL::max( max, countIncident( count, e.second ) );
    }
    return max/2;
}

template< typename Triangulation >
size_t degree( const Triangulation& T ) {
    // fill a vector with edges so we can call the range-based degree function
    std::vector<pair<Point,Point>> edges;
    edges.reserve( T.number_of_vertices() );

    for( auto e=T.finite_edges_begin(); e!=T.finite_edges_end(); ++e ) {
        auto p = make_pair(
            e->first->vertex( (e->second+1)%3 )->point(),
            e->first->vertex( (e->second+2)%3 )->point()
        );
        // Add both in and out edges
        forBoth( p, [&]( Point a, Point b ) {
            edges.emplace_back(a,b);
        });
    }
    return degree( edges.begin(), edges.end() );
}

template< typename RandomAccessIterator >
size_t weight( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
    double w = 0;
    for( auto e=edgesBegin; e!=edgesEnd; ++e ) {
        w += CGAL::sqrt( CGAL::squared_distance( e->first, e->second ) );
    }
    return w;
}
template< typename Triangulation >
size_t weight( const Triangulation& T ) {
    double w = 0;
    for( auto e=T.finite_edges_begin(); e!=T.finite_edges_end(); ++e ) {
        auto p = make_pair(
            e->first->vertex( (e->second+1)%3 )->point(),
            e->first->vertex( (e->second+2)%3 )->point()
        );
        w += CGAL::sqrt( CGAL::squared_distance( p.first, p.second ) );
    }
    return w;
}

} // namespace gsnunf


#endif // GSNUNF_STRETCHFACTOR_H


