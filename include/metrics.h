#ifndef GSNUNF_METRICS_H
#define GSNUNF_METRICS_H

#include <algorithm> // swap
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/squared_distance_2.h> //for 2D functions

#include <omp.h>

#include "DelaunayGraph.h"
#include "FloydWarshall.h"
#include "utilities.h"

namespace gsnunf {

using namespace std;

template< typename Point >
struct PointHasher {
    std::size_t operator()(const Point& p) const noexcept {
        size_t seed = 31;
        boost::hash_combine( seed, p.x() );
        boost::hash_combine( seed, p.y() );
        return seed;
    }
};

template< typename T1, typename T2, typename F >
void forBoth( const std::pair<T1,T2>& p, F func ) {
    func( p.first, p.second );
    func( p.second, p.first );
}

namespace metrics {

//inline number_t getDistance( const Vertex_handle a, const Vertex_handle b ) {
//    return a == b ? 0 : CGAL::sqrt( CGAL::squared_distance( a->point(), b->point() ) );
//}


void createVertexToIndexMaps( const DelaunayGraph& G,
                              VertexMap<size_t>& handleToIndex,
                              vector<Vertex_handle>& indexToHandle ) {
    handleToIndex.clear();
    indexToHandle.clear();
    handleToIndex.reserve( G.n() );
    indexToHandle.reserve( G.n() );
    size_t i=0;
    for( auto it = G.m_DT.finite_vertices_begin();
        it != G.m_DT.finite_vertices_end();
        ++it ) {
        indexToHandle.emplace_back( it );
        handleToIndex.emplace( it, i );
        ++i;
    }
}

} // namespace metrics

void EuclideanDistanceMatrix( const DelaunayGraph& G,
                              const VertexMap<size_t>& index,
                              vector< vector< optional<number_t> > >& euclidean ) {
    using namespace metrics;
    size_t N = G.n();

    // Create an NxN table to hold distances.
    vector< vector< optional<number_t> > > eucl( N, vector< optional<number_t> >(N, nullopt) );

    for( auto i = G.m_DT.finite_vertices_begin(); i != G.m_DT.finite_vertices_end(); ++i )
        for( auto j = G.m_DT.finite_vertices_begin(); j != G.m_DT.finite_vertices_end(); ++j )
            eucl.at(index.at(i)).at(index.at(j)) =
                make_optional(
                    distance( i->point(), j->point() )
                );

    // Make sure we added distances for all pairs, none should be nullopt (infinite)
    for( size_t i=0; i<N; ++i )
        for( size_t j=0; j<N; ++j )
            assert( eucl.at(i).at(j) );

    swap( eucl, euclidean );
}

using StretchFactorIndexEntry = pair<pair<size_t,size_t>, number_t>;
using StretchFactorVertexHandleEntry = pair<pair<Vertex_handle,Vertex_handle>, number_t>;

StretchFactorVertexHandleEntry StretchFactorFloydWarshall( const DelaunayGraph& G ) {
    using namespace metrics;
    vector< vector< optional<number_t> > > stretch;
    size_t N = G.n();

    // First, create a vertex-to-index map
    // Add all vertices to a vertex map and assign an index

    VertexMap< size_t > handleToIndex;
    vector<Vertex_handle> indexToHandle;
    createVertexToIndexMaps( G, handleToIndex, indexToHandle );

    // Next, conduct Floyd-Warshall to determine all paths' cost
    FloydWarshall( G, handleToIndex, stretch );
    // Next, determine Euclidean distance between all vertices
    vector< vector< optional<number_t> > > euclidean;
    EuclideanDistanceMatrix( G, handleToIndex, euclidean );

    vector< vector< optional<number_t> > > quotient( N, vector< optional<number_t> >(N) );

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
StretchFactorVertexHandleEntry StretchFactorFloydWarshall( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
    DelaunayGraph G;
    G.buildFromEdgeList( edgesBegin, edgesEnd );
    return StretchFactorFloydWarshall(G);
}

template< typename Point >
inline size_t countIncident( std::unordered_map< Point,size_t,PointHasher<Point> >& count, const Point& p ) {
    if( !contains( count, p ) )
        count.emplace( p, 1 );
    else
        count[p]++;

    return count.at(p);
}

template< typename RandomAccessIterator >
size_t degree( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
    typedef typename RandomAccessIterator::value_type Edge;
    typedef typename Edge::first_type Point_2;

    std::vector<pair<Point_2,Point_2>> edges( edgesBegin, edgesEnd );
    std::unordered_map<Point_2,unordered_set<Point_2>> adj;
    // for each edge
    for( auto e : edges ) {
        auto first = adj.begin();
        tie(first,ignore) = adj.emplace( e.first, unordered_set<Point_2>() );
        (*first).second.insert(e.second);

        auto second = adj.begin();
        tie(second,ignore) = adj.emplace( e.second, unordered_set<Point_2>() );
        (*second).second.insert(e.first);
    }
    auto max_el = max_element( adj.begin(), adj.end(), [&] ( const auto& lhs, const auto& rhs ) {
        return lhs.second.size() < rhs.second.size();
    });

    return max_el->second.size();
}

template< typename Triangulation >
size_t degree( const Triangulation& T ) {
    typedef typename Triangulation::Point
        Point_2;

    // fill a vector with edges so we can call the range-based degree function
    std::vector<pair<Point_2,Point_2>> edges;
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
number_t weight( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
    number_t w = 0;
    for( auto e=edgesBegin; e!=edgesEnd; ++e ) {
        w += distance( e->first, e->second );
    }
    return w;
}
template< typename Triangulation >
number_t weight( const Triangulation& T ) {
    number_t w = 0;
    for( auto e=T.finite_edges_begin(); e!=T.finite_edges_end(); ++e ) {
        auto p = make_pair(
            e->first->vertex( (e->second+1)%3 )->point(),
            e->first->vertex( (e->second+2)%3 )->point()
        );
        w += distance( p.first, p.second );
    }
    return w;
}

template< typename Point_2 >
struct EuclideanDistanceToPoint {
    Point_2 goal;
    double operator()( Point_2 p ) {
        return distance(p,goal);
    }
};

template<typename T>
struct MinHeapCompare {
    bool operator()( const T &n1, const T &n2 ) const {
        return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
    }
};

template<typename T>
struct MaxHeapCompare {
    bool operator()( const T &n1, const T &n2 ) const {
        return (n1.first < n2.first) || ((n1.first == n2.first) && (n1.second < n2.second));
    }
};

template< typename VertexContainer, typename VertexMap, typename AdjacencyList >
optional<number_t> AStar( VertexContainer V, VertexMap vMap, AdjacencyList G_prime, size_t start, size_t goal ) {
    typedef pair<number_t,size_t>
        DistanceIndexPair;
    typedef boost::heap::fibonacci_heap< DistanceIndexPair,boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
        Heap;
    typedef Heap::handle_type
        HeapHandle;

    size_t n = V.size();
    auto startPoint = V.at(start)->point();
    auto goalPoint = V.at(goal)->point();
    EuclideanDistanceToPoint h = { V.at(goal)->point() }; // initialize heuristic functor

    Heap open;
    unordered_map<size_t,HeapHandle> handleToHeap(n);
    handleToHeap[start] = open.emplace( h( startPoint ), start );

    //unordered_set<size_t> closed(n);
    vector<size_t> parents(n);

    vector<number_t> g( n, INF );
    g[start] = 0;

    vector<number_t> f( n, INF );
    f[start] = h( startPoint );

    DistanceIndexPair current = open.top(); // initialize current vertex to start
    size_t u_index = current.second;
    auto currentPoint = startPoint;
    auto neighborPoint = currentPoint;
//    cout<<"\n    A* start:"<<startPoint;
//    cout<<",";
//    cout<<" goal:"<<V.at(goal)->point();
//    cout<<",";

    do {
        current = open.top();
        open.pop();

        u_index = current.second;
        currentPoint = V.at(u_index)->point();
//        cout<<"\n      current:"<<currentPoint;
//        cout<<",";
        if( u_index == goal ) return make_optional( g.at(goal) );
//        cout<<" no goal, ";
        // loop through neighbors of current
        for( size_t neighbor : G_prime.at(u_index) ) {
            neighborPoint = V.at(neighbor)->point();
//            cout<<"\n        n:"<<neighborPoint;
//            cout<<",";
            number_t newScore = g.at(u_index)
                + d( currentPoint, neighborPoint );
//            cout<<"g:"<<newScore;
//            cout<<",";
            if( newScore < g.at( neighbor ) ) {
                parents[neighbor] = u_index;
                g[neighbor] = newScore;
                f[neighbor] = g.at(neighbor) + h(neighborPoint);
                DistanceIndexPair q = make_pair( f.at(neighbor), neighbor );

                if( contains( handleToHeap, neighbor ) ) {
                    HeapHandle neighborHandle = handleToHeap.at(neighbor);
                    open.update(neighborHandle,q);
                    open.update(neighborHandle);
                } else {
                    handleToHeap[neighbor] = open.push(q);
                }
            }
        }
    } while( !open.empty() );

    return nullopt;
}



template< typename VertexContainer, typename AdjacencyList >
void Dijkstra( const size_t i, const VertexContainer& V, const AdjacencyList& G, vector<number_t>& ShortestPaths, vector<size_t>& Parents ) {

//    typedef pair<double,size_t>
//        DistanceIndexPair;
    //typedef boost::heap::fibonacci_heap< DistanceIndexPair,boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
    typedef std::map<number_t,size_t>
        Heap;
    typedef Heap::iterator
        HeapHandle;

    const size_t n = V.size();
    auto startPoint = V.at(i);

    Heap open;
    unordered_map<size_t,HeapHandle> handleToHeap(n);
    handleToHeap[i] = open.emplace( 0, i ).first;

    //unordered_set<size_t> closed(n);


    ShortestPaths[i] = 0;

    auto current = open.begin(); // initialize current vertex to start
    size_t u_index = current->second;
    auto currentPoint = startPoint;
    auto neighborPoint = currentPoint;
    number_t newScore = 0;
//    cout<<"\n    A* start:"<<startPoint;
//    cout<<",";
//    cout<<" goal:"<<V.at(goal)->point();
//    cout<<",";

    do {
        current = open.begin();

        u_index = current->second;
        currentPoint = V[u_index];
//        cout<<"\n      current:"<<currentPoint;
//        cout<<",";
        // loop through neighbors of current
        for( size_t neighbor : G.at(u_index) ) {
            neighborPoint = V[neighbor];
//            cout<<"\n        n:"<<neighborPoint;
//            cout<<",";
            newScore = ShortestPaths[u_index]
                + distance( currentPoint, neighborPoint );
//            cout<<"g:"<<newScore;
//            cout<<",";
            if( newScore < ShortestPaths[neighbor] ) {
                Parents[neighbor] = u_index;
                ShortestPaths[neighbor] = newScore;
                //DistanceIndexPair q = make_pair( ShortestPaths.at(neighbor), neighbor );
                if( contains( handleToHeap, neighbor ) ) {
                    open.erase(handleToHeap[neighbor]);
                }
                handleToHeap[neighbor] = open.emplace(ShortestPaths[neighbor], neighbor).first;
            }
        }
        //closed.insert(current.second);
        open.erase(current);
    } while( !open.empty() );
}

//template< typename RandomAccessIterator >
//double StretchFactorDijkstra( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
//    typedef typename RandomAccessIterator::value_type Edge;
//    typedef typename Edge::first_type Point;
//
//    vector<Point> V; // container for vertices
//    unordered_map< Point, size_t, PointHasher<Point> > vMap; // map point to index in V
//    unordered_map< size_t, unordered_set<size_t> > G; // adjacency list
//    size_t index = 0;
//
//    // Create list of vertices, map to their indices, and adjacency list
//    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) {
//        // If vMap doesn't contain p, put it in V
//        Point p = eit->first;
//        size_t i_p = index;
//        bool inserted = false;
//        auto vMapIt = vMap.begin();
//        tie( vMapIt, inserted ) = vMap.emplace( p, i_p ); // map for reverse lookup
//        if( inserted ) {
//            V.push_back(p);
//            ++index;
//        }
//        i_p = vMapIt->second;
//
//        // If vMap doesn't contain q, put it in V
//        Point q = eit->second;
//        size_t i_q = index;
//        tie( vMapIt, inserted ) = vMap.emplace( q, i_q ); // map for reverse lookup
//        if( inserted ) {
//            V.push_back(q);
//            ++index;
//        }
//        i_q = vMapIt->second;
//
//        G[i_p].insert(i_q); // add edge to adjacency list
//    }
//    size_t n = V.size();
//    vector<double> T( n, INF );
//
//    // calculate euclidean distance between all pairs
//    //#pragma omp parallel for
//    for( size_t i=0; i<n; ++i ) {
//        vector<double> ShortestPaths( n, INF );
//        vector<double> D( n, INF );     // Euclidean distances
//
//        for( size_t j=0; j<n; ++j ) {
//            D.at(j) =
//                i==j ? 0 : d( V.at(i), V.at(j) );
//        }
//
//        Dijkstra( i, V, G, ShortestPaths );
//
//        // Divide each shortest path distance by the euclidean distance between the vertices.
//        for( size_t j=0; j<n; ++j ) {
//            ShortestPaths.at(j) = ( // avoid /0
//                i==j ? 0 : ShortestPaths.at(j)/D.at(j)
//            );
//        }
//        // Find max t and place in T
//        T.at(i) = *max_element(
//            begin( ShortestPaths ),
//            end(   ShortestPaths )
//        );
//    }
//    // Find the big mac daddy t aka big money
//    return *max_element( T.begin(), T.end() );
//}
//
//template< typename RandomAccessIterator >
//double StretchFactorDijkstraParallel( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
//    typedef typename RandomAccessIterator::value_type Edge;
//    typedef typename Edge::first_type Point;
//
//    vector<Point> V; // container for vertices
//    unordered_map< Point, size_t, PointHasher<Point> > vMap; // map point to index in V
//    unordered_map< size_t, unordered_set<size_t> > G; // adjacency list
//    size_t index = 0;
//
//    // Create list of vertices, map to their indices, and adjacency list
//    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) {
//        // If vMap doesn't contain p, put it in V
//        Point p = eit->first;
//        size_t i_p = index;
//        bool inserted = false;
//        auto vMapIt = vMap.begin();
//        tie( vMapIt, inserted ) = vMap.emplace( p, i_p ); // map for reverse lookup
//        if( inserted ) {
//            V.push_back(p);
//            ++index;
//        }
//        i_p = vMapIt->second;
//
//        // If vMap doesn't contain q, put it in V
//        Point q = eit->second;
//        size_t i_q = index;
//        tie( vMapIt, inserted ) = vMap.emplace( q, i_q ); // map for reverse lookup
//        if( inserted ) {
//            V.push_back(q);
//            ++index;
//        }
//        i_q = vMapIt->second;
//
//        G[i_p].insert(i_q); // add edge to adjacency list
//    }
//    size_t n = V.size();
//    vector<double> T( n, INF );
//
//    // calculate euclidean distance between all pairs
//    #pragma omp parallel for
//    for( size_t i=0; i<n; ++i ) {
//        vector<double> ShortestPaths( n, INF );
//        vector<double> D( n, INF );     // Euclidean distances
//
//        for( size_t j=0; j<n; ++j ) {
//            D.at(j) =
//                i==j ? 0 : d( V.at(i), V.at(j) );
//        }
//
//        Dijkstra( i, V, G, D, ShortestPaths );
//
//    // Divide each shortest path distance by the euclidean distance between the vertices.
//        for( size_t j=0; j<n; ++j ) {
//            ShortestPaths.at(j) = ( // avoid /0
//                i==j ? 0 : ShortestPaths.at(j)/D.at(j)
//            );
//        }
//        // Find max t and place in T
//        T.at(i) = *max_element(
//            begin( ShortestPaths ),
//            end(   ShortestPaths )
//        );
//    }
//    // Find the big mac daddy t aka big money
//    return *max_element( T.begin(), T.end() );
//}
template< typename VertexIterator, typename EdgeIterator >
number_t StretchFactorDijkstraReduction( VertexIterator pointsBegin,
                                       VertexIterator pointsEnd,
                                       EdgeIterator edgesBegin,
                                       EdgeIterator edgesEnd )
{
    typedef typename VertexIterator::value_type Point_2;

    const vector<Point_2> V(pointsBegin, pointsEnd); // container for vertices
    const size_t n = V.size();
//    unordered_map< size_t, size_t > vMap; // map point to index in V
    vector< unordered_set<size_t> > G( n, unordered_set<size_t>() ); // adjacency list
    //size_t index = 0;

    // Create list of vertices, map to their indices, and adjacency list
    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit )
    {
        auto p = eit->first,
             q = eit->second;

        G[p].insert(q);
        G[q].insert(p);
    }
    //vector<double> T( n, INF );
    double t_max = 0.0;

    // calculate euclidean distance between all pairs
    #pragma omp parallel for reduction( max: t_max )
    for( size_t i=0; i<n; ++i ) {
        // Euclidean distances
        vector<number_t> D( n, INF );
        for( size_t j=0; j<n; ++j ) {
            D.at(j) =
                i==j ? 0 : distance( V.at(i), V.at(j) );
        }
        // Shortest paths
        vector<number_t> ShortestPaths( n, INF );
        vector<size_t> Parents(n);
        Dijkstra( i, V, G, ShortestPaths, Parents );

        // Divide each shortest path distance by the euclidean distance between the vertices.
        for( size_t j=0; j<n; ++j ) {
            ShortestPaths.at(j) = ( // avoid /0
                i==j ? 0 : ShortestPaths.at(j)/D.at(j)
            );
        }
        // Find max_t
        auto t_local = max_element(
            begin( ShortestPaths ),
            end(   ShortestPaths )
        );
        if( *t_local > t_max ) {
            t_max = *t_local;
        }
    }
    // Find the big mac daddy t aka big money
    return t_max;
}
template< typename K, typename RandomAccessIterator>
double StretchFactorDijkstraConvexHull( RandomAccessIterator edgesBegin,
                    RandomAccessIterator edgesEnd ) {
//    typedef typename RandomAccessIterator::value_type Edge;
    typedef typename K::Point_2 Point_2;

    vector<Point_2> V; // container for vertices
    unordered_map< Point_2, size_t, PointHasher<Point_2> > vMap; // map point to index in V
    unordered_map< size_t, unordered_set<size_t> > G; // adjacency list
    size_t index = 0;

    // Create list of vertices, map to their indices, and adjacency list
    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) {
        // If vMap doesn't contain p, put it in V
        Point_2 p = eit->first;
        size_t i_p = index;
        bool inserted = false;
        auto vMapIt = vMap.begin();
        tie( vMapIt, inserted ) = vMap.emplace( p, i_p ); // map for reverse lookup
        if( inserted ) {
            V.push_back(p);
            ++index;
        }
        i_p = vMapIt->second;

        // If vMap doesn't contain q, put it in V
        Point_2 q = eit->second;
        size_t i_q = index;
        tie( vMapIt, inserted ) = vMap.emplace( q, i_q ); // map for reverse lookup
        if( inserted ) {
            V.push_back(q);
            ++index;
        }
        i_q = vMapIt->second;

        G[i_p].insert(i_q); // add edge to adjacency list
    }
    size_t n = V.size();

    // Convex hull with indices
    typedef CGAL::Convex_hull_traits_adapter_2<K,
        typename CGAL::Pointer_property_map<Point_2>::type >
            Convex_hull_traits;

    std::vector<std::size_t> indices(n), ConvexHull;
    std::iota(indices.begin(), indices.end(),0);

    // Find the convex hull of the point set
    CGAL::convex_hull_2( indices.begin(), indices.end(), back_inserter(ConvexHull),
                         Convex_hull_traits(CGAL::make_property_map(V)) );

    //vector<double> T( n, INF );
    number_t t_max = 0.0;
    size_t i = 0;

    // Loop over the vertices of the convex hull
    //#pragma omp parallel for reduction( max: t_max )
    for( auto it=ConvexHull.begin(); it<ConvexHull.end(); ++it ) {
        i = *it;
        //cout<<endl<<endl<<i<<endl;
        // Euclidean distances
        vector<number_t> D( n, INF );
        for( size_t j=0; j<n; ++j ) {
            D.at(j) =
                i==j ? 0 : distance( V.at(i), V.at(j) );
        }
//        cout<<"  D"<<endl;
//        for( auto d : D ) {
//            cout<<"  "<<d<<endl;
//        }

        // Shortest paths
        vector<number_t> ShortestPaths( n, INF );
        vector<size_t> Parents(n);
        Dijkstra( i, V, G, ShortestPaths, Parents );

        // Divide each shortest path distance by the euclidean distance between the vertices.
        for( size_t j=0; j<n; ++j ) {
            ShortestPaths.at(j) = ( // avoid /0
                i==j ? 0 : ShortestPaths.at(j)/D.at(j)
            );
        }
//        cout<<"  ShortestPaths"<<endl;
//        for( auto d : ShortestPaths ) {
//            cout<<"  "<<d<<endl;
//        }
        // Find max_t
        auto t_local = max_element(
            begin( ShortestPaths ),
            end(   ShortestPaths )
        );
        if( *t_local > t_max ) {
            t_max = *t_local;
        }
    }
    // Find the big mac daddy t aka big money
    return t_max;
}

template< typename RandomAccessIterator, typename OutputIterator >
double SFWorstPath( RandomAccessIterator edgesBegin,
                    RandomAccessIterator edgesEnd,
                    std::optional<OutputIterator> out = std::nullopt ) {
    typedef typename RandomAccessIterator::value_type Edge;
    typedef typename Edge::first_type Point_2;

    vector<Point_2> V; // container for vertices
    unordered_map< Point_2, size_t, PointHasher<Point_2> > vMap; // map point to index in V
    unordered_map< size_t, unordered_set<size_t> > G; // adjacency list
    size_t index = 0;

    // Create list of vertices, map to their indices, and adjacency list
    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) {
        // If vMap doesn't contain p, put it in V
        Point_2 p = eit->first;
        size_t i_p = index;
        bool inserted = false;
        auto vMapIt = vMap.begin();
        tie( vMapIt, inserted ) = vMap.emplace( p, i_p ); // map for reverse lookup
        if( inserted ) {
            V.push_back(p);
            ++index;
        }
        i_p = vMapIt->second;

        // If vMap doesn't contain q, put it in V
        Point_2 q = eit->second;
        size_t i_q = index;
        tie( vMapIt, inserted ) = vMap.emplace( q, i_q ); // map for reverse lookup
        if( inserted ) {
            V.push_back(q);
            ++index;
        }
        i_q = vMapIt->second;

        G[i_p].insert(i_q); // add edge to adjacency list
    }
    size_t n = V.size();
    //vector<double> T( n, INF );
    number_t t_max = 0.0;
    vector<size_t> MaxParents;
    size_t i_max = 0, j_max=1;

    // calculate euclidean distance between all pairs
    //#pragma omp parallel for reduction( max: t_max )
    for( size_t i=0; i<n; ++i ) {
        // Euclidean distances
        vector<number_t> D( n, INF );
        for( size_t j=0; j<n; ++j ) {
            D.at(j) =
                i==j ? 0 : distance( V.at(i), V.at(j) );
        }
        // Shortest paths
        vector<number_t> ShortestPaths( n, INF );
        vector<size_t> Parents(n);
        Dijkstra( i, V, G, ShortestPaths, Parents );

        // Divide each shortest path distance by the euclidean distance between the vertices.
        for( size_t j=0; j<n; ++j ) {
            ShortestPaths.at(j) = ( // avoid /0
                i==j ? 0 : ShortestPaths.at(j)/D.at(j)
            );
        }
        // Find max_t
        auto t_local = max_element(
            begin( ShortestPaths ),
            end(   ShortestPaths )
        );
        if( *t_local > t_max ) {
            t_max = *t_local;
            // remove the following for parallel reduction function
            if(out){
                std::swap(Parents,MaxParents);
                i_max = i;
                j_max = t_local - ShortestPaths.begin();
            }
        }
    }
    if(out) {
        size_t walk = j_max;
        do {
            *(*out)++ = make_pair( V.at(walk), V.at(MaxParents.at(walk)) );
            walk = MaxParents.at(walk);
        } while( walk != i_max );
    }
    // Find the big mac daddy t aka big money
    return t_max;
}



template< typename VertexContainer, typename VertexMap, typename AdjacencyList >
optional<number_t> ShortestPath( VertexContainer V, VertexMap vMap, AdjacencyList G_prime, size_t start, size_t goal ) {
    return AStar( V, vMap, G_prime, start, goal );
}

// edgelists are range-supporting containers containing pairs of Points
template< typename RandomAccessIterator, typename Triangulation >
double StretchFactor( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd, const Triangulation& superGraph ) {
    using GraphVertex = typename Triangulation::Vertex_handle;
    using GraphCirculator = typename Triangulation::Vertex_circulator;

    typedef typename RandomAccessIterator::value_type Edge;
    typedef typename Edge::first_type Point_2;

    //using GPoint = typename Triangulation::Point_2;
    Triangulation G(superGraph);
    // create vector of points from given range and map points to indices
    vector<Point_2> P;
    P.reserve( G.number_of_vertices() );
    vector<GraphVertex> V;
    V.reserve( P.capacity() );
    unordered_map< GraphVertex, size_t, boost::hash<GraphVertex> > vMap( P.capacity() );
    size_t incidenceListInitialSize = 3;
    unordered_map< size_t, unordered_set<size_t> > G_prime( P.capacity() ); // E = O(n) so this is a decent enough guess
    size_t i = 0;

    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) { // Add all vertices. If insertion successful, put vertex in V and increment i.
        // first vertex in pair
        Point_2 p = eit->first;
        GraphVertex u = G.insert(p);
        bool inserted = false; // we need to know if the vertex is inserted
        auto it = vMap.begin();
        tie( it, inserted ) = vMap.emplace( u, i );
        size_t firstIndex = it->second;
        if( inserted ) { // only add and increment for unique vertices
            V.push_back(u);
            P.push_back(p);
            ++i;
        }
        // second vertex in pair
        Point_2 q = eit->second;
        GraphVertex v = G.insert(q);
        tie( it, inserted ) = vMap.emplace( v, i );
        size_t secondIndex = it->second;
        if( inserted ) { // only add and increment for unique vertices
            V.push_back(v);
            P.push_back(q);
            ++i;
        }
        assert( G.is_edge( u, v ) ); // ensure each edge in G_prime is in G
        G_prime[firstIndex].emplace( secondIndex ); // add edge to list
    }
    assert( G.number_of_vertices() == P.size() ); // ensure the insert functions above did not add a new vertex to the triangulation

    //double MAX_DOUBLE = numeric_limits<double>::max();
    double t_max = 0;

    //cout<<"V.size():"<<V.size()<<" G_prime.size():"<<G_prime.size()<<" ";

    // BFS containers
    vector<bool> status( V.size(), true ); // true == unknown

    for( GraphVertex u : V ) {
        // reset containers
        queue< pair<size_t,size_t> > level;
        fill( status.begin(), status.end(), true );
        double t_lvl = 0,
            t_v,
            euclideanDistance;
        optional<double> shortestPathLength = nullopt;
        size_t lvl=0;
        size_t u_index = vMap.at(u);
        // For testing, output info about u
//        cout<<"Linear:"<< u_index <<" ("<<u->point()<<") \n";
//        cout<<"  bfs\n";
        level.emplace( u_index, lvl); // start with u
        status[u_index] = false;
        pair<size_t,size_t> v;
        do { // BFS
            v = level.front();
            level.pop();

            // For testing, output info about v
//            cout<<"   n:"<<v.first<<" (lvl:"<<v.second<<", "<<P.at(v.first)<<")";
            // Check if we are done or on a new level
            if( v.second > lvl ) { // are we popping a new level from the queue?
                if( t_lvl < t_max ) { // yes, is this level's t less than the last level's?
                    break;
                }
                t_max = t_lvl;
//                cout<<" t_lvl:"<<t_lvl;
                ++lvl;
            }
            // process the vertex
            shortestPathLength = ShortestPath( V, vMap, G_prime, u_index, v.first );
            //cout<<" SPL:"<<*shortestPathLength<<" ";
            if( !shortestPathLength ) return numeric_limits<double>::max();
            Point_2 p_u = P.at(u_index ),
                  p_v = P.at( v.first );
            euclideanDistance = d( p_u, p_v );
            t_v = ( v.first != u_index ) ? // avoid divide-by-0
                *shortestPathLength / euclideanDistance : 0;
            t_lvl = CGAL::max( t_lvl, t_v );
//            cout<<*shortestPathLength<< "/"<<euclideanDistance<<"="<<t_v;
//            cout<<",";

            // BFS ADMIN
            // Add each neighbor in G to level
            //cout<<"     searching neighbors for BFS... ";
            GraphCirculator v_n = G.incident_vertices( V.at(v.first) ),
                done(v_n);
            do if( !G.is_infinite(v_n) ) {
                    size_t j = vMap.at(v_n);
                    //cout<<j<<" "<<P.at(j)<<" "<<(status.at(j) ? "unknown" : "known")<<"... ";
                    if( status.at(j) ) {
                        level.emplace(j, lvl + 1);
                        status[j] = false;
                    }
            } while( --v_n != done );
            // END BFS ADMIN
//            cout<<"\n";
        } while( !level.empty() );
//        cout<<"\n";
    }

    return t_max;
}

/*template< typename VertexContainer, typename VertexMap, typename AdjacencyList, typename Matrix, typename H>
void AStar( const VertexContainer& V, const VertexMap& vMap, AdjacencyList& G_prime, Matrix<double>& ShortestKnownPaths, const Matrix<double>& EuclideanDistances, Matrix<H>& upperBoundHandles, size_t start, size_t goal ) {
    typedef pair<double,size_t>
        DistanceIndexPair;
    typedef boost::heap::fibonacci_heap< DistanceIndexPair,boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
        Heap;
    typedef Heap::handle_type
        HeapHandle;

    size_t n = V.size();
    size_t inf = numeric_limits<size_t>::max();
    Point startPoint = V.at(start)->point();
    Point goalPoint = V.at(goal)->point();
    EuclideanDistance h = { V.at(goal)->point() }; // initialize heuristic functor

    Heap open;
    unordered_map<size_t,HeapHandle> handleToHeap(n);
    handleToHeap[start] = open.emplace( h( startPoint ), start );

    //unordered_set<size_t> closed(n);
    vector<size_t> parents(n);

    vector<double>& g = ShortestKnownPaths.at(i);

    vector<double> f( n, inf );
    f[start] = h( startPoint );

    DistanceIndexPair current = open.top(); // initialize current vertex to start
    size_t u_index = current.second;
    Point currentPoint = startPoint;
    Point neighborPoint;
//    cout<<"\n    A* start:"<<startPoint;
//    cout<<",";
//    cout<<" goal:"<<V.at(goal)->point();
//    cout<<",";

    do {
        current = open.top();
        open.pop();

        u_index = current.second;
        currentPoint = V.at(u_index)->point();
//        cout<<"\n      current:"<<currentPoint;
//        cout<<",";
        if( u_index == goal ) return;
//        cout<<" no goal, ";
        double t_new = 0;
        // loop through neighbors of current
        for( size_t neighbor : G_prime.at(u_index) ) {
            neighborPoint = V.at(neighbor)->point();
//            cout<<"\n        n:"<<neighborPoint;
//            cout<<",";
            double newScore = g.at(u_index)
                + d( currentPoint, neighborPoint );
//            cout<<"g:"<<newScore;
//            cout<<",";
            if( newScore < g.at( neighbor ) ) {
                parents[neighbor] = u_index;
                g[neighbor] = newScore;
                f[neighbor] = g.at(neighbor) + h(neighborPoint);
                DistanceIndexPair q = make_pair( f.at(neighbor), neighbor );

                // calculate the new path's t
                t_new = newScore / EuclideanDistances.at(i).at(u_index);
                // update t_upper in t-Heap
                auto tValue = make_pair( t_new, make_pair(i,u_index) );
                H tHandle = upperBoundHandles.at(i).at(u_index);


                if( contains( handleToHeap, neighbor ) ) {
                    HeapHandle neighborHandle = handleToHeap.at(neighbor);
                    open.update(neighborHandle,q);
                    open.update(neighborHandle);
                } else {
                    handleToHeap[neighbor] = open.push(q);
                }
            }
        }
    } while( !open.empty() );

    return;
}*/

// edgelists are range-supporting containers containing pairs of Points
//template< typename RandomAccessIterator >
//double StretchFactorExperimental( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
//    // First, parse input to structures that are convenient for our purposes
//    vector<Point> V; // container for vertices
//    unordered_map< Point, size_t, PointHasher > vMap; // map point to index in V
//    size_t index = 0;
//
//    // Create list of vertices and map to their indices
//    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) {
//        // If vMap doesn't contain p, put it in V
//        Point p = eit->first;
//        size_t i_p = index;
//        bool inserted = false;
//        auto vMapIt = vMap.begin();
//        tie( vMapIt, inserted ) = vMap.emplace( p, i_p ); // map for reverse lookup
//        if( inserted ) {
//            V.push_back(p);
//            ++index;
//        }
//        i_p = vMapIt->second;
//
//        // If vMap doesn't contain q, put it in V
//        Point q = eit->second;
//        size_t i_q = index;
//        tie( vMapIt, inserted ) = vMap.emplace( q, i_q ); // map for reverse lookup
//        if( inserted ) {
//            V.push_back(q);
//            ++index;
//        }
//        i_q = vMapIt->second;
//    }
//
//    // Fill adjacency list E
//    size_t n = V.size();
//    vector<unordered_set<size_t> > E(n); // adjacency list
//    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) {
//        size_t i_p = vMap.at( eit->first ),
//               i_q = vMap.at( eit->second );
//        E[i_p].insert(i_q); // add edge to adjacency list
//    }
//
//    const double INF = numeric_limits<double>::max();
//
//    // Fill euclidean distance matrix
//    vector< vector< double > > EuclideanDistances( n, vector<double>(n,INF) );
//    for( size_t i=0; i<n; ++i ) {
//        for( size_t j=0; j<n; ++j ) {
//            EuclideanDistances.at(i).at(j) = ( i==j ?
//                0 : d( V.at(i), V.at(j) )
//            );
//        }
//    }
//
//    // Step 1. Prepare the PQ (as a Fibonacci maxheap)
//    typedef pair< size_t, size_t >
//        StretchPair;
//    typedef pair< double, StretchPair >
//        tWithStretchPair;
//    typedef boost::heap::fibonacci_heap< tWithStretchPair,boost::heap::compare<MaxHeapCompare<tWithStretchPair>>>
//        UpperBoundHeap;
//    typedef UpperBoundHeap::handle_type
//        UpperBoundHeapHandle;
//    UpperBoundHeap upperBounds;
//
//    vector< vector<UpperBoundHeapHandle> > upperBoundHandles( n, vector<UpperBoundHeapHandle>(n) );
//
//    // Place all unique pairs that are not edges in E in upperBounds with t=inf
//    for( size_t i=0; i<n; ++i ) {
//        for( size_t j=0; j<n; ++j ) {
//            double t = INF;
//            t = (i==j ? 0 : t);
//            t = (contains( E.at(i), j ) ? 1 : t);
//            upperBoundHandles[i][j] = upperBounds.emplace( t, make_pair(i,j) );
//        }
//    }
//    // Step 2. Prepare the shortest known paths matrix
//    vector< vector< double > > ShortestKnownPaths(n, vector<double>(n, INF) );
//    for( size_t i=0; i<n; ++i ) { // Set self-loops to 0
//        ShortestKnownPaths[i][i] = 0;
//    }
////    for( size_t i=0; i<n; ++i ) { // Set edges to their euclidean distance
////        for( auto j : E.at(i) ) {
////            ShortestKnownPaths[i][j] = EuclideanDistances.at(i).at(j);
////        }
////    }
//
//    double t_currentPair = 1,
//           t_max = t_currentPair,
//           shortestKnownForCurrentPair = INF,
//           distance = INF;
//
//    tWithStretchPair currentPair = upperBounds.top();
//    size_t start = currentPair.second.first,
//           goal = currentPair.second.second;
//
//    while( upperBounds.top().first > t_max ) {
//        currentPair = upperBounds.top();
//        start = currentPair.second.first,
//        goal = currentPair.second.second;
//        cout<<"upperBounds.top():";
//        cout<<start<<","<<goal;
//        cout<<"("<<currentPair.first<<")";
//        cout<<",";
//
//        // START A*
//
//        typedef pair<double,size_t>
//            DistanceIndexPair;
//        typedef boost::heap::fibonacci_heap< DistanceIndexPair,boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
//            DistanceIndexHeap;
//        typedef DistanceIndexHeap::handle_type
//            DistanceIndexHeapHandle;
//
//        Point startPoint = V.at(start);
//        Point goalPoint = V.at(goal);
//        EuclideanDistance h = { V.at(goal) }; // initialize heuristic functor
//
//        DistanceIndexHeap open;
//        unordered_map<size_t,DistanceIndexHeapHandle> openHeapHandle(n);
//        openHeapHandle[start] = open.emplace( h( startPoint ), start );
//
//        //unordered_set<size_t> closed(n);
//        //vector<size_t> parents(n);
//
//        vector<double>& g = ShortestKnownPaths.at(start);
//
//        vector<double> f( n, INF );
//        f[start] = h( startPoint );
//
//        DistanceIndexPair current = open.top(); // initialize current vertex to start
//        size_t u_index = current.second;
//        Point currentPoint = startPoint;
//        Point neighborPoint;
//        cout<<"\n    A* start:"<<startPoint;
//        cout<<",";
//        cout<<" goal:"<<V.at(goal);
//        cout<<",";
//
//        do {
//            current = open.top();
//            open.pop();
//
//            u_index = current.second;
//            currentPoint = V.at(u_index);
//            cout<<"\n      current:"<<currentPoint;
//            cout<<",";
//            if( u_index == goal ) { // found the optimal path to the goal, quit
//                open.clear();
//                cout<<"found optimal path to goal\n";
//                break;
//            }
////            cout<<" no goal, ";
//            double t_new = 0;
//            // loop through neighbors of current
//            for( size_t neighbor : E.at(u_index) ) {
//                neighborPoint = V.at(neighbor);
//                cout<<"\n        n:"<<neighborPoint;
//                cout<<",";
//                double newScore = g.at(u_index)
//                    + EuclideanDistances.at(u_index).at(neighbor);
//                cout<<"g_old:"<<g.at(neighbor);
//                cout<<",";
//                cout<<"g_new:"<<newScore;
//                cout<<",";
//                if( newScore < g.at( neighbor ) ) {
//                    //parents[neighbor] = u_index;
//                    g[neighbor] = newScore;
//                    f[neighbor] = g.at(neighbor) + h(neighborPoint);
//                    DistanceIndexPair q = make_pair( f.at(neighbor), neighbor );
//
//                    // calculate the new path's t
//                    t_new = newScore / EuclideanDistances.at(start).at(u_index);
//                    cout<<"t_new:";
//                    cout<<t_new;
//                    cout<<",";
//                    // update t_upper in t-Heap
//                    UpperBoundHeapHandle tHandle = upperBoundHandles.at(start).at(u_index);
//                    auto tValue = make_pair( t_new, make_pair(start,u_index) );
//                    upperBounds.update( tHandle, tValue );
//
//                    if( contains( openHeapHandle, neighbor ) ) {
//                        DistanceIndexHeapHandle neighborHandle = openHeapHandle.at(neighbor);
//                        open.update( neighborHandle, q );
//                        open.update( neighborHandle );
//                    } else {
//                        openHeapHandle[neighbor] = open.push(q);
//                    }
//                    // if we found the goal and the t is good, quit
//                    // this may not be a good idea, honestly, but it might save a lot of time
//                    if( neighbor == goal && t_new < t_max ) {
//                        open.clear();
//                        cout<<"found a path to goal\n";
//                        break;
//                    }
//                }
//            }
//        } while( !open.empty() );
//
//        // END ASTAR
//
//        // Calculate new t for this pair and set new t_max if necessary
//        t_currentPair = ShortestKnownPaths.at(start).at(goal) / EuclideanDistances.at(start).at(goal);
//        cout<<"t_current:";
//        cout<<t_currentPair;
//        cout<<",";
//
//        t_max = CGAL::max( t_max, t_currentPair );
//        cout<<"t_max:";
//        cout<<t_max;
//        cout<<",";
//
//        cout<<"\n";
//    }
//
//    return t_max;
//}

class Timer {
  public:
    explicit Timer( std::string delimiter = "," ) : m_delimiter(delimiter) {
        m_startTime = std::chrono::high_resolution_clock::now();
    }
    ~Timer() {
        stop();
    }
    void stop() {
        auto endTime = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startTime).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTime).time_since_epoch().count();
        auto duration = end - start;

        std::cout << duration << m_delimiter;
    }
  private:
    std::chrono::time_point< std::chrono::high_resolution_clock > m_startTime;
    std::string m_delimiter;
};

} // namespace gsnunf


#endif // GSNUNF_METRICS_H


