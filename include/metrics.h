#ifndef GSNUNF_METRICS_H
#define GSNUNF_METRICS_H

#include <algorithm> // swap
#include <functional>
#include <limits>
#include <optional>
#include <unordered_map>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h> //for 2D functions

#include <omp.h>

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

inline double d( Point p, Point q ) {
    return CGAL::sqrt( CGAL::squared_distance(p,q) );
}

struct EuclideanDistance {
    Point goal;
    double operator()( Point p ) {
        return d(p,goal);
    }
};

template<typename T>
struct MinHeapCompare {
    bool operator()( const T &n1, const T &n2 ) const {
        return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
    }
};

template< typename VertexContainer, typename VertexMap, typename AdjacencyList >
optional<double> AStar( VertexContainer V, VertexMap vMap, AdjacencyList G_prime, size_t start, size_t goal ) {
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

    vector<double> g( n, inf );
    g[start] = 0;

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
        if( u_index == goal ) return make_optional( g.at(goal) );
//        cout<<" no goal, ";
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
void Djikstra( const size_t i, const VertexContainer& V, const AdjacencyList& G, vector<double>& ShortestPaths ) {
    typedef pair<double,size_t>
        DistanceIndexPair;
    typedef boost::heap::fibonacci_heap< DistanceIndexPair,boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
        Heap;
    typedef Heap::handle_type
        HeapHandle;

    size_t n = V.size();
    size_t inf = numeric_limits<size_t>::max();
    Point startPoint = V.at(i);

    Heap open;
    unordered_map<size_t,HeapHandle> handleToHeap(n);
    handleToHeap[i] = open.emplace( 0, i );

    //unordered_set<size_t> closed(n);
    vector<size_t> parents(n);

    ShortestPaths[i] = 0;

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
        currentPoint = V.at(u_index);
//        cout<<"\n      current:"<<currentPoint;
//        cout<<",";
        // loop through neighbors of current
        for( size_t neighbor : G.at(u_index) ) {
            neighborPoint = V.at(neighbor);
//            cout<<"\n        n:"<<neighborPoint;
//            cout<<",";
            double newScore = ShortestPaths.at(u_index)
                + d( currentPoint, neighborPoint );
//            cout<<"g:"<<newScore;
//            cout<<",";
            if( newScore < ShortestPaths.at( neighbor ) ) {
                parents[neighbor] = u_index;
                ShortestPaths[neighbor] = newScore;
                DistanceIndexPair q = make_pair( ShortestPaths.at(neighbor), neighbor );

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
}

template< typename RandomAccessIterator >
double StretchFactorDjikstra( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
    vector<Point> V; // container for vertices
    unordered_map< Point, size_t, PointHasher > vMap; // map point to index in V
    unordered_map< size_t, unordered_set<size_t> > G; // adjacency list
    size_t index = 0;

    // Create list of vertices, map to their indices, and adjacency list
    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) {
        // If vMap doesn't contain p, put it in V
        Point p = eit->first;
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
        Point q = eit->second;
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
    const double INF = numeric_limits<double>::max();
    vector<vector<double> > ShortestPaths(n, vector<double>(n, INF) );
    vector<vector<double> > D(n, vector<double>(n, INF) ); // Euclidean distances
    vector<double> T( n, INF );

    // calculate euclidean distance between all pairs
    //#pragma omp parallel for
    for( size_t i=0; i<n; ++i ) {
        for( size_t j=0; j<n; ++j ) {
            D.at(i).at(j) =
                i==j ? 0 : d( V.at(i), V.at(j) );
        }
    }

    // linear scan over vertices, perform single-source shortest paths
    //#pragma omp parallel for
    for( size_t i=0; i<n; ++i ) {
        Djikstra( i, V, G, ShortestPaths.at(i) );
    }

    // Divide each shortest path distance by the euclidean distance between the vertices.
    //#pragma omp parallel for
    for( size_t i=0; i<n; ++i ) {
        for( size_t j=0; j<n; ++j ) {
            ShortestPaths.at(i).at(j) = ( // avoid /0
                i==j ? 0 : ShortestPaths.at(i).at(j)/D.at(i).at(j)
            );
        }
        // Find max t and place in T
        T.at(i) = *max_element(
            begin( ShortestPaths.at(i) ),
            end(   ShortestPaths.at(i) )
        );
    }
    // Find the big mac daddy t aka big money
    return *max_element( T.begin(), T.end() );
}

template< typename RandomAccessIterator >
double StretchFactorDjikstraParallel( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
    vector<Point> V; // container for vertices
    unordered_map< Point, size_t, PointHasher > vMap; // map point to index in V
    unordered_map< size_t, unordered_set<size_t> > G; // adjacency list
    size_t index = 0;

    // Create list of vertices, map to their indices, and adjacency list
    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) {
        // If vMap doesn't contain p, put it in V
        Point p = eit->first;
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
        Point q = eit->second;
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
    const double INF = numeric_limits<double>::max();
    vector<vector<double> > ShortestPaths(n, vector<double>(n, INF) );
    vector<vector<double> > D(n, vector<double>(n, INF) ); // Euclidean distances
    vector<double> T( n, INF );

    // calculate euclidean distance between all pairs
    #pragma omp parallel for
    for( size_t i=0; i<n; ++i ) {
        for( size_t j=0; j<n; ++j ) {
            D.at(i).at(j) =
                i==j ? 0 : d( V.at(i), V.at(j) );
        }
    }

    // linear scan over vertices, perform single-source shortest paths
    #pragma omp parallel for
    for( size_t i=0; i<n; ++i ) {
        Djikstra( i, V, G, ShortestPaths.at(i) );
    }

    // Divide each shortest path distance by the euclidean distance between the vertices.
    #pragma omp parallel for
    for( size_t i=0; i<n; ++i ) {
        for( size_t j=0; j<n; ++j ) {
            ShortestPaths.at(i).at(j) = ( // avoid /0
                i==j ? 0 : ShortestPaths.at(i).at(j)/D.at(i).at(j)
            );
        }
        // Find max t and place in T
        T.at(i) = *max_element(
            begin( ShortestPaths.at(i) ),
            end(   ShortestPaths.at(i) )
        );
    }
    // Find the big mac daddy t aka big money
    return *max_element( T.begin(), T.end() );
}

template< typename VertexContainer, typename VertexMap, typename AdjacencyList >
optional<double> ShortestPath( VertexContainer V, VertexMap vMap, AdjacencyList G_prime, size_t start, size_t goal ) {
    return AStar( V, vMap, G_prime, start, goal );
}

// edgelists are range-supporting containers containing pairs of Points
template< typename RandomAccessIterator, typename Triangulation >
double StretchFactor( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd, const Triangulation& superGraph ) {
    using GraphVertex = typename Triangulation::Vertex_handle;
    using GraphCirculator = typename Triangulation::Vertex_circulator;
    //using GPoint = typename Triangulation::Point_2;
    Triangulation G(superGraph);
    // create vector of points from given range and map points to indices
    vector<Point> P;
    P.reserve( G.number_of_vertices() );
    vector<GraphVertex> V;
    V.reserve( P.capacity() );
    unordered_map< GraphVertex, size_t, boost::hash<GraphVertex> > vMap( P.capacity() );
    size_t incidenceListInitialSize = 3;
    unordered_map< size_t, unordered_set<size_t> > G_prime( P.capacity() ); // E = O(n) so this is a decent enough guess
    size_t i = 0;

    for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) { // Add all vertices. If insertion successful, put vertex in V and increment i.
        // first vertex in pair
        Point p = eit->first;
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
        Point q = eit->second;
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
            Point p_u = P.at( u_index ),
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
            Vertex_circulator v_n = G.incident_vertices( V.at(v.first) ),
                done(v_n);
            do if( !G.is_infinite(v_n) ) {
                    size_t i = vMap.at(v_n);
                    //cout<<i<<" "<<P.at(i)<<" "<<(status.at(i) ? "unknown" : "known")<<"... ";
                    if( status.at(i) ) {
                        level.emplace(i,lvl+1);
                        status[i] = false;
                    }
            } while( --v_n != done );
            // END BFS ADMIN
//            cout<<"\n";
        } while( !level.empty() );
//        cout<<"\n";
    }

    return t_max;
}

} // namespace gsnunf


#endif // GSNUNF_STRETCHFACTOR_H


