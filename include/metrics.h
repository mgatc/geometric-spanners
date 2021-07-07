#ifndef GSNUNF_METRICS_H
#define GSNUNF_METRICS_H

#include <algorithm> // swap
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
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

namespace unf_planespanners {

using namespace std;


    template<typename RandomAccessIterator>
    size_t degree(RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd) {
        typedef typename RandomAccessIterator::value_type EdgeType;
        typedef typename EdgeType::first_type VertexType;

        std::vector<EdgeType> edges(edgesBegin, edgesEnd);
        std::unordered_map<VertexType, unordered_set<VertexType>> adj;
        // for each edge
        for (auto e : edges) {
            auto first = adj.begin();
            tie(first, ignore) = adj.emplace(e.first, unordered_set<VertexType>());
            (*first).second.insert(e.second);

            auto second = adj.begin();
            tie(second, ignore) = adj.emplace(e.second, unordered_set<VertexType>());
            (*second).second.insert(e.first);
        }
        auto max_el = max_element(adj.begin(), adj.end(), [&](const auto &lhs, const auto &rhs) {
            return lhs.second.size() < rhs.second.size();
        });

        return max_el->second.size();
    }

    template<typename Triangulation>
    size_t degree(const Triangulation &T) {
        typedef typename Triangulation::Point
                Point_2;

        // fill a vector with edges so we can call the range-based degree function
        std::vector<pair<Point_2, Point_2>> edges;
        edges.reserve(T.number_of_vertices());

        for (auto e = T.finite_edges_begin(); e != T.finite_edges_end(); ++e) {
            auto p = make_pair(
                    e->first->vertex((e->second + 1) % 3)->point(),
                    e->first->vertex((e->second + 2) % 3)->point()
            );
            // Add both in and out edges
            forBoth(p, [&](Point a, Point b) {
                edges.emplace_back(a, b);
            });
        }
        return degree(edges.begin(), edges.end());
    }

    template<typename RandomAccessIterator>
    number_t weight(RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd) {
        number_t w = 0;
        for (auto e = edgesBegin; e != edgesEnd; ++e) {
            w += getDistance(e->first, e->second);
        }
        return w;
    }

    template<typename Triangulation>
    number_t weight(const Triangulation &T) {
        number_t w = 0;
        for (auto e = T.finite_edges_begin(); e != T.finite_edges_end(); ++e) {
            auto p = make_pair(
                    e->first->vertex((e->second + 1) % 3)->point(),
                    e->first->vertex((e->second + 2) % 3)->point()
            );
            w += getDistance(p.first, p.second);
        }
        return w;
    }

    template<typename VertexContainer, typename AdjacencyList>
    optional<number_t> AStar(VertexContainer V, AdjacencyList G_prime, index_t start, index_t goal) {
        typedef pair<number_t, index_t>
                DistanceIndexPair;
        typedef boost::heap::fibonacci_heap<DistanceIndexPair, boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
                Heap;
        typedef Heap::handle_type
                HeapHandle;

        index_t n = V.size();
        auto startPoint = V.at(start)->point();
        auto goalPoint = V.at(goal)->point();
        EuclideanDistanceToPoint h = {V.at(goal)->point()}; // initialize heuristic functor

        Heap open;
        unordered_map<index_t, HeapHandle> handleToHeap(n);
        handleToHeap[start] = open.emplace(h(startPoint), start);

        //unordered_set<size_t> closed(n);
        vector<index_t> parents(n);

        vector<number_t> g(n, INF);
        g[start] = 0;

        vector<number_t> f(n, INF);
        f[start] = h(startPoint);

        DistanceIndexPair current = open.top(); // initialize current vertex to start
        index_t u_index = current.second;
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
            if (u_index == goal) return make_optional(g.at(goal));
//        cout<<" no goal, ";
            // loop through neighbors of current
            for (size_t neighbor : G_prime.at(u_index)) {
                neighborPoint = V.at(neighbor)->point();
//            cout<<"\n        n:"<<neighborPoint;
//            cout<<",";
                number_t newScore = g.at(u_index)
                                    + d(currentPoint, neighborPoint);
//            cout<<"g:"<<newScore;
//            cout<<",";
                if (newScore < g.at(neighbor)) {
                    parents[neighbor] = u_index;
                    g[neighbor] = newScore;
                    f[neighbor] = g.at(neighbor) + h(neighborPoint);
                    DistanceIndexPair q = make_pair(f.at(neighbor), neighbor);

                    if (contains(handleToHeap, neighbor)) {
                        HeapHandle neighborHandle = handleToHeap.at(neighbor);
                        open.update(neighborHandle, q);
                        open.update(neighborHandle);
                    } else {
                        handleToHeap[neighbor] = open.push(q);
                    }
                }
            }
        } while (!open.empty());

        return nullopt;
    }


    template<typename VertexContainer, typename AdjacencyList>
    void Dijkstra(const index_t i,
                  const VertexContainer &V,
                  const AdjacencyList &G,
                  vector<number_t> &ShortestPaths,
                  vector<index_t> &Parents) {

//    typedef pair<double,size_t>
//        DistanceIndexPair;
        //typedef boost::heap::fibonacci_heap< DistanceIndexPair,boost::heap::compare<MinHeapCompare<DistanceIndexPair>>>
        typedef std::map<number_t, index_t>
                Heap;
        typedef Heap::iterator
                HeapHandle;

        const index_t n = V.size();
        auto startPoint = V.at(i);

        Heap open;
        unordered_map<index_t, HeapHandle> handleToHeap(n);
        handleToHeap[i] = open.emplace(0, i).first;

        ShortestPaths[i] = 0;

        auto current = open.begin(); // initialize current vertex to start
        index_t u_index = current->second;
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
            for (index_t neighbor : G.at(u_index)) {
                neighborPoint = V[neighbor];
//            cout<<"\n        n:"<<neighborPoint;
//            cout<<",";
                newScore = ShortestPaths[u_index]
                           + getDistance(currentPoint, neighborPoint);
//            cout<<"g:"<<newScore;
//            cout<<",";
                if (newScore < ShortestPaths[neighbor]) {
                    Parents[neighbor] = u_index;
                    ShortestPaths[neighbor] = newScore;
                    //DistanceIndexPair q = make_pair( ShortestPaths.at(neighbor), neighbor );
                    if (contains(handleToHeap, neighbor)) {
                        open.erase(handleToHeap[neighbor]);
                    }
                    handleToHeap[neighbor] = open.emplace(ShortestPaths[neighbor], neighbor).first;
                }
            }
            //closed.insert(current.second);
            open.erase(current);
        } while (!open.empty());
    }

    template<typename VertexIterator, typename EdgeIterator>
    number_t StretchFactorDijkstraReduction(VertexIterator pointsBegin,
                                            VertexIterator pointsEnd,
                                            EdgeIterator edgesBegin,
                                            EdgeIterator edgesEnd) {
        typedef typename VertexIterator::value_type Point_2;

        const vector<Point_2> V(pointsBegin, pointsEnd); // container for vertices
        const size_t n = V.size();
//    unordered_map< size_t, size_t > vMap; // map point to index in V
        vector<unordered_set<size_t> > G(n, unordered_set<size_t>()); // adjacency list
        //size_t index = 0;

        // Create list of vertices, map to their indices, and adjacency list
        for (auto eit = edgesBegin; eit != edgesEnd; ++eit) {
            auto p = eit->first,
                    q = eit->second;

            G[p].insert(q);
            G[q].insert(p);
        }
        //vector<double> T( n, INF );
        double t_max = 0.0;

        // calculate euclidean getDistance between all pairs
#pragma omp parallel for reduction( max: t_max ) default( shared )
        for (size_t i = 0; i < n; ++i) {
            // Euclidean distances
            vector<number_t> D(n, INF);
            for (size_t j = 0; j < n; ++j) {
                D.at(j) =
                        i == j ? 0 : getDistance(V.at(i), V.at(j));
            }
            // Shortest paths
            vector<number_t> ShortestPaths(n, INF);
            vector<size_t> Parents(n);
            Dijkstra(i, V, G, ShortestPaths, Parents);

            // Divide each shortest path getDistance by the euclidean distance between the vertices.
            for (size_t j = 0; j < n; ++j) {
                ShortestPaths.at(j) = ( // avoid /0
                        i == j ? 0 : ShortestPaths.at(j) / D.at(j)
                );
            }
            // Find max_t
            auto t_local = max_element(
                    begin(ShortestPaths),
                    end(ShortestPaths)
            );
            if (*t_local > t_max) {
                t_max = *t_local;
            }
        }
        // Find the big mac daddy t aka big money
        return t_max;
    }

    template<typename RandomAccessIterator, typename OutputIterator>
    number_t SFWorstPath(RandomAccessIterator edgesBegin,
                         RandomAccessIterator edgesEnd,
                         std::optional<OutputIterator> out = std::nullopt) {
        typedef typename RandomAccessIterator::value_type EdgeType;
        typedef typename EdgeType::first_type VertexType;

        vector<VertexType> V; // container for vertices
        unordered_map<VertexType, index_t, PointHasher<VertexType> > vMap; // map point to index in V
        unordered_map<index_t, unordered_set<index_t> > G; // adjacency list
        index_t index = 0;

        // Create list of vertices, map to their indices, and adjacency list
        for (auto eit = edgesBegin; eit != edgesEnd; ++eit) {
            // If vMap doesn't contain p, put it in V
            VertexType p = eit->first;
            index_t i_p = index;
            bool inserted = false;
            auto vMapIt = vMap.begin();
            tie(vMapIt, inserted) = vMap.emplace(p, i_p); // map for reverse lookup
            if (inserted) {
                V.push_back(p);
                ++index;
            }
            i_p = vMapIt->second;

            // If vMap doesn't contain q, put it in V
            VertexType q = eit->second;
            index_t i_q = index;
            tie(vMapIt, inserted) = vMap.emplace(q, i_q); // map for reverse lookup
            if (inserted) {
                V.push_back(q);
                ++index;
            }
            i_q = vMapIt->second;

            G[i_p].insert(i_q); // add edge to adjacency list
        }
        index_t n = V.size();
        //vector<double> T( n, INF );
        number_t t_max = 0.0;
        vector<index_t> MaxParents;
        index_t i_max = 0, j_max = 1;

        // calculate euclidean getDistance between all pairs
        //#pragma omp parallel for reduction( max: t_max )
        for (index_t i = 0; i < n; ++i) {
            // Euclidean distances
            vector<number_t> D(n, INF);
            for (index_t j = 0; j < n; ++j) {
                D.at(j) =
                        i == j ? 0 : getDistance(V.at(i), V.at(j));
            }
            // Shortest paths
            vector<number_t> ShortestPaths(n, INF);
            vector<index_t> Parents(n);
            Dijkstra(i, V, G, ShortestPaths, Parents);

            // Divide each shortest path distance by the euclidean getDistance between the vertices.
            for (index_t j = 0; j < n; ++j) {
                ShortestPaths.at(j) = ( // avoid /0
                        i == j ? 0 : ShortestPaths.at(j) / D.at(j)
                );
            }
            // Find max_t
            auto t_local = max_element(
                    begin(ShortestPaths),
                    end(ShortestPaths)
            );
            if (*t_local > t_max) {
                t_max = *t_local;
                // remove the following for parallel reduction function
                if (out) {
                    std::swap(Parents, MaxParents);
                    i_max = i;
                    j_max = t_local - ShortestPaths.begin();
                }
            }
        }
        if (out) {
            size_t walk = j_max;
            do {
                *(*out)++ = make_pair(V.at(walk), V.at(MaxParents.at(walk)));
                walk = MaxParents.at(walk);
            } while (walk != i_max);
        }
        // Find the big mac daddy t aka big money
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

    class Timer {
    public:
        explicit Timer(std::string delimiter = ",") : m_delimiter(std::move(delimiter)) {
            m_startTime = std::chrono::high_resolution_clock::now();
        }

        ~Timer() {
            stop();
        }

        void stop() {
            auto endTime = std::chrono::high_resolution_clock::now();
            auto start = std::chrono::time_point_cast<std::chrono::microseconds>(
                    m_startTime).time_since_epoch().count();
            auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTime).time_since_epoch().count();
            auto duration = end - start;

            std::cout << duration << m_delimiter;
        }

    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> m_startTime;
        std::string m_delimiter;
    };

} // namespace unf_planespanners


#endif // GSNUNF_METRICS_H


