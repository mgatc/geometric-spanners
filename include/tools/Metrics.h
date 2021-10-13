#ifndef UNF_SPANNERS_METRICS_H
#define UNF_SPANNERS_METRICS_H

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
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/Real_timer.h>
#include <CGAL/squared_distance_2.h> //for 2D functions

#include <omp.h>

#include "tools/DelaunayL2.h"
#include "tools/FloydWarshall.h"
#include "tools/StretchFactorExact.h"
#include "tools/StretchFactorExperimental.h"
#include "tools/Utilities.h"

namespace spanners {

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
//        cout<<"Largest degree vertex="<<max_el->first<<endl;

        return max_el->second.size();
    }
    template<typename RandomAccessIterator>
    number_t degreeAvg(RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd) {
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
        auto avg = std::accumulate(adj.begin(), adj.end(), 0.0, [&](const number_t &sum, const auto &current) {
            return sum + current.second.size();
        }) / number_t(adj.size());

        return avg;
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

    template< class VertexIterator, class EdgeIterator>
    number_t weight( VertexIterator pointsBegin,
                     VertexIterator pointsEnd,
                     EdgeIterator edgesBegin,
                     EdgeIterator edgesEnd ) {
        vector<Point> P(pointsBegin,pointsEnd);

        number_t w = 0;
        index_t p,q;
        for (auto e = edgesBegin; e != edgesEnd; ++e) {
            tie(p,q) = *e;
            w += getDistance(P[p],P[q]);
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

                // calculate the new path's stretchFactor
                t_new = newScore / EuclideanDistances.at(i).at(u_index);
                // update t_upper in stretchFactor-Heap
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
    template< class VertexIterator, class EdgeIterator, class EdgeOutputIterator >
    void getMST( VertexIterator pointsBegin,
                     VertexIterator pointsEnd,
                     EdgeIterator edgesBegin,
                     EdgeIterator edgesEnd,
                     EdgeOutputIterator out)
    {

        using namespace boost;
        typedef adjacency_list<vecS, vecS, undirectedS,
            Point,
            property<edge_weight_t,number_t>
            > Graph;

        Graph G;

        for( auto pit=pointsBegin; pit!=pointsEnd; ++pit ) {
            auto v = add_vertex( *pit, G );
        }

        index_t p, q;
        for( auto eit=edgesBegin; eit!=edgesEnd; ++eit ) {
            tie( p, q ) = *eit;
            number_t wt = getDistance( G[p], G[q] );
            auto e = add_edge( p, q, wt, G );
        }

        typedef Graph::vertex_descriptor VertexDescriptor;
        typedef Graph::edge_descriptor EdgeDescriptor;

        std::list<EdgeDescriptor> mst;
        boost::kruskal_minimum_spanning_tree(G,std::back_inserter(mst));

        for(auto it = mst.begin(); it != mst.end(); ++it){
            EdgeDescriptor ed = *it;
            VertexDescriptor p = source(ed, G),
                             q = target(ed, G);
            *out = make_pair(p,q);
        }
    }
    template< class VertexIterator, class EdgeIterator>
    number_t getLightness( VertexIterator pointsBegin,
                           VertexIterator pointsEnd,
                           EdgeIterator edgesBegin,
                           EdgeIterator edgesEnd ) {
        vector<Point> P(pointsBegin,pointsEnd);
        vector<Edge> E(edgesBegin,edgesEnd);
        list<Edge> MST;
        getMST( P.begin(), P.end(), E.begin(), E.end(), back_inserter(MST) );
        number_t weightOfMST = weight(P.begin(), P.end(), MST.begin(), MST.end() ),
                 weightOfG   = weight(P.begin(), P.end(), E.begin(), E.end() ),
                 lightness = weightOfG / weightOfMST;
        return lightness;
    }

    class Timer {
    public:
        explicit Timer(std::string delimiter = ",") : m_delimiter(std::move(delimiter)) {
            m_clock.start();
        }

        ~Timer() {
            if( m_clock.is_running() ) {
                std::cout << stop() << m_delimiter;
            }
        }

        double stop() {
            m_clock.stop();
            return m_clock.time();
        }

    private:
        std::string m_delimiter;
        CGAL::Real_timer m_clock;
    };

} // namespace spanners


#endif // UNF_SPANNERS_METRICS_H


