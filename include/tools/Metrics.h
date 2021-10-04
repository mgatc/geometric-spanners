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

        do {
            current = open.begin();

            u_index = current->second;
            currentPoint = V[u_index];

            // loop through neighbors of current
            for (index_t neighbor : G.at(u_index)) {
                neighborPoint = V[neighbor];
                newScore = ShortestPaths[u_index]
                           + getDistance(currentPoint, neighborPoint);
                if (newScore < ShortestPaths[neighbor]) {
                    Parents[neighbor] = u_index;
                    ShortestPaths[neighbor] = newScore;

                    if (contains(handleToHeap, neighbor)) {
                        open.erase(handleToHeap[neighbor]);
                    }
                    handleToHeap[neighbor] = open.emplace(ShortestPaths[neighbor], neighbor).first;
                }
            }
            open.erase(current);
        } while (!open.empty());
    }

    template<typename VertexIterator, typename EdgeIterator>
    number_t StretchFactorDijkstraReduction(VertexIterator pointsBegin,
                                            VertexIterator pointsEnd,
                                            EdgeIterator edgesBegin,
                                            EdgeIterator edgesEnd,
                                            const size_t numberOfThreads = 4) {
        typedef typename VertexIterator::value_type Point_2;

        const vector<Point_2> V(pointsBegin, pointsEnd); // container for vertices
        const size_t n = V.size();

        vector<unordered_set<size_t> > G(n, unordered_set<size_t>()); // adjacency list

        // Create list of vertices, map to their indices, and adjacency list
        for (auto eit = edgesBegin; eit != edgesEnd; ++eit) {
            auto p = eit->first,
                 q = eit->second;

            G[p].insert(q);
            G[q].insert(p);
        }
        double t_max = 0.0;

        #pragma omp parallel for reduction(max: t_max) shared(G) default(none) num_threads(numberOfThreads)
        for (size_t i = 0; i < n; ++i) {
            // Euclidean distances
            vector<number_t> D(n);
            for (size_t j = 0; j < n; ++j) {
                D[j] = i == j ? 0 : getDistance(V[i], V[j]);
            }
            // Shortest paths
            vector<number_t> ShortestPaths(n, INF);
            vector<size_t> Parents(n);
            Dijkstra(i, V, G, ShortestPaths, Parents);

            // Divide each shortest path getDistance by the euclidean distance between the vertices.
            for (size_t j = 0; j < n; ++j) {
                ShortestPaths[j] = ( // avoid /0
                        i == j ? 0 : ShortestPaths[j] / D[j]
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
        // Find the big mac daddy stretchFactor aka big money
        return t_max;
    }

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
            add_vertex( *pit, G );
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

        for(auto ed : mst){
            VertexDescriptor u = source(ed, G),
                             v = target(ed, G);
            *out = make_pair(u,v);
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
            if( m_clock.is_running() )
                std::cout << stop() << m_delimiter;
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


