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
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
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

namespace unf_spanners {

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

    template<typename VertexIterator, typename EdgeIterator, typename OutputIterator>
    number_t SFWorstPath(VertexIterator pointsBegin,
                         VertexIterator pointsEnd,
                         EdgeIterator edgesBegin,
                         EdgeIterator edgesEnd,
                         std::optional<OutputIterator> out = std::nullopt) {
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
        //vector<double> T( n, INF );
        double t_max = 0.0;

        vector<index_t> MaxParents;
        index_t i_max = 0, j_max = 1;

        // calculate euclidean getDistance between all pairs
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
                *(*out)++ = make_pair(walk, MaxParents.at(walk));
                walk = MaxParents.at(walk);
            } while (walk != i_max);
        }
        // Find the big mac daddy t aka big money
        return t_max;
    }

    class DelaunayTriangulationSFH {
        typedef CGAL::Triangulation_vertex_base_with_info_2<index_t, K> Vb;
        typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
        typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
        class PointSFH : public K::Point_2 {
        public:
            index_t id = numeric_limits<index_t>::max(); // dummy value for id
            string color = "white";

            PointSFH() = default;

            PointSFH(number_t X, number_t Y) : K::Point_2(X, Y) { }
            PointSFH(number_t X, number_t Y, index_t k) : K::Point_2(X, Y) { id = k; }

            friend ostream &operator<<(ostream &strm, const PointSFH &p) {
                return strm << p.id << ": (" << p.x() << ", " << p.y() << ")";
            }
        };

        Delaunay T;
    public:
        DelaunayTriangulationSFH(const vector<Point> &P, vector<Edge> &E) {
            insertPoints(P);
            getEdges(E);
        }

        DelaunayTriangulationSFH(const vector<Point> &P) {
            insertPoints(P);
        }

        void getEdges(vector<Edge> &E) {
            for(auto it = T.finite_edges_begin(); it != T.finite_edges_end(); ++it) {
                Delaunay::Edge e=*it;
                index_t i1 = e.first->vertex( (e.second+1)%3 )->info();
                index_t i2 = e.first->vertex( (e.second+2)%3 )->info();
                E.emplace_back(i1,i2);
            }
        }

        void insertPoints(const vector<Point> &P) {
            vector<pair<PointSFH, index_t> > points;
            index_t id = 0;

            for( Point p : P )
                points.emplace_back(make_pair(PointSFH(p.x(),p.y()),id++));

            T.insert(points.begin(), points.end());
        }

        index_t findClosestPoint(const Point &p) {
            auto handle = T.nearest_vertex(p);
            return handle->info();
        }
    };

    struct priorityComparator{
        bool operator()(const pair<number_t,index_t> &p1, const pair<number_t,index_t> &p2) const{
            return p1.first > p2.first;
        }
    };

    typedef boost::heap::fibonacci_heap<pair<number_t,index_t>, boost::heap::compare<priorityComparator>> FibonacciHeap;

    number_t aStar(const vector<Point> &P, const vector<unordered_set<index_t>> &G, const index_t startVertex, const index_t goalVertex){

        FibonacciHeap openSet;
        vector<index_t> cameFrom(P.size(),-1);
        vector<bool> isInOpenSet(P.size(),false);
        vector<FibonacciHeap::handle_type> handleOfVertex(P.size());
        vector<number_t>  g(P.size(), INF);

        handleOfVertex[startVertex] = openSet.push(make_pair(getDistance(P[startVertex],P[goalVertex]),startVertex));
        isInOpenSet[startVertex] = true;
        g[startVertex] = 0;

        //unsigned count = 0;
        //bool flag = false;

        while( !openSet.empty() ){
            // count++;

//        if( count == 400 ) {
//            flag = true;
//            break;
//        }

            pair<number_t,index_t> current = openSet.top();
            openSet.pop();
            isInOpenSet[current.second] = false;

            if( current.second == goalVertex ) {
                number_t pathLength = 0;
                index_t currentVertex = current.second;
                while ( currentVertex !=  startVertex ) {
                    pathLength += getDistance( P[currentVertex], P[cameFrom[currentVertex]]);
                    currentVertex = cameFrom[currentVertex];
                }
                return pathLength;
            }

            for( index_t neighbor : G[current.second] ) {

                number_t tentativeGscore = g[current.second] + getDistance(P[current.second],P[neighbor]);
                if( tentativeGscore < g[neighbor]){
                    cameFrom[neighbor] = current.second;
                    g[neighbor] = tentativeGscore;
                    number_t fOfNeighbor = g[neighbor] + getDistance(P[neighbor], P[goalVertex]);

                    if (!isInOpenSet[neighbor]) {
                        handleOfVertex[neighbor] = openSet.push(make_pair(fOfNeighbor, neighbor));
                        isInOpenSet[neighbor] = true;
                    }
                    else
                        openSet.decrease(handleOfVertex[neighbor],make_pair(fOfNeighbor, neighbor));
                }
            }
        }

        return INF;
    }

    template<typename VertexIterator, typename EdgeIterator>
    number_t StretchFactorUsingHeuristic(VertexIterator pointsBegin,
                                                VertexIterator pointsEnd,
                                                EdgeIterator edgesBegin,
                                                EdgeIterator edgesEnd,
                                                const size_t numberOfThreads = 4) {
        const vector<Point> P(pointsBegin, pointsEnd);
        const vector<Edge> E(edgesBegin, edgesEnd);
        vector<unordered_set<index_t>> G(P.size()),
            DelG(P.size());

        for (Edge e : E) {
            G[e.first].insert(e.second);
            G[e.second].insert(e.first);
        }

        vector<Edge> edgesOfDT;
        DelaunayTriangulationSFH DT(P, edgesOfDT);
        for(Edge e : edgesOfDT ) {
            DelG[e.first].insert(e.second);
            DelG[e.second].insert(e.first);
        }

        vector<number_t> stretchFactorOfG(numberOfThreads,0);
        vector<pair<index_t, index_t>> worstPairOfG(numberOfThreads);
        vector<unordered_map<index_t,number_t>> tracker(P.size());

#pragma omp parallel for num_threads(numberOfThreads)
        for( index_t u = 0; u < P.size(); u++ ) {

            number_t largestSfFromU = 0, currentStretchFactor = 0;
            pair<index_t, index_t> worstPairSoFar, worstPairCurrent;
            unordered_set<index_t> localPoints;

            vector<bool> isConsidered(P.size(),false);
            isConsidered[u] = true;
            localPoints.insert(u);

            while( true ) {

                unordered_set<index_t> underConsideration;
                for (auto i : localPoints) {
                    for (auto j : DelG[i]) {
                        if (!isConsidered[j]) {
                            underConsideration.insert(j);
                            isConsidered[j] = true;
                        }
                    }
                }

                for (index_t v : underConsideration) {

                    unordered_map<index_t,number_t>::iterator result;

#pragma omp critical
                    result = tracker[std::min(u,v)].find(std::max(u,v));

                    if (result == tracker[u].end()) {
                        number_t stretchFactorUV = aStar(P, G, u, v) / getDistance(P[u], P[v]);

                        if (stretchFactorUV > currentStretchFactor) {
                            worstPairCurrent.first  = u;
                            worstPairCurrent.second = v;
                        }

                        currentStretchFactor = std::max(currentStretchFactor, stretchFactorUV);

#pragma omp critical
                        tracker[std::min(u,v)].insert(make_pair(std::max(u,v), stretchFactorUV));
                    }
                    else {
                        if( (*result).second > currentStretchFactor ){
                            worstPairCurrent.first = u;
                            worstPairCurrent.second = v;
                            currentStretchFactor = (*result).second;
                        }
                    }
                }

                if (currentStretchFactor > largestSfFromU) {
                    largestSfFromU = currentStretchFactor;
                    worstPairSoFar.first  = worstPairCurrent.first;
                    worstPairSoFar.second = worstPairCurrent.second;
                }
                else
                    break;

                localPoints.clear();
                for( index_t i : underConsideration )
                    localPoints.insert(i);
            }

            if (largestSfFromU > stretchFactorOfG[omp_get_thread_num()]) {
                worstPairOfG[omp_get_thread_num()].first  = worstPairSoFar.first;
                worstPairOfG[omp_get_thread_num()].second = worstPairSoFar.second;
                stretchFactorOfG[omp_get_thread_num()]    = largestSfFromU;
            }
        }

        index_t worstPairU = worstPairOfG[0].first, worstPairV = worstPairOfG[0].second;
        number_t finalSf = stretchFactorOfG[0];
        for(index_t i = 1; i < numberOfThreads; i++)
            if( stretchFactorOfG[i] >  finalSf ) {
                finalSf = stretchFactorOfG[i];
                worstPairU = worstPairOfG[i].first;
                worstPairV = worstPairOfG[i].second;
            }

        //cout << "Heuristic worst pair: " <<  worstPairU << ", " << worstPairV << endl;
        return *std::max_element(stretchFactorOfG.begin(),stretchFactorOfG.end());
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
        explicit Timer(std::string delimiter = ",") : running(true), m_delimiter(std::move(delimiter)) {
            m_startTime = std::chrono::high_resolution_clock::now();
        }

        ~Timer() {
            if( running ) {
                std::cout << stop() << m_delimiter;
            }
        }

        size_t stop() {
            auto endTime = std::chrono::high_resolution_clock::now();
            auto start = std::chrono::time_point_cast<std::chrono::microseconds>(
                    m_startTime).time_since_epoch().count();
            auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTime).time_since_epoch().count();
            auto duration = end - start;
            running = false;

            return duration;
        }

    private:
        bool running;
        std::chrono::time_point<std::chrono::high_resolution_clock> m_startTime;
        std::string m_delimiter;
    };

} // namespace unf_spanners


#endif // GSNUNF_METRICS_H


