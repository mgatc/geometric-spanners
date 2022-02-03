#ifndef STRETCHFACTOR_STRETCHFACTOREXACT_H
#define STRETCHFACTOR_STRETCHFACTOREXACT_H

#include <algorithm>
#include <map>
#include <unordered_set>
#include <vector>

#include <omp.h> // <-- UNCOMMENT FOR PARALLEL, ALSO LINE 98

#include "Utilities.h"

namespace spanners {

    template<typename VertexContainer, typename AdjacencyList>
    void Dijkstra(const index_t i,
                  const VertexContainer &V,
                  const AdjacencyList &G,
                  std::vector<number_t> &ShortestPaths,
                  std::vector<index_t> &Parents) {

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
        std::unordered_map<index_t, HeapHandle> handleToHeap(n);
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

        const std::vector<Point_2> V(pointsBegin, pointsEnd); // container for vertices
        const index_t n = V.size();
//    unordered_map< size_t, size_t > vMap; // map point to index in V
        std::vector<std::unordered_set<index_t> > G(n, std::unordered_set<index_t>()); // adjacency list
        //size_t index = 0;

        // Create list of vertices, map to their indices, and adjacency list
        for (auto eit = edgesBegin; eit != edgesEnd; ++eit) {
            auto p = eit->first,
                    q = eit->second;

            G[p].insert(q);
            G[q].insert(p);
        }
        //vector<double> T( n, INF );
        number_t t_max = 0.0;

        // calculate euclidean getDistance between all pairs
#pragma omp parallel for reduction( max: t_max ) num_threads(numberOfThreads) default( shared ) // <-- UNCOMMENT FOR PARALLEL, ALSO LINE 9
        for (index_t i = 0; i < n; ++i) {
            // Euclidean distances
            std::vector<number_t> D(n, INF);
            for (index_t j = 0; j < n; ++j) {
                D.at(j) =
                        i == j ? 0 : getDistance(V.at(i), V.at(j));
            }
            // Shortest paths
            std::vector<number_t> ShortestPaths(n, INF);
            std::vector<index_t> Parents(n);
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
        // Find the big mac daddy stretchFactor aka big money
        return t_max;
    }

    template<typename VertexIterator, typename EdgeIterator, typename OutputIterator>
    number_t SFWorstPath(VertexIterator pointsBegin,
                         VertexIterator pointsEnd,
                         EdgeIterator edgesBegin,
                         EdgeIterator edgesEnd,
                         std::optional<OutputIterator> out = std::nullopt) {
        typedef typename VertexIterator::value_type Point_2;

        using std::unordered_set;
        using std::vector;

        vector<Point_2> V(pointsBegin, pointsEnd); // container for vertices
        const index_t n = V.size();

        vector<unordered_set<index_t> > G(n, unordered_set<index_t>()); // adjacency list

        // Create list of vertices, map to their indices, and adjacency list
        for (auto eit = edgesBegin; eit != edgesEnd; ++eit) {
            auto p = eit->first,
                    q = eit->second;

            G[p].insert(q);
            G[q].insert(p);
        }
        //vector<double> T( n, INF );
        number_t t_max = 0.0;

        vector<index_t> MaxParents;
        index_t i_max = 0, j_max = 1;

        // calculate euclidean getDistance between all pairs
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

            // Divide each shortest path getDistance by the euclidean distance between the vertices.
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
        // Find the big mac daddy stretchFactor aka big money
        return t_max;
    }

} // spanners

#endif //STRETCHFACTOR_STRETCHFACTOREXACT_H
