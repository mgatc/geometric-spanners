//Needs optimizing currently testing.
#ifndef GSNUNF_KPT2017_H
#define GSNUNF_KPT2017_H

//Base libraries.
#include <cmath>         // ceil, floor, isinf
#include <limits>
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

//CGAL library
#include <CGAL/algorithm.h>

//Project library
#include "GeometricSpannerPrinter.h"
#include "metrics.h"
#include "TDDelaunay.h"
#include "utilities.h"


namespace planespanners {

    using namespace std;

    namespace kpt2017 {

        typedef HalfThetaTriangulation<K> TDDelaunayGraph;
        typedef TDDelaunayGraph::VertexDescriptor VertexDescriptor;

        enum Color {
            Blue, White
        };

        //Compute max of getCone(p,q) and (getCone(q,p)+3)%6, is used to make sure cones are calculated correctly.
        inline Color getColor(const index_t p, const index_t q, const vector<Point> &H) {
            //cout<<"Getting color of "<<p<<"-"<<q<<endl;
            return getCone(p, q, H) % 3 == 1 ? Blue : White;
        }
        inline Color getColor(const Edge &e, const vector<Point> &H) {
            return getColor(e.first, e.second, H);
        }

        enum ConePolarity {
            Positive = 0, // even cones
            Negative = 1  // odd cones
        };
        inline ConePolarity getConePolarity(const index_t p, const index_t q, const vector<Point> &H) {
            return ConePolarity(getCone(p, q, H) % 2);
        }

        //Finds the bisector length of a given edge.
        number_t bisectorLength(const Edge &e, const vector<Point> &H) {
            cone_t cone = getCone(e.first, e.second, H);
            //assert(cone < 6);
            //assert(e.first < H.size());

            number_t xCord = H[e.first].x() - orthBisectorSlopes[cone];
            number_t yCord = H[e.first].y() + 1;

            Point bisectorPoint(xCord, yCord);
            K::Line_2 bisectorLine(H[e.first], bisectorPoint);
            Point intersectionPoint = bisectorLine.projection(H[e.second]);
            number_t bisectorLen = getDistance(H[e.first], intersectionPoint);

            return bisectorLen;
        }

        template<class AnchorListMap, class Triangulation, class PointContainer>
        void findAnchors(AnchorListMap &anchors,
                         Triangulation &D,
                         const PointContainer &P) {
            anchors =
                    {
                            {Blue,  {}},
                            {White, {}}
                    };

            // find and classify anchors
            for (auto vit = D.finite_vertices_begin();
                 vit != D.finite_vertices_end(); ++vit) {
                auto w = *vit;
                //cout<<w<<") "<<D.point(w)<<endl;

                // negative cones
                VertexDescriptor localAnchors[3];
                number_t localMinimumBisectors[3] = {INF, INF, INF};

                // get the edges for which w is the target
                for (auto eit = D.negative_cone_edges_begin(w);
                     eit != D.negative_cone_edges_end(w); ++eit) {
                    auto e = *eit;
                    auto v = D.source(e);
//            bool same_y = CGAL::compare_y(P[w], P[v]) == CGAL::EQUAL;
                    cone_t cone = getCone(w, v, P);
                    cone_t flattened_cone = cone / 2;
                    auto e_bisector_length = bisectorLength(make_pair(w, v), P);
                    if (e_bisector_length < localMinimumBisectors[flattened_cone]) {
                        localAnchors[flattened_cone] = v;
                        localMinimumBisectors[flattened_cone] = e_bisector_length;
                    }
                    //cout<<"  -"<<D.source(e)<<" is in cone "<<cone<<endl;
                }

                // set the anchor for each negative cone
                for (cone_t i = 0; i < 3; ++i) {
                    if (localMinimumBisectors[i] < INF) {
//                bool same_y = CGAL::compare_y(P[localAnchors[i]], P[w]) == CGAL::EQUAL;
                        Color c = getColor(localAnchors[i], w,
                                           P);//same but less clear -> static_cast<Color>(int(i!=0));
                        anchors[c].emplace_back(localAnchors[i], w);
                        //cout<<"  added "<< (c==Blue?"blue":"white")<<" anchor "<<localAnchors[i]<<"-"<< w<<endl;
                    }
                }

            }
        }

        template<class Set>
        inline bool
        vertexIsEmpty(const VertexDescriptor vertex, const bool isHigherVertex, const Set &AdjacentWhiteConeIsEmpty) {
            bool vertexIsKnown = contains(AdjacentWhiteConeIsEmpty, vertex);
            return // vertex with higher y value needs its negative (second) cones checked
                    (!vertexIsKnown)
                    || (isHigherVertex && AdjacentWhiteConeIsEmpty.at(vertex).negative)
                    || (!isHigherVertex && AdjacentWhiteConeIsEmpty.at(vertex).positive);
        }

        template<class Set>
        inline void
        fillVertex(const VertexDescriptor vertex, const bool isHigherVertex, Set &AdjacentWhiteConeIsEmpty) {
            typedef typename Set::mapped_type PositiveAndNegativeWhiteConesForVertex;
            AdjacentWhiteConeIsEmpty.try_emplace(vertex, PositiveAndNegativeWhiteConesForVertex());
            AdjacentWhiteConeIsEmpty[vertex].negative = !isHigherVertex
                                                        && AdjacentWhiteConeIsEmpty[vertex].negative;
            AdjacentWhiteConeIsEmpty[vertex].positive = isHigherVertex
                                                        && AdjacentWhiteConeIsEmpty[vertex].positive;
        }

        template<class AnchorList, class PointContainer, class EdgeList>
        void addWhiteAnchors(AnchorList &whiteAnchors, const PointContainer &P, EdgeList &A) {
            std::sort(whiteAnchors.begin(), whiteAnchors.end(),
                      [&](const auto &lhs, const auto &rhs) {
                          return bisectorLength(lhs, P) < bisectorLength(rhs, P);
                      }
            );

            // Add all white anchors whose adjacent white cone is empty
            struct PositiveAndNegativeWhiteConesForVertex {    // each vertex can have zero or one of both positive and negative white anchors
                bool positive = true;
                bool negative = true;
            };

            typedef unordered_map<VertexDescriptor,
                    PositiveAndNegativeWhiteConesForVertex> PolarVertexStatusMap;
            PolarVertexStatusMap AdjacentWhiteConeIsEmpty;

            for (auto whiteAnchor : whiteAnchors) {
                //cout<<"WhiteAnchor:"<<whiteAnchor.first<<"-"<<whiteAnchor.second<<"\n";

                VertexDescriptor source = whiteAnchor.first,
                        target = whiteAnchor.second;
                bool sourceIsHigher = CGAL::compare_y(P[source], P[target]) == CGAL::LARGER;

                if (vertexIsEmpty(source, sourceIsHigher, AdjacentWhiteConeIsEmpty)
                    && vertexIsEmpty(target, !sourceIsHigher, AdjacentWhiteConeIsEmpty)) {
                    //cout<<"Both agree\n";
                    A.insert(whiteAnchor);

                    fillVertex(source, sourceIsHigher, AdjacentWhiteConeIsEmpty);
                    fillVertex(target, !sourceIsHigher, AdjacentWhiteConeIsEmpty);
                }
            }
        }

        template<class Triangulation, class EdgeList, class AnchorList, class AdjacencyList>
        void addBlueCanonicalEdges(const Triangulation &D,
                                   const EdgeList &A,
                                   const AnchorList &blueAnchors,
                                   AdjacencyList &S) {
            for (auto blueAnchor : blueAnchors) {
                typedef typename Triangulation::Edge_descriptor EdgeDescriptor;
                vector<EdgeDescriptor> fan;
                //auto source = blueAnchor.first;
                auto target = blueAnchor.second;
                D.fanOfCone(target, 1, fan);

                auto eit = fan.begin();
                auto previousVertex = D.source(*eit);
                // get fan of blueAnchor's target vertex cone 1 (negative/blue)
                for (++eit; eit != fan.end(); ++eit) {
                    auto thisVertex = D.source(*eit);

                    pair<VertexDescriptor, VertexDescriptor> canonicalEdge;
                    bool edgeExists;
                    boost::tie(canonicalEdge, edgeExists)
                            = D.eitherEdge(thisVertex, previousVertex);

                    if (edgeExists && !contains(A, canonicalEdge)) {
                        //cout<<"Canonical edge "<<canonicalEdge.first<<"-"<<canonicalEdge.second<<"\n";
                        S.try_emplace(canonicalEdge.second, typename AdjacencyList::mapped_type());
                        S[canonicalEdge.second].insert(canonicalEdge.first);
                    }
                    previousVertex = thisVertex;
                }
            }
        }

        template<class Triangulation, class EdgeList, class PointContainer, class AnchorList, class AdjacencyList>
        void addWhiteCanonicalEdges(const Triangulation &D,
                                    const EdgeList &A,
                                    const PointContainer &P,
                                    const AnchorList &whiteAnchors,
                                    AdjacencyList &S) {
            for (auto whiteAnchor : whiteAnchors) {
                typedef typename Triangulation::Edge_descriptor EdgeDescriptor;

                auto source = whiteAnchor.first;
                auto target = whiteAnchor.second;
                cone_t cone = D.getCone(target, source);
                vector<EdgeDescriptor> fan;
                D.fanOfCone(target, cone, fan);

                auto boundaryVertex = D.source(fan.back()),
                        previousVertex = D.source(fan.front());
                int direction = 1;

                if ((cone + 1) % 3 == 0) { //white side of anchor is CW
                    direction = -1;
                    swap(previousVertex, boundaryVertex);
                }

                auto eit = fan.begin();
                while (D.source(*eit) != source)
                    ++eit;

                previousVertex = D.source(*eit);

                for (eit += direction; previousVertex != boundaryVertex; eit += direction) {
                    auto thisVertex = D.source(*eit);
                    //cout<<"Considering white canonical edge "<<previousVertex<<"-"<<thisVertex<<endl;

                    pair<VertexDescriptor, VertexDescriptor> canonicalEdge;
                    bool edgeExists;
                    boost::tie(canonicalEdge, edgeExists)
                            = D.eitherEdge(thisVertex, previousVertex);

                    // For white cones, add can. edge if it is white and its anchor is not in A
                    if (edgeExists
                        && getColor(canonicalEdge, P) == White
                        && !contains(A, canonicalEdge)) {
                        //cout<<"Canonical edge "<<canonicalEdge.first<<"-"<<canonicalEdge.second<<"\n";
                        S.try_emplace(canonicalEdge.second, typename AdjacencyList::mapped_type());
                        S[canonicalEdge.second].insert(canonicalEdge.first);
                    }
                    previousVertex = thisVertex;
                }
            }
        }

        template<class AdjacencyList>
        void createShortcut(const VertexDescriptor &p, const VertexDescriptor &q, const VertexDescriptor &r,
                            AdjacencyList &adj) {
            adj[q].erase(p);
            adj[q].erase(r);

            adj.try_emplace(r, typename AdjacencyList::mapped_type());
            adj[r].insert(p);
        }

        template<class Triangulation, class AdjacencyList>
        void addBlueShortcuts(const Triangulation &D, AdjacencyList &S_not_A) {
            for (auto v : S_not_A) {
                // check if there are two edges that share a target
                //cout<<"  |"<<v.first<<"|="<<v.second.size()<<"\n";
                if (v.second.size() > 1) {
                    auto p = *v.second.begin(),
                            q = v.first,
                            r = *(--v.second.end());

//            cout<< "    Blue shortcut edge "<<*v.second.begin()
//                <<"-"<<v.first<<"-"<<*(--v.second.end())<< " --> "
//                <<*v.second.begin()<<"-"<<*(--v.second.end())<<"\n";
                    if(!(D.edgeExists(make_pair(p, q)) && D.edgeExists(make_pair(r, q))))
                        cout<<p<<" "<<q<<" "<<r<<"\n";
                    //assert(D.edgeExists(make_pair(p, q)) && D.edgeExists(make_pair(r, q)));
                    createShortcut(p, q, r, S_not_A);
                }
            }
        }


        template<class Triangulation, class PointContainer, class AnchorList, class AdjacencyList>
        void addWhiteShortcuts(const Triangulation &D,
                               const PointContainer &P,
                               const AnchorList &whiteAnchors,
                               AdjacencyList &S) {
            typedef typename Triangulation::Edge_descriptor EdgeDescriptor;

            for (auto whiteAnchor : whiteAnchors) {
                auto source = whiteAnchor.first;
                auto target = whiteAnchor.second;
                //cout<<"Looking for shortcuts of white anchor "<<source<<"-"<<target<<endl;
                cone_t cone = D.getCone(target, source);
                vector<EdgeDescriptor> fan;
                D.fanOfCone(target, cone, fan);

                auto boundaryEdgeOnWhiteSide = --(fan.end()),
                        boundaryEdgeOnBlueSide = fan.begin();
                int direction = 1;

                if ((cone + 1) % 3 == 0) { //white side of anchor is CW
                    direction = -1;
                    swap(boundaryEdgeOnWhiteSide,
                         boundaryEdgeOnBlueSide);
                }
                //cout<<"direction="<<direction<<endl;

                auto i = std::find_if(fan.begin(), fan.end(),
                                      [&](const auto &e) {
                                          return D.source(e) == source;
                                      });

                if (i != boundaryEdgeOnWhiteSide) {
                    auto iPlusOne = i + direction,
                            iLessOne = i;
//            while( D.source(*i) != source )
//                ++i;

                    VertexDescriptor v_i = D.source(*i),
                            v_iPlusOne = v_i,
                            v_iLessOne = v_i,
                            boundaryVertex = D.source(*boundaryEdgeOnWhiteSide);

                    //cout<<"index in fan="<<(i-fan.begin())<<endl;

                    for (; i != boundaryEdgeOnWhiteSide;) {
                        v_iLessOne = v_i,
                        v_i = v_iPlusOne,
                        v_iPlusOne = D.source(*iPlusOne);
                        //cout<<"iLessOne="<<v_iLessOne<<" i="<<v_i<<" iPlusOne="<<v_iPlusOne<<" "<<endl;

                        if (getColor(v_iPlusOne, v_i, P) == White) {
                            //cout<<"color is white"<<endl;
                            iLessOne = i,
                            i = iPlusOne,
                            iPlusOne += direction;
                        } else {
                            //cout<<"color is blue"<<endl;
                            auto j = iPlusOne;

                            auto j_vertex = D.source(*j),
                                    j_prev = D.source(*i),
                                    min_j = j,
                                    max_j = j;

                            auto j_angle = getAngle(P[target], P[v_i], P[j_vertex]),
                                    maxAngle = j_angle,
                                    minAngle = j_angle;

                            for (; j_prev !=
                                   boundaryVertex; j += direction) {    // j is in a white cone AND all vertices from i+1 to boundaryVertex are on the same side of i,j
                                j_vertex = D.source(*j);
                                if (getColor(v_i, j_vertex, P) == White
                                    && getConePolarity(v_i, j_vertex, P) == Negative) {
                                    j_angle = getAngle(P[target], P[v_i], P[j_vertex]);

                                    if (j_angle > maxAngle) {
                                        maxAngle = j_angle;
                                        max_j = j;
                                    } else if (j_angle < minAngle) {
                                        minAngle = j_angle;
                                        min_j = j;
                                    }
                                }
                                j_prev = j_vertex;
                            }
                            auto min_dist = min_j - i,
                                    max_dist = max_j - i;

                            auto best_j = min_dist < max_dist ? min_j : max_j;
                            j_vertex = D.source(*best_j);
                            j_prev = D.source(*(best_j - direction));

                            //cout<<"White shortcut edge "<<j_vertex<<"-"<<v_i<<"\n";
                            S.try_emplace(v_i, typename AdjacencyList::mapped_type());
                            S[v_i].insert(j_vertex);

                            if (S.find(j_prev) != S.end()) {
                                S[j_prev].erase(j_vertex);
                            }

                            i = best_j;
                            iLessOne = i - direction;
                            iPlusOne = i + direction;
                        }
                    }
                }
            }
        }

    } // namespace KPT2017


// Main algorithm.
    template<typename RandomAccessIterator, typename OutputIterator>
    void KPT2017(RandomAccessIterator pointsBegin,
                 RandomAccessIterator pointsEnd,
                 OutputIterator result) {
        using namespace kpt2017;

        vector<Point> P(pointsBegin, pointsEnd);

        TDDelaunayGraph D(P.begin(), P.end());
        {
            //Timer tim;

            map<Color, vector<Edge>> Anchors;
            findAnchors(Anchors, D, P);

            auto AnchorComp = [&P](const Edge &lhs, const Edge &rhs) {
                return lhs.second < rhs.second
                       || (lhs.second == rhs.second
                           && getCone(lhs.second, lhs.first, P) < getCone(rhs.second, rhs.first, P));
            };
            // Step 1. add all blue anchors to A
            set<Edge, decltype(AnchorComp)> A(
                    Anchors[Blue].begin(),
                    Anchors[Blue].end(),
                    AnchorComp);

            // Step 2.
            addWhiteAnchors(Anchors[White], P, A);

            // Step 3.
            typedef map<index_t, set<index_t>> AdjacencyListMap;
            AdjacencyListMap S_not_A;
            // Add to S every canonical edge in negative blue cones (cone 1) if the edge isn'stretchFactor in A
            addBlueCanonicalEdges(D, A, Anchors[Blue], S_not_A);

            //Step 4.
            addBlueShortcuts(D, S_not_A);

            //Step 5.
            addWhiteCanonicalEdges(D, A, P, Anchors[White], S_not_A);

            //Step 6.
            addWhiteShortcuts(D, P, Anchors[White], S_not_A);




            // Since we didn'stretchFactor set S to A in step 3, we need to combine them now
            set<Edge> S(A.begin(), A.end());

            //cout<<"Adding S_not_A to S...\n";
            for (const auto &v : S_not_A) {
                if (!v.second.empty()) {
                    for (auto u : v.second) {
                        S.emplace(u, v.first);
                    }
                }
            }

            // Send resultant graph to output iterator
            std::copy(S.begin(), S.end(), result);
//        for(auto e : S)
//        {
//            *result = e;
//            ++result;
////            *result = reverse_pair(e);
////            ++result;
//        }





            // START PRINTER NONSENSE
//        if(printLog)
//        {
//            vector<pair<Point,Point>> edgeList;
//
//            for(auto e : S)
//            {
//                edgeList.emplace_back(P.at(e.first), P.at(e.second));
//            }
//
//            GraphPrinter printer(0.01);
//            GraphPrinter::OptionsList options;
//
//            options = {
//                {"color", printer.inactiveEdgeColor},
//                {"line width", to_string(printer.inactiveEdgeWidth)}
//            };
//            printer.drawEdgesOfHalfTheta(D, options);
//
//            options = { // active edge options
//                {"color", printer.activeEdgeColor},
//                {"line width", to_string(printer.activeEdgeWidth)}
//            };
//            printer.drawEdges(S.begin(), S.end(), P, options);
//
//            options = {
//                {"vertex", make_optional(to_string(printer.vertexRadius))}, // vertex width
//                {"color", make_optional(printer.backgroundColor)}, // text color
//                {"fill", make_optional(printer.activeVertexColor)}, // vertex color
//                {"line width", make_optional(to_string(0))} // vertex border (same color as text)
//            };
//            GraphPrinter::OptionsList borderOptions = {
//                {"border", make_optional(to_string(printer.vertexRadius))}, // choose shape of vertex
//                {"color", printer.activeEdgeColor}, // additional border color
//                {"line width", to_string(printer.inactiveEdgeWidth)}, // additional border width
//            };
//            printer.drawVerticesWithInfo(D.points_begin(), D.points_end(), options, borderOptions);
//
//            printer.print("KPT2017");
//            cout << "\n";
//        }
            // END PRINTER NONSENSE
        }
    } // function KPT2017

} // namespace planespanners

#endif // GSNUNF_KPT2017_H