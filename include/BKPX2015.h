//Needs optimizing currently testing.
#ifndef GSNUNF_BKPX2015_H
#define GSNUNF_BKPX2015_H

#include <array>
#include <iostream>
#include <fstream>
#include <list>
#include <cassert>
#include <string>
#include <limits>
#include <cmath>         // ceil, floor
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles
#include <boost/functional/hash.hpp> // size_t pair hash : used in Yao_inf_4

#include <CGAL/algorithm.h> //
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/circulator.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>

#include "GeometricSpannerPrinter.h"
#include "metrics.h"

namespace gsnunf {

    using namespace std;

    namespace bkpx2015 {

    // CGAL objects

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<K,CGAL::Field_with_sqrt_tag>  Gt;
    typedef size_t id_type;

    template<typename T>
    struct InfoConvert
    {
        typedef T Info;
        typedef const Info& result_type;

        inline const Info& operator()(const Info& info0, bool) const {
            // just return the info of the supporting segment
            return info0;
        }
        inline const Info& operator()(const Info& info0, const Info& , bool) const {
            // just return the info of the supporting segment
            return info0;
        }
    };
    // functor that defines how to merge color info when a site (either
    // point or segment) corresponds to point(s) on plane belonging to
    // more than one input site
    template<typename T>
    struct InfoMerge
    {
        typedef T    Info;
        typedef Info result_type;

        inline Info operator()(const Info& info0, const Info& info1) const {
            // just return the info of the supporting segment
            return info0;
        }
    };

    typedef CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt,
                                                                    id_type,
                                                                    InfoConvert<id_type>,
                                                                    InfoMerge<id_type>> ST;
    typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt,ST> SDG2;
    typedef SDG2::Site_2                   Site_2;
    typedef Site_2::Point_2                Point_2;
    typedef SDG2::Vertex_circulator        Vertex_circulator;
    typedef SDG2::Edge                     LinfEdge;
    typedef SDG2::Vertex_handle            Vertex_handle;

    typedef CGAL::Spatial_sort_traits_adapter_2<K,
        CGAL::Pointer_property_map<Point_2>::type > Search_traits_2;


    enum AnchorType { None, Weak, Strong, StrongSelected, WeakSelected, StartOddChain };


    // project objects
    typedef vector<pair<Vertex_handle, Vertex_handle>>                  fanCones;
    typedef vector<pair<Vertex_handle, size_t>>                         yaoCones;
    typedef vector<size_t>                                              numYaoEdges;
    typedef vector<pair<Vertex_handle, AnchorType>>                     anchorCones;
    typedef pair<Vertex_handle, Vertex_handle>                          spannerEdge;
    typedef vector<vector<spannerEdge>>                                 spannerCones;

    // inline functions required for bkpx2015

    inline size_t getSingleCone( const Vertex_handle &u, const Vertex_handle &v) {
//       bool X = (v->site().point().x()- u->site().point().x()) > 0;
//       bool Y = (v->site().point().y()- u->site().point().y()) > 0;
//       return !Y*2 + (!X*Y+X*!Y)*1;
        const number_t alpha = PI/2;
        const Point_2 refPoint( u->site().point().x(), u->site().point().y() + 1 );
        //Point refPoint(h[p]->point().x(), h[p] ->point().y() + 1);

        number_t theta = get_angle(refPoint, u->site().point(), v->site().point());

        size_t cone = (theta / alpha);

        return cone;

    }

    // get the cone of v wrt u (where u is at the center)
    inline size_t getCone( const Vertex_handle &u, const Vertex_handle &v) {

        if (u > v) {
            return getSingleCone(u, v);
        }

        else {
            return (getSingleCone(v, u) + 2) % 4;
        }

    }

    // add yao edges
    inline void addYaoEdges(vector<yaoCones> &yaoEdges,
                            vector<fanCones> &pointFans,
                            vector<numYaoEdges> &yaoEdgeCount,
                            const vector<Vertex_handle> &handles,
                            const SDG2 &sdg)
    {

        Vertex_circulator Circ = sdg.incident_vertices(handles[0]); // default value
        vector<number_t> distances(4);


        size_t cone = 0;
        size_t index = 0;

        for (auto point : handles) {

            Circ = sdg.incident_vertices(point);
            distances = {INF, INF, INF, INF};
            index = point->storage_site().info();

            yaoCones &edges = yaoEdges[index];
            fanCones &fans = pointFans[index];

            while( sdg.is_infinite(Circ) )
                ++Circ;

            size_t previousCone = getCone(point, Circ);

            auto startPoint = Circ++;

            while (!sdg.is_infinite(Circ) && getCone(point, Circ) == previousCone && Circ != startPoint)
                ++Circ;

            while (sdg.is_infinite(Circ))
                ++Circ;

            auto endpoint = Circ;
            cone = getCone(point, Circ);
            previousCone = 10; // invalid value, default as 10

            size_t point_id = point->storage_site().info();

            do {

                if (!sdg.is_infinite(Circ))
                {
                    cone = getCone(point,Circ);
                    ++(edges[cone].second);

                    if (cone != previousCone)
                        fans[cone].first = Circ;

                    fans[cone].second = Circ;
                    previousCone = cone;


            //        double proposedDistance = std::max(std::abs(point->site().point().x() - Circ->site().point().x()), std::abs(point->site().point().y() - Circ->site().point().y()));

                    number_t proposedDistance = CGAL::l_infinity_distance(point->site().point(), Circ->site().point());

//                    if (proposedDistance == distances[cone]) {
//                        cout << point->storage_site().info() << " " << Circ->storage_site().info() << " tied with " << point->storage_site().info() << " " << edges[cone].first->storage_site().info() << endl;
//                    }

                    if (proposedDistance < distances[cone])
                    {
                        distances[cone] = proposedDistance;
                        edges[cone].first = Circ;
                    }

                }

            } while(++Circ != endpoint); // finished determining the Yao edges + how many points are in fan of u's cone i

            for (size_t cone = 0; cone < 4; cone++) {

                if (yaoEdges[index][cone].second == 0)
                    continue;

                auto v = yaoEdges[index][cone].first;
                size_t v_id = v->storage_site().info();

                if (yaoEdges[v_id][(cone+2)%4].first != point)
                {
                    ++(yaoEdgeCount[index][cone]);
                    ++(yaoEdgeCount[v_id][(cone+2)%4]);
                }
            }
        } // finished moving through points
    } // yaoEdges are added


    inline void determineAnchors( vector<anchorCones> &anchorEdges,
                                  vector<yaoCones> &yaoEdges,
                                  vector<fanCones> &pointFans,
                                  vector<numYaoEdges> &yaoEdgeCount,
                                  vector<Vertex_handle> &handles,
                                  SDG2 &sdg)
    {

        for (auto u : handles) {

            size_t u_id = u->storage_site().info();

            anchorCones &anchors = anchorEdges[u_id];

            for (size_t cone = 0; cone < 4; cone++) {

                // no neighbors in cone
                if (yaoEdges[u_id][cone].second == 0) {
                    continue;
                }

                auto v = yaoEdges[u_id][cone].first;

                // only one neighbor in cone --> check if mutually single
                if (yaoEdgeCount[u_id][cone] == 1)
                {
                    if (yaoEdgeCount[v->storage_site().info()][(cone+2)%4] == 1)
                    {
                        anchorEdges[u_id][cone].first = v;
                        anchorEdges[u_id][cone].second = Weak;
                    }
                }

                if (yaoEdgeCount[u_id][cone] >= 2)
                {

                    auto v1 = pointFans[u_id][cone].first;
                    auto vk = pointFans[u_id][cone].second;
                    Vertex_circulator current = sdg.incident_vertices(u);
                    while (current != v1) {++current;}
                    size_t numNeighbors = yaoEdges[u_id][cone].second;
                    size_t lcount = 1;

                    while (current != v)
                    {
                        ++lcount;
                        ++current;
                    }

                    // establish vlower and vhigher and recalibrate Circ accordingly
                    auto vlower = current;
                    --vlower;
                    auto vhigher = current;
                    ++vhigher;

                    bool ccwCanonical = (lcount >= 2 && !(sdg.is_infinite(vlower)) && (yaoEdges[vlower->storage_site().info()][getCone(vlower, v)].first == v &&
                        yaoEdges[v->storage_site().info()][getCone(v, vlower)].first != vlower));

                    bool cwCanonical = (!ccwCanonical && lcount <= (numNeighbors - 1) && !(sdg.is_infinite(vhigher)) && (yaoEdges[vhigher->storage_site().info()][getCone(vhigher, v)].first == v) &&
                        yaoEdges[v->storage_site().info()][getCone(v, vhigher)].first != vhigher);

                    bool inCanonical = ccwCanonical || cwCanonical;
                    size_t position = lcount;
                    auto previous = current;
                    auto crown = previous;

                    int direction = 1 - 2 * ccwCanonical;

                    if (inCanonical) {
                        while (inCanonical)
                        {
                            previous = current;

                            if (ccwCanonical)
                                --current;
                            else
                                ++current;

                            position += direction;

                            bool inFan = !((ccwCanonical && position < 1) || (cwCanonical && position > numNeighbors));
                            bool unidirectional = (!(sdg.is_infinite(current)) && yaoEdges[current->storage_site().info()][getCone(current, previous)].first == previous &&
                                yaoEdges[previous->storage_site().info()][getCone(previous, current)].first != current);

                            inCanonical = inFan && unidirectional;

                            bool yaoConnected = inCanonical && (yaoEdges[u_id][cone].first == current || yaoEdges[current->storage_site().info()][(cone+2)%4].first == u);

                            if (yaoConnected)
                                crown = current;

                        }

                        assert(yaoEdges[u_id][cone].first == crown || yaoEdges[crown->storage_site().info()][(cone+2)%4].first == u);

                        anchorEdges[u_id][cone].first = crown;
                        anchorEdges[u_id][cone].second = Weak;
                    }
                    else
                    {
                        anchorEdges[u_id][cone].first = v;
                        anchorEdges[u_id][cone].second = Weak;
                    }
                } // when there are more than 2 edges!
            }
        }

        for (auto u : handles) {
            for (size_t cone = 0; cone < 4; cone++)
            {
                size_t u_id = u->storage_site().info();

                Vertex_handle v = anchorEdges[u_id][cone].first;

                if( !sdg.is_infinite(v)
                 && (   anchorEdges[v->storage_site().info()][(cone+2)%4].first == u
                     || sdg.is_infinite(anchorEdges[v->storage_site().info()][(cone+2)%4].first) ) )
                {
                    anchorEdges[u_id][cone].second = Strong;
                }
            }

        } // strong anchors are identified


        // now it is time to select anchors
        for (auto u : handles)
        {
            size_t u_id = u->storage_site().info();

            for (size_t cone = 0; cone < 4; cone++)
            {
                if (sdg.is_infinite(anchorEdges[u_id][cone].first))
                    continue;

                if (anchorEdges[u_id][cone].second == Strong)
                    anchorEdges[u_id][cone].second = StrongSelected;

                if (anchorEdges[u_id][cone].second == Weak)
                {
                    Vertex_circulator Circ = sdg.incident_vertices(u);
                    auto v1 = pointFans[u_id][cone].first,
                         vk = pointFans[u_id][cone].second;
                    while (Circ != v1) {++Circ;}

                    bool found = false;

                    do {
                        found = anchorEdges[Circ->storage_site().info()][(cone+2)%4].first == u &&
                            (anchorEdges[Circ->storage_site().info()][(cone+2)%4].second == Weak ||
                             anchorEdges[Circ->storage_site().info()][(cone+2)%4].second == WeakSelected);
                    } while (!found && Circ++ != vk);

                    if (found) {continue;}

//                    cout << endl << endl << u_id << " is the beginning of the weak anchor chain. " << endl;

                    vector<Vertex_handle> visited;
                    auto previous = u;
                    auto current = anchorEdges[u_id][cone].first;
                    size_t previous_id = u_id;
                    size_t currentCone = cone;
                    size_t localCone = (cone+2)%4;

                    do {

                        visited.push_back(previous);

                        previous = current;
                        previous_id = previous->storage_site().info();

                        currentCone = (visited.size() % 2) ? localCone : cone;
                        current = anchorEdges[previous_id][currentCone].first;

                    } while (anchorEdges[previous_id][currentCone].second == Weak);

                    bool oddChain = (visited.size() % 2);

                    size_t position = 0;
                    current = visited.at(position);
                    size_t current_id = current->storage_site().info();

                    currentCone = cone;

                    if (oddChain)
                    {
                        anchorEdges[current_id][cone].second = StartOddChain;
                        currentCone = localCone;
                        ++position;
                    }

                    for (size_t i = position; i < visited.size(); i += 2)
                    {
                        current = visited.at(i);
                        current_id = current->storage_site().info();
                        anchorEdges[current_id][currentCone].second = WeakSelected;
                    }
                }
            }
        }
    } // function complete

    inline bool inEdgeList(const vector<spannerEdge> &edgeList,
                           const Vertex_handle u,
                           const Vertex_handle v)
    {
        for (auto edge : edgeList)
            if ((u == edge.first || u == edge.second) && (v == edge.first || v == edge.second))
                return true;

        return false;
    }


    inline void degreeEightSpanner(vector<spannerCones> &H8,
                                 vector<anchorCones> &anchorEdges,
                                 vector<yaoCones> &yaoEdges,
                                 vector<fanCones> &pointFans,
                                 vector<numYaoEdges> &yaoEdgeCount,
                                 vector<Vertex_handle> &handles,
                                 SDG2 &sdg)
    {
        // put the edges conal vector
        for (auto w : handles)
        {
            spannerCones edges(4);
            H8[w->storage_site().info()] = edges;
        }

        for (auto w : handles) {

            size_t w_id = w->storage_site().info();

            for (size_t cone = 0; cone < 4; cone++) {

                if (yaoEdges[w_id][cone].second == 0) {continue;}

                if (anchorEdges[w_id][cone].second == StrongSelected || anchorEdges[w_id][cone].second == WeakSelected) {

                    auto v = anchorEdges[w_id][cone].first;
                    bool original = !inEdgeList(H8[w_id][cone], w, v);

                    if (original) {
                        spannerEdge edge = std::make_pair(w, v);
                        H8[w_id][cone].push_back(edge);
                        H8[v->storage_site().info()][(cone+2)%4].push_back(edge);
                    }

                }

                if (yaoEdges[w_id][cone].second > 1) {

                    Vertex_circulator current = sdg.incident_vertices(w);
                    while (current != pointFans[w_id][cone].first) {++current;}
                    auto previous = current;
                    ++current;

                    size_t position = 1;
                    size_t total = yaoEdges[w_id][cone].second;

                    while (position < total) {

                        size_t cwCone = getCone(current, previous);
                        size_t ccwCone = getCone(previous, current);

                        size_t current_id = current->storage_site().info();
                        size_t previous_id = previous->storage_site().info();

                        bool yaoConnected = (yaoEdges[w_id][cone].first == current || yaoEdges[current_id][(cone+2)%4].first == w)
                                         && (yaoEdges[w_id][cone].first == previous || yaoEdges[previous_id][(cone+2)%4].first == w);

                        bool uniCW = yaoEdges[current_id][cwCone].first == previous
                        && yaoEdges[previous_id][ccwCone].first != current;

                        bool uniCCW = yaoEdges[current_id][cwCone].first != previous
                        && yaoEdges[previous_id][ccwCone].first == current;

                        bool nonAnchor = anchorEdges[current_id][cwCone].first != previous
                            && anchorEdges[previous_id][ccwCone].first != current;

                        bool addEdge = (uniCW || uniCCW) && yaoConnected;

                        if (yaoConnected && position == 1 && uniCW && nonAnchor) {

                            bool dual = yaoEdgeCount[w_id][cone] > 1 && yaoEdgeCount[previous_id][(cone+2)%4] > 1;

                            bool startOdd = anchorEdges[previous_id][(cone+2)%4].first == w
                                && anchorEdges[previous_id][(cone+2)%4].second == StartOddChain;

                            addEdge = !(dual && !startOdd);

                        }

                        if (yaoConnected && position == total-1 && uniCCW && nonAnchor) {

                            bool dual = yaoEdgeCount[w_id][cone] > 1 && yaoEdgeCount[current_id][(cone+2)%4] > 1;

                            bool startOdd = anchorEdges[current_id][(cone+2)%4].first == w
                                && anchorEdges[current_id][(cone+2)%4].second == StartOddChain;

                            addEdge = !(dual && !startOdd);

                        }

                        if (addEdge) {

                            auto source = previous;
                            auto target = current;
                            size_t directedCone = ccwCone;

                            if (uniCW) {

                                source = current;
                                target = previous;
                                directedCone = cwCone;

                            }

                            size_t source_id = source->storage_site().info();
                            size_t target_id = target->storage_site().info();

                            if (!inEdgeList(H8[source_id][directedCone], source, target)) {

                                H8[source_id][directedCone].emplace_back(source, target);

                                if (nonAnchor) {
                                    H8[target_id][(cone+2)%4].emplace_back(source, target);
                                }

                                else {
                                    H8[target_id][(directedCone+2)%4].emplace_back(source, target);
                                }

                            }

                        }

                        ++position;
                        previous = current;
                        ++current;

                    }
                }
            }
        }
    }


    inline void processSpanner(vector<spannerCones> &H8,
                               const vector<anchorCones> &anchorEdges,
                               const vector<yaoCones> &yaoEdges,
                               const vector<fanCones> &pointFans,
                               const vector<numYaoEdges> &yaoEdgeCount,
                               const vector<Vertex_handle> &handles,
                               const SDG2 &sdg)
    {

        for (auto u : handles) {

            size_t u_id = u->storage_site().info();

            for (size_t cone = 0; cone < 4; cone++) {

                size_t charge = H8[u_id][cone].size();

                if (charge == 1) {

                    spannerEdge edge = H8[u_id][cone][0];

                    auto source = edge.first;
                    auto target = edge.second;

                    size_t source_id = source->storage_site().info();
                    size_t target_id = target->storage_site().info();

                    bool unidirectional = yaoEdges[source_id][cone].first == target && yaoEdges[target_id][(cone+2)%4].first != source;
                    bool nonanchor = anchorEdges[source_id][cone].first != target && anchorEdges[target_id][(cone+2)%4].first != source;

                    if (!(unidirectional && nonanchor)) { continue; }

//                    cout << endl << endl << "duplicate edge chain begins at " << source_id << endl;

                    for (size_t i = 1; i <= 3; i += 2) {

                        size_t localCone = (cone+i)%4;

                        if (inEdgeList(H8[target_id][localCone], source, target) && H8[target_id][localCone].size() == 2) {

                            vector<Vertex_handle> visited;
                            auto previous = source;
                            size_t previous_id = previous->storage_site().info();
                            auto current = target;
                            size_t current_id = current->storage_site().info();
                            size_t currentCone = cone;
                            size_t indicatorCone = localCone;

                            do  {

                                visited.push_back(previous);

                                previous = current;
                                previous_id = previous->storage_site().info();

                                currentCone = (visited.size() % 2) ? localCone : cone;
                                current = yaoEdges[previous_id][currentCone].first;

                            } while (inEdgeList(H8[previous_id][currentCone], previous, current) && H8[previous_id][currentCone].size() == 2);

                            size_t total = visited.size()-1;

                            // it's time to process edges

                            for (size_t i = 0; i < total; i += 2) {

                                size_t l = (total-i);

                                auto finish = visited.at(l);
                                auto start = visited.at(l-1);

                                size_t finish_id = finish->storage_site().info();
                                size_t start_id = start->storage_site().info();

                                spannerEdge edge = std::make_pair(start, finish);

                                for (size_t i = 0; i < 4; i++) {

                                    while (inEdgeList(H8[start_id][i], start, finish))
                                        H8[start_id][i].erase(find(H8[start_id][i].begin(), H8[start_id][i].end(), edge));

                                    while (inEdgeList(H8[finish_id][i], start, finish))
                                        H8[finish_id][i].erase(find(H8[finish_id][i].begin(), H8[finish_id][i].end(), edge));

                                }
                            }
                        }
                    }
                }
            }
        }

        for (auto u : handles) {

            size_t u_id = u->storage_site().info();
            Vertex_circulator middle = sdg.incident_vertices(u);

            for (size_t cone = 0; cone < 4; cone++) {

                if (yaoEdges[u_id][cone].second >= 3) {

                    auto v1 = pointFans[u_id][cone].first,
                         vk = pointFans[u_id][cone].second;

                    while (middle != v1) {++middle;}
                    ++middle;

                    // will continue so long as the charge is not 2 on middle in cone i+2 and middle is not vk

                    size_t total = yaoEdges[u_id][cone].second,
                           position = 1;

                    bool found = false;

                    while (!found && position < total-1) {

                        if (H8[middle->storage_site().info()][(cone+2)%4].size() == 2 && middle != vk) {
                            auto previous = middle,
                            next = middle;

                            --previous;
                            ++next;

                            size_t previous_id = previous->storage_site().info();
                            size_t middle_id = middle->storage_site().info();
                            size_t next_id = next->storage_site().info();

                            size_t previousCone = getCone(previous, middle);
                            size_t nextCone = getCone(next, middle);

                            bool edgePair = (yaoEdges[previous_id][previousCone].first == middle && yaoEdges[middle_id][(previousCone+2)%4].first != previous
                                      && anchorEdges[previous_id][previousCone].first != middle && anchorEdges[middle_id][(previousCone+2)%4].first != previous)
                                      && (yaoEdges[next_id][nextCone].first == middle && yaoEdges[middle_id][(nextCone+2)%4].first != next
                                      && anchorEdges[next_id][nextCone].first != middle && anchorEdges[middle_id][(nextCone+2)%4].first != next);

                            if (edgePair) {

                            //                        cout << endl << "removing <" << previous_id << ", " << middle_id << "> and <" << next_id << ", " << middle_id << "> " << endl;
                            //                        cout << "adding <" << previous_id << ", " << next_id << "> " << endl;

                            spannerEdge previousEdge = std::make_pair(previous, middle);
                            spannerEdge nextEdge = std::make_pair(next, middle);

                            //                         remove previousEdge
                            while (inEdgeList(H8[previous_id][previousCone], previous, middle)) {
                                H8[previous_id][previousCone].erase(find(H8[previous_id][previousCone].begin(), H8[previous_id][previousCone].end(), previousEdge));
                            }

                            while (inEdgeList(H8[middle_id][(cone+2)%4], previous, middle)) {
                                H8[middle_id][(cone+2)%4].erase(find(H8[middle_id][(cone+2)%4].begin(), H8[middle_id][(cone+2)%4].end(), previousEdge));
                            }


                            // remove nextEdge
                            while (inEdgeList(H8[next_id][nextCone], next, middle)) {
                                H8[next_id][nextCone].erase(find(H8[next_id][nextCone].begin(), H8[next_id][nextCone].end(), nextEdge));
                            }

                            while (inEdgeList(H8[middle_id][(cone+2)%4], next, middle)) {
                                H8[middle_id][(cone+2)%4].erase(find(H8[middle_id][(cone+2)%4].begin(), H8[middle_id][(cone+2)%4].end(), nextEdge));
                            }

                            //                        cout << "adding shortcut: <" << previous_id << ", " << next_id << "> " << endl;

                            spannerEdge shortcut = std::make_pair(previous, next);
                            H8[previous_id][previousCone].push_back(shortcut);
                            H8[next_id][nextCone].push_back(shortcut);

                            }
                        }

                        ++position;
                        ++middle;

                    }
                }
            }
        }
    }
}



template<typename RandomAccessIterator, typename OutputIterator>
void BKPX2015(RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false) {

    using namespace bkpx2015;

    using bkpx2015::Vertex_handle, bkpx2015::Vertex_circulator;

    // construct Linf Delaunay triangulation
    vector<Point_2> P(pointsBegin, pointsEnd);
    const size_t n = P.size();

    //    N is the number of vertices in the delaunay triangulation.
    if(n > SIZE_T_MAX - 1) return;

    SDG2 sdg;
    Site_2 site;
    id_type id = 0;

    vector<size_t> indices(n);

    for (size_t i = 0; i < n; i++) {
        indices.at(i) = i;
    }

    CGAL::spatial_sort( indices.begin(),
                      indices.end(),
                      Search_traits_2(CGAL::make_property_map(P)) );

   // store the vertex handles
    vector<Vertex_handle> handles(n);

    for (size_t i = 0; i < n; i++) {
        size_t index = indices.at(i);
        Point_2 point = P.at(index);
        site = Site_2::construct_site_2(point);
        Vertex_handle v = sdg.insert(site, index);
        handles.at(index) = v;
    }



    //construct YaoEdges
    vector<yaoCones> yaoEdges(n, yaoCones(4));
    vector<fanCones> pointFans(n, fanCones(4));
    vector<numYaoEdges> yaoEdgeCount(n, numYaoEdges(4));
    addYaoEdges(yaoEdges, pointFans, yaoEdgeCount, handles, sdg);

    // identify the anchors
    vector<anchorCones> anchorEdges(n, anchorCones(4, std::make_pair(sdg.infinite_vertex(), None)));
    determineAnchors(anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, sdg);

   // construct H8 --> degree 8 spanner
    vector<spannerCones> H8(n, spannerCones(4));
    degreeEightSpanner(H8, anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, sdg);


    processSpanner(H8, anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, sdg);

    vector<size_tPair> edgeList;

    for (auto u : handles) {
        for (size_t cone = 0; cone < 4; cone++) {
            for (auto edge : H8[u->storage_site().info()][cone]) {
                edgeList.emplace_back( (edge.first)->storage_site().info(),
                                       (edge.second)->storage_site().info() );
            }
        }
    }

    // Send resultant graph to output iterator
    for(auto e : edgeList) {
        // Edge list is only needed for printing. Remove for production.
        //edgeList.emplace_back(handles.at(e.first)->point(), handles.at(e.second)->point());

        *result = e;
        ++result;
    }

    // START PRINTER NONSENSE
//    if(printLog) {
//        GraphPrinter printer(1.25); // argument number is scaling factor --> manipulate based on size of point set
//        GraphPrinter::OptionsList options;
//
//        options = {
//            {"color", printer.inactiveEdgeColor},
//            {"line width", to_string(printer.inactiveEdgeWidth)}
//        };
////        printer.drawEdgesOfSDG(sdg, options);
//
//        options = { // active edge options
//            {"color", printer.activeEdgeColor},
//            {"line width", to_string(printer.activeEdgeWidth)}
//        };
//
//        vector<pair<Point_2, Point_2>> pointEdgeList;
//        pointEdgeList.reserve(edgeList.size());
//
//        for (auto e : edgeList) {
//
//            pointEdgeList.emplace_back(P.at(e.first), P.at(e.second));
//
//        }
//
//        printer.drawEdges(edgeList.begin(), edgeList.end(), P, options);
//
//
//        options = {
//            {"vertex", make_optional(to_string(printer.vertexRadius))}, // vertex width
//            {"color", make_optional(printer.backgroundColor)}, // text color
//            {"fill", make_optional(printer.activeVertexColor)}, // vertex color
//            {"line width", make_optional(to_string(0))} // vertex border (same color as text)
//        };
//        GraphPrinter::OptionsList borderOptions = {
//            {"border", make_optional(to_string(printer.vertexRadius))}, // choose shape of vertex
//            {"color", printer.activeEdgeColor}, // additional border color
//            {"line width", to_string(printer.inactiveEdgeWidth)}, // additional border width
//        };
//        printer.drawVerticesWithInfoSDG(sdg, options, borderOptions);
//
//        printer.print("BKPX2015");
//        cout << "\n";
//    }
    // END PRINTER NONSENSE

}


}



#endif // GSNUNF_BKPX2015_H
