//Needs optimizing currently testing.
#ifndef GSNUNF_BKPX2015_H
#define GSNUNF_BKPX2015_H

#include <array>
#include <iostream>
#include <fstream>
#include <list>
#include <cassert>
#include <string>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <limits>
#include <cmath>         // ceil, floor
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles
#include <boost/functional/hash.hpp> // size_t pair hash : used in Yao_inf_4
#include <CGAL/algorithm.h> //
#include <CGAL/circulator.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>

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


    enum AnchorType { None, Weak, Strong, StrongSelected, WeakSelected, StartOddChain };


    // project objects
    typedef pair<size_t,size_t>                                         size_tPair;
    typedef vector<pair<Vertex_handle, Vertex_handle>>                  fanCones;
    typedef boost::hash<size_tPair>                                     size_tPairHash;
    typedef unordered_map<size_tPair,bool,size_tPairHash>               size_tPairMap;
    typedef vector<pair<Vertex_handle, size_t>>                         yaoCones;
    typedef vector<size_t>                                              numYaoEdges;
    typedef vector<pair<Vertex_handle, AnchorType>>                     anchorCones;
    typedef pair<Vertex_handle, Vertex_handle>                          spannerEdge;
    typedef vector<vector<spannerEdge>>                                 spannerCones;

    // inline functions required for bkpx2015

    // get the cone of v wrt u (where u is at the center)
    inline size_t getCone( const Vertex_handle &u, const Vertex_handle &v) {
       bool X = (v->site().point().x()- u->site().point().x()) > 0;
       bool Y = (v->site().point().y()- u->site().point().y()) > 0;
       return !Y*2 + (!X*Y+X*!Y)*1;
    }

    // add yao edges
    inline void addYaoEdges(vector<yaoCones> &yaoEdges, vector<fanCones> &pointFans,
                            vector<numYaoEdges> &yaoEdgeCount,
                            const vector<Vertex_handle> &handles,
                            const SDG2 &sdg) {

        Vertex_circulator Circ = sdg.incident_vertices(handles[0]); // default value
        vector<double> distances(4);


        size_t cone = 0;
        size_t index = 0;

        for (auto point : handles) {

            yaoCones edges(4);
            fanCones fans(4);

            Circ = sdg.incident_vertices(point);
            distances = {INF, INF, INF, INF};
            index = point->storage_site().info();

            while( sdg.is_infinite(++Circ) );

            size_t previousCone = getCone(point, Circ);

            while( !sdg.is_infinite(++Circ) && getCone(point, Circ) == previousCone );

            if( sdg.is_infinite(Circ) )
                ++Circ;

            auto endpoint = Circ;

            fans[cone].first = Circ;

            do {

                if (!sdg.is_infinite(Circ)) {
                    cone = getCone(point,Circ);
                    ++(edges[cone].second);

                    if (cone != previousCone) {
                        fans[cone].first = Circ;
                    }
                    fans[cone].second = Circ;
                    previousCone = cone;

                    double proposedDistance = std::max(std::abs(point->site().point().x() - Circ->site().point().x()), std::abs(point->site().point().y() - Circ->site().point().y()));

                    if (proposedDistance < distances[cone]) {
                        distances[cone] = proposedDistance;
                        edges[cone].first = Circ;
                    }

                }

            } while(++Circ != endpoint); // finished determining the Yao edges + how many points are in fan of u's cone i


            yaoEdges[index] = edges;
            pointFans[index] = fans;

        } // finished moving through points

        // account for how many yao edges are incident upon each vertex u in cone i
        for (auto point : handles) {

            Circ = sdg.incident_vertices(point);
            numYaoEdges yaos(4);
            index = point->storage_site().info();

            cout << point->storage_site().info() << endl;

            for (size_t cone = 0; cone < 4; cone++) {

                if (yaoEdges[index][cone].second == 0) {
                    cout << "skip" << endl;
                    continue;
                }

                auto v1 = (pointFans[index][cone]).first;
                auto vk = (pointFans[index][cone]).second;

                cout << index << "] <" << v1->storage_site().info() << ", " << vk->storage_site().info() << "> " << endl;

                while (Circ != v1) {
                    ++Circ;
                }

                size_t position = 1;
                size_t total = yaoEdges[index][cone].second;

                while (position < (total+1)) {

                    cout << "test" << endl;

                    if (yaoEdges[index][cone].first == Circ || yaoEdges[Circ->storage_site().info()][(cone+2)%4].first == point) {
                        ++(yaos[cone]);
                    }

                    ++Circ;
                    ++position;

                }

            }

            cout << point->storage_site().info() << " is complete." << endl;

        yaoEdgeCount[index] = yaos;

        } // finished determining incident yao edges in cone i of vertex u

    } // yaoEdges are added


    inline void determineAnchors( vector<anchorCones> &anchorEdges,
                            vector<yaoCones> &yaoEdges, vector<fanCones> &pointFans,
                            vector<numYaoEdges> &yaoEdgeCount, vector<Vertex_handle> &handles,
                            SDG2 &sdg) {

        for (auto u : handles) {

            anchorCones anchors(4);

            for (size_t i = 0; i < 4; i++) {
                anchors[i].second = None;
            }

            size_t u_id = u->storage_site().info();
            anchorEdges[u_id] = anchors;

            for (size_t cone = 0; cone < 4; cone++) {

                bool anchorRetrieved = false;
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
                        anchorRetrieved = true;
                    }
                }

                if (yaoEdgeCount[u_id][cone] >= 2)
                {

                    auto v1 = pointFans[u_id][cone].first;
                    auto vk = pointFans[u_id][cone].second;

                    Vertex_circulator Circ = sdg.incident_vertices(u);
                    while (Circ != v1) {++Circ;}
                    size_t moves = 1;
                    size_t numNeighbors = yaoEdges[u_id][cone].second;
                    size_t lcount = 1;
                    size_t yaoTotal = 0;

                    while (moves < numNeighbors + 1)
                    {

                        if (yaoEdges[u_id][cone].first == Circ || yaoEdges[Circ->storage_site().info()][(cone+2)%4].first == u) {++yaoTotal;}

                        if (Circ == v)
                        {
                            lcount = moves;
                        }

                        ++moves;
                        ++Circ;

                    }

                    if (yaoTotal >= 2)
                    {
                        Vertex_circulator startCirc = sdg.incident_vertices(u);
                        while (startCirc != v) {++startCirc;}

                        // establish vlower and vhigher and recalibrate Circ accordingly
                        auto vlower = --startCirc;
                        ++startCirc;
                        auto vhigher = ++startCirc;
                        --startCirc;

                        bool cwCanonical = false;
                        bool ccwCanonical = false;

                        if (lcount >= 2 && (yaoEdges[vlower->storage_site().info()][getCone(vlower, v)].first == v &&
                            yaoEdges[v->storage_site().info()][getCone(v, vlower)].first != vlower)) {
                                ccwCanonical = true;
                            }

                        if (lcount <= (numNeighbors - 1) && (yaoEdges[vhigher->storage_site().info()][getCone(vhigher, v)].first == v) &&
                            yaoEdges[v->storage_site().info()][getCone(v, vhigher)].first != vhigher) {
                                cwCanonical = true;
                            }


                        if (ccwCanonical) {

                            bool inCanonical = true;
                            auto higher = startCirc;
                            size_t position = lcount;

                            while (inCanonical) {
                                higher = startCirc;
                                --startCirc;
                                --position;

                                if (position < 1) {
                                    anchorEdges[u_id][cone].first = higher;
                                    anchorEdges[u_id][cone].second = Weak;
                                    inCanonical = false;
                                    continue;
                                }

                                bool inCanonical = false;

                                if (yaoEdges[startCirc->storage_site().info()][getCone(startCirc, higher)].first == higher &&
                                    yaoEdges[higher->storage_site().info()][getCone(higher, startCirc)].first != startCirc) {
                                        inCanonical = true;
                                    }

                                if (!inCanonical) {
                                    anchorEdges[u_id][cone].first = higher;
                                    anchorEdges[u_id][cone].second = Weak;
                                }
                            }
                        }

                        else if (cwCanonical) {

                            bool inCanonical = true;
                            auto lower = startCirc;
                            size_t position = lcount;

                            while (inCanonical) {

                                lower = startCirc;
                                ++startCirc;
                                ++position;

                                if (position > numNeighbors) {

                                    anchorEdges[u_id][cone].first = lower;
                                    anchorEdges[u_id][cone].second = Weak;
                                    inCanonical = false;
                                    continue;

                                }

                                inCanonical = false;

                                if (yaoEdges[startCirc->storage_site().info()][getCone(startCirc, lower)].first == lower &&
                                    yaoEdges[lower->storage_site().info()][getCone(lower, startCirc)].first != startCirc) {
                                        inCanonical = true;
                                    }

                                if (!inCanonical) {
                                    anchorEdges[u_id][cone].first = lower;
                                    anchorEdges[u_id][cone].second = Weak;
                                }
                            }
                        }

                        else {
                            anchorEdges[u_id][cone].first = v;
                            anchorEdges[u_id][cone].second = Weak;
                        }
                    }
                } // when there are more than 2 edges!

            }

        }

        for (auto u : handles) {

            for (size_t cone = 0; cone < 4; cone++) {

                size_t u_id = u->storage_site().info();

                if (anchorEdges[u_id][cone].first == nullptr) {continue;}

                auto v = anchorEdges[u_id][cone].first;
                if (anchorEdges[v->storage_site().info()][(cone+2)%4].first == u || anchorEdges[v->storage_site().info()][(cone+2)%4].first == nullptr) {
                    anchorEdges[u_id][cone].second = Strong;
                }
            }

        } // strong anchors are identified


        // now it is time to select anchors
        for (auto u : handles) {

            size_t u_id = u->storage_site().info();

            for (size_t cone = 0; cone < 4; cone++) {

                if (anchorEdges[u_id][cone].first == nullptr) {continue;}

                if (anchorEdges[u_id][cone].second == Strong) {
                    anchorEdges[u_id][cone].second = StrongSelected;
                }

                if (anchorEdges[u_id][cone].second == Weak) {

                    auto current = u;
                    size_t current_id = current->storage_site().info();
                    size_t currentCone = cone;
                    bool inWeakAnchorChain = true;
                    bool previouslyAnalyzed = false;

                    // navigate to the end of the weak anchor chain (terminus = wk)
                    while (inWeakAnchorChain && !previouslyAnalyzed) {

                        if (anchorEdges[current_id][currentCone].second == WeakSelected) {
                            previouslyAnalyzed = true;
                        }

                        if (anchorEdges[current_id][currentCone].second == Strong || anchorEdges[current_id][currentCone].second == StrongSelected) {
                            inWeakAnchorChain = false;
                            continue;
                        }

                        current = anchorEdges[current_id][currentCone].first;
                        currentCone = (currentCone+2)%4;
                        current_id = current->storage_site().info();
                    }

                    inWeakAnchorChain = true;
                    size_t position = 0;

                    // select edges in the weak anchor chain
                    while (inWeakAnchorChain && !previouslyAnalyzed) {

                        Vertex_circulator Circ = sdg.incident_vertices(current);

                        auto v1 = pointFans[current_id][currentCone].first;

                        while (Circ != v1) {++Circ;}

                        size_t total = yaoEdges[current_id][currentCone].second;
                        size_t fanCount = 1;

                        size_t Circ_id = 0;

                        bool chainContinues = false;

                        while (fanCount < (total+1)) {

                           /* cout <<
                            cout << "Testing: ";
                            cout << Circ->storage_site().info() << endl;
                            cout << "Testing!" << endl;*/

                            Circ_id = Circ->storage_site().info();

                            if(anchorEdges[Circ_id][(currentCone+2)%4].first == current) {
                                fanCount = total+1;
                                chainContinues = true;
                                continue;
                            }

                            ++Circ;
                            ++fanCount;

                        }
                        if (chainContinues) {

                            cout << current_id << ", " << currentCone << "] " << endl;

                            current = Circ;
                            current_id = current->storage_site().info();

                            if (anchorEdges[current_id][(currentCone+2)%4].second == Strong) {
                                inWeakAnchorChain = false;
                            }
                            if (!inWeakAnchorChain && (position % 2 == 1)) {
                                auto startOdd = anchorEdges[current_id][(currentCone+2)%4].first;
                                anchorEdges[startOdd->storage_site().info()][currentCone].second = StartOddChain;
                            }

                            if (!inWeakAnchorChain) {continue;}

                            ++position;

                            if (position % 2 == 0) {
                                anchorEdges[current_id][(currentCone+2)%4].second = WeakSelected;
                            }

                            currentCone = (currentCone+2)%4;
                        }

                        else {

                            if (position % 2 == 1) {
                                anchorEdges[current_id][currentCone].second = StartOddChain;
                            }

                            cout << endl << "the chain's been processed at " << current_id << ", " << currentCone << endl;
                            inWeakAnchorChain = false;

                        }

                    }

                }

            }

        }

    } // function complete


  inline void degreeEightSpanner(vector<spannerCones> &H8,
                            vector<anchorCones> &anchorEdges,
                            vector<yaoCones> &yaoEdges, vector<fanCones> &pointFans,
                            vector<numYaoEdges> &yaoEdgeCount, vector<Vertex_handle> &handles,
                            SDG2 &sdg) {


    // put the edges conal vector
    for (auto w : handles) {

        spannerCones edges(4);
        H8[w->storage_site().info()] = edges;

    }

    for (auto w : handles) {

        size_t w_id = w->storage_site().info();
        Vertex_circulator Circ = sdg.incident_vertices(w);

        for (size_t cone = 0; cone < 4; cone++) {

            if (yaoEdges[w_id][cone].second == 0) {continue;}

            if (anchorEdges[w_id][cone].second == StrongSelected || anchorEdges[w_id][cone].second == WeakSelected) {

                auto v = anchorEdges[w_id][cone].first;
                bool original = true;

                for (auto path : H8[w_id][cone]) {
                    if (original && (w == path.first || w == path.second) && (v == path.first || v == path.second)) {original = false;}
                }

                if (original) {
                    spannerEdge edge = std::make_pair(w, v);
                    H8[w_id][cone].push_back(edge);
                    H8[v->storage_site().info()][(cone+2)%4].push_back(edge);
                }

            }

            if (yaoEdges[w_id][cone].second == 2) {

                auto v1 = pointFans[w_id][cone].first;
                auto vk = pointFans[w_id][cone].second;

                bool cw = false;
                bool ccw = false;

                if (yaoEdges[v1->storage_site().info()][getCone(v1, vk)].first == vk && yaoEdges[vk->storage_site().info()][getCone(vk, v1)].first != v1) {
                    ccw = true;
                }

                if (yaoEdges[v1->storage_site().info()][getCone(v1, vk)].first != vk && yaoEdges[vk->storage_site().info()][getCone(vk, v1)].first == v1) {
                    cw = true;
                }

                if (cw) {

                    // assess if both vertices are connected to u via yao edges

                    bool yaoConnection = false;

                    if ((yaoEdges[w_id][cone].first == v1 || yaoEdges[v1->storage_site().info()][(cone+2)%4].first == w)
                    && (yaoEdges[w_id][cone].first == vk || yaoEdges[vk->storage_site().info()][(cone+2)%4].first == w)) {
                        yaoConnection = true;
                    }

                    if (!yaoConnection) {continue;}

                    bool nonAnchor = true;

                    if (anchorEdges[v1->storage_site().info()][getCone(v1, vk)].first == vk || anchorEdges[vk->storage_site().info()][getCone(vk, v1)].first == v1) {
                        nonAnchor = false;
                    }

                    bool dual = false;

                    if (yaoEdgeCount[w_id][cone] > 1 && yaoEdgeCount[v1->storage_site().info()][(cone+2)%4] > 1) {
                        dual = true;

                        if (anchorEdges[v1->storage_site().info()][(cone+2)%4].first == w && anchorEdges[v1->storage_site().info()][(cone+2)%4].second == StartOddChain) {
                            dual = false;
                        }
                    }

                    if (nonAnchor && dual) {continue;}

             //       if (nonAnchor && yaoEdges[v1->storage_site().info()][(cone+2)%4].first != w) {continue;}

                    bool original = true;

                    for (auto path : H8[vk->storage_site().info()][getCone(vk, v1)]) {
                        if (original && (v1 == path.first || v1 == path.second) && (vk == path.first || vk == path.second)) {original = false;}
                    }

                    if (original) {


                        if (nonAnchor) {
                            spannerEdge edge = std::make_pair(vk, v1);
                            H8[vk->storage_site().info()][getCone(vk, v1)].push_back(edge);
                            H8[v1->storage_site().info()][getCone(v1, w)].push_back(edge);
                        }
                        else if (!nonAnchor){
                            spannerEdge edge = std::make_pair(vk, v1);
                            H8[vk->storage_site().info()][getCone(vk, v1)].push_back(edge);
                            H8[v1->storage_site().info()][getCone(v1, vk)].push_back(edge);
                        }
                    }

                }

                if (ccw) {

                    // assess if both vertices are connected to u via yao edges

                    bool yaoConnection = false;

                    if ((yaoEdges[w_id][cone].first == v1 || yaoEdges[v1->storage_site().info()][(cone+2)%4].first == w)
                    && (yaoEdges[w_id][cone].first == vk || yaoEdges[vk->storage_site().info()][(cone+2)%4].first == w)) {
                        yaoConnection = true;
                    }

                    if (!yaoConnection) {continue;}

                    bool nonAnchor = true;

                    if (anchorEdges[v1->storage_site().info()][getCone(v1, vk)].first == vk || anchorEdges[vk->storage_site().info()][getCone(vk, v1)].first == v1) {
                        nonAnchor = false;
                    }

                    bool dual = false;

                    if (yaoEdgeCount[w_id][cone] > 1 && yaoEdgeCount[vk->storage_site().info()][(cone+2)%4] > 1) {
                        dual = true;

                        if (anchorEdges[vk->storage_site().info()][(cone+2)%4].first == w && anchorEdges[vk->storage_site().info()][(cone+2)%4].second == StartOddChain) {
                            dual = false;
                        }
                    }

                    if (nonAnchor && dual) {continue;}

                    bool original = true;

                    for (auto path : H8[v1->storage_site().info()][getCone(v1, vk)]) {
                        if (original && (v1 == path.first || v1 == path.second) && (vk == path.first || vk == path.second)) {original = false;}
                    }

                    if (original) {

                        if (nonAnchor) {
                            spannerEdge edge = std::make_pair(v1, vk);
                            H8[v1->storage_site().info()][getCone(v1, vk)].push_back(edge);
                            H8[vk->storage_site().info()][getCone(vk, w)].push_back(edge);
                        }
                        else if (!nonAnchor){
                            spannerEdge edge = std::make_pair(v1, vk);
                            H8[v1->storage_site().info()][getCone(v1, vk)].push_back(edge);
                            H8[vk->storage_site().info()][getCone(vk, v1)].push_back(edge);
                        }
                    }

                }

            }




            if (yaoEdges[w_id][cone].second > 2) {

                auto v1 = pointFans[w_id][cone].first;
                auto vk = pointFans[w_id][cone].second;

                while (Circ != v1) {++Circ;}

                size_t total = yaoEdges[w_id][cone].second;
                size_t position = 1;

                while (position < total) {

                    if (position == 1) {

                        auto vnext = ++Circ;

                        // assess if both vertices are connected to u via yao edges

                        bool yaoConnection = false;

                        if ((yaoEdges[w_id][cone].first == v1 || yaoEdges[v1->storage_site().info()][(cone+2)%4].first == w)
                        && (yaoEdges[w_id][cone].first == vnext || yaoEdges[vnext->storage_site().info()][(cone+2)%4].first == w)) {
                            yaoConnection = true;
                        }

                        if (!yaoConnection) {continue;}


                        bool nonAnchor = true;

                        if (anchorEdges[v1->storage_site().info()][getCone(v1,vnext)].first == vnext ||
                            anchorEdges[vnext->storage_site().info()][getCone(vnext, v1)].first == v1) {

                                nonAnchor = false;

                            }

                        bool cw = false;
                        bool ccw = false;

                        if (yaoEdges[vnext->storage_site().info()][getCone(vnext,v1)].first == v1) {
                            cw = true;
                        }

                        if (yaoEdges[v1->storage_site().info()][getCone(v1, vnext)].first == vnext) {
                            ccw = true;
                        }

                        bool dual = false;

                        if (yaoEdgeCount[w_id][cone] > 1 && yaoEdgeCount[v1->storage_site().info()][(cone+2)%4] > 1) {
                            dual = true;

                            if (anchorEdges[v1->storage_site().info()][(cone+2)%4].first == w && anchorEdges[v1->storage_site().info()][(cone+2)%4].second == StartOddChain) {
                                dual = false;
                            }
                        }

                        if (nonAnchor && cw && !ccw && dual) {
                            ++position;
                            continue;
                        }

                        else {

                            // unidirectional going clockwise --> (v2,v1) in Yao
                            if (cw && !ccw) {

                             //   if (nonAnchor && yaoEdges[v1->storage_site().info()][(cone+2)%4].first != w) {continue;}

                                bool original = true;

                                for (auto path : H8[vnext->storage_site().info()][getCone(vnext, v1)]) {
                                    if (original && (v1 == path.first || v1 == path.second) && (vnext == path.first || vnext == path.second)) {original = false;}
                                }

                                if (original) {


                                    if (nonAnchor) {
                                        spannerEdge edge = std::make_pair(vnext, v1);
                                        H8[vnext->storage_site().info()][getCone(vnext, v1)].push_back(edge);
                                        H8[v1->storage_site().info()][getCone(v1, w)].push_back(edge);
                                    }
                                    else if (!nonAnchor) {
                                        spannerEdge edge = std::make_pair(vnext, v1);
                                        H8[vnext->storage_site().info()][getCone(vnext, v1)].push_back(edge);
                                        H8[v1->storage_site().info()][getCone(v1, vnext)].push_back(edge);
                                    }
                                }
                            }

                            if (!cw && ccw) {

                                if (nonAnchor && anchorEdges[w_id][cone].first == vnext) {
                                    ++position;
                                    continue;
                                }

                                bool original = true;

                                for (auto path : H8[v1->storage_site().info()][getCone(v1, vnext)]) {
                                    if (original && (v1 == path.first || v1 == path.second) && (vnext == path.first || vnext == path.second)) {original = false;}
                                }

                                if (original) {

                                    if (nonAnchor) {
                                        spannerEdge edge = std::make_pair(v1, vnext);
                                        H8[v1->storage_site().info()][getCone(v1, vnext)].push_back(edge);
                                        H8[vnext->storage_site().info()][getCone(vnext, w)].push_back(edge);
                                    }
                                    else if (!nonAnchor) {
                                        spannerEdge edge = std::make_pair(v1, vnext);
                                        H8[v1->storage_site().info()][getCone(v1, vnext)].push_back(edge);
                                        H8[vnext->storage_site().info()][getCone(vnext, v1)].push_back(edge);
                                    }
                                }
                            }
                        }

                        ++position;
                        continue;

                    }

                    if (position > 1 && position < (total-1)) {

                        auto vcurrent = Circ;
                        auto vnext = ++Circ;

                        // assess if both vertices are connected to u via yao edges

                        bool yaoConnection = false;

                        if ((yaoEdges[w_id][cone].first == vcurrent || yaoEdges[vcurrent->storage_site().info()][(cone+2)%4].first == w)
                        && (yaoEdges[w_id][cone].first == vnext || yaoEdges[vnext->storage_site().info()][(cone+2)%4].first == w)) {
                            yaoConnection = true;
                        }

                        if (!yaoConnection) {continue;}

                        bool nonAnchor = true;

                        if (anchorEdges[vcurrent->storage_site().info()][getCone(vcurrent, vnext)].first == vnext
                        || anchorEdges[vnext->storage_site().info()][getCone(vnext, vcurrent)].first == vcurrent) {
                            nonAnchor = false;
                        }

                        bool cw = false;
                        bool ccw = false;

                        if (yaoEdges[vnext->storage_site().info()][getCone(vnext, vcurrent)].first == vcurrent) {
                            cw = true;
                        }

                        if (yaoEdges[vcurrent->storage_site().info()][getCone(vcurrent, vnext)].first == vnext) {
                            ccw = true;
                        }



                        if (ccw && !cw) {

                            if (nonAnchor && anchorEdges[w_id][cone].first == vnext) {
                                ++position;
                                continue;
                            }

                            bool original = true;

                            for (auto path : H8[vcurrent->storage_site().info()][getCone(vcurrent, vnext)]) {
                                if (original && (vnext == path.first || vnext == path.second) && (vcurrent == path.first || vcurrent == path.second)) {original = false;}
                            }

                            if (original) {
                                spannerEdge edge = std::make_pair(vcurrent, vnext);
                                H8[vcurrent->storage_site().info()][getCone(vcurrent, vnext)].push_back(edge);

                                if (nonAnchor) {
                                    H8[vnext->storage_site().info()][getCone(vnext, w)].push_back(edge);
                                }
                                else {
                                    H8[vnext->storage_site().info()][getCone(vnext,vcurrent)].push_back(edge);
                                }
                            }
                        }

                        if (!ccw && cw) {

                            if (nonAnchor && anchorEdges[w_id][cone].first == vcurrent) {
                                ++position;
                                continue;
                            }

                            bool original = true;

                            for (auto path : H8[vnext->storage_site().info()][getCone(vnext, vcurrent)]) {
                                if (original && (vnext == path.first || vnext == path.second) && (vcurrent == path.first || vcurrent == path.second)) {original = false;}
                            }

                            if (original) {
                                spannerEdge edge = std::make_pair(vnext, vcurrent);
                                H8[vnext->storage_site().info()][getCone(vnext, vcurrent)].push_back(edge);

                                if (nonAnchor) {
                                    H8[vcurrent->storage_site().info()][getCone(vcurrent, w)].push_back(edge);
                                }
                                else {
                                    H8[vcurrent->storage_site().info()][getCone(vcurrent, vnext)].push_back(edge);
                                }
                            }
                        }

                        ++position;
                        continue;

                    }

                    if (position == (total-1)) {

                        auto vcurrent = Circ;
                        ++Circ;

                        // assess if both vertices are connected to u via yao edges

                        bool yaoConnection = false;

                        if ((yaoEdges[w_id][cone].first == vcurrent || yaoEdges[vcurrent->storage_site().info()][(cone+2)%4].first == w)
                        && (yaoEdges[w_id][cone].first == vk || yaoEdges[vk->storage_site().info()][(cone+2)%4].first == w)) {
                            yaoConnection = true;
                        }

                        bool nonAnchor = true;

                        if (anchorEdges[vcurrent->storage_site().info()][getCone(vcurrent, vk)].first == vk ||
                            anchorEdges[vk->storage_site().info()][getCone(vk, vcurrent)].first == vcurrent) {

                                nonAnchor = false;

                            }

                        bool cw = false;
                        bool ccw = false;

                        if (yaoEdges[vk->storage_site().info()][getCone(vk, vcurrent)].first == vcurrent) {
                            cw = true;
                        }

                        if (yaoEdges[vcurrent->storage_site().info()][getCone(vcurrent, vk)].first == vk) {
                            ccw = true;
                        }

                        bool dual = false;

                        if (yaoEdgeCount[vk->storage_site().info()][(cone+2)%4] > 1 && yaoEdgeCount[w->storage_site().info()][cone] > 1) {

                            dual = true;

                            if (anchorEdges[vk->storage_site().info()][(cone+2)%4].first == w && anchorEdges[vk->storage_site().info()][(cone+2)%4].second == StartOddChain) {
                                dual = false;
                            }

                        }

                        if (nonAnchor && !cw && ccw && dual) {
                            ++position;
                            continue;
                        }

                        else {

                            if (cw && !ccw) {

                                if (nonAnchor && anchorEdges[w_id][cone].first == vcurrent) {
                                    ++position;
                                    continue;
                                }

                                bool original = true;

                                for (auto path : H8[vk->storage_site().info()][getCone(vk, vcurrent)]) {
                                    if (original && (vk == path.first || vk == path.second) && (vcurrent == path.first || vcurrent == path.second)) {original = false;}
                                }

                                if (original) {

                                    if (nonAnchor) {
                                        spannerEdge edge = std::make_pair(vk, vcurrent);
                                        H8[vk->storage_site().info()][getCone(vk, vcurrent)].push_back(edge);
                                        H8[vcurrent->storage_site().info()][getCone(vcurrent, w)].push_back(edge);
                                    }
                                    else if (!nonAnchor) {
                                        spannerEdge edge = std::make_pair(vk, vcurrent);
                                        H8[vk->storage_site().info()][getCone(vk, vcurrent)].push_back(edge);
                                        H8[vcurrent->storage_site().info()][getCone(vcurrent, vk)].push_back(edge);
                                    }
                                }
                            }

                            if (!cw && ccw) {

                                if (nonAnchor) {continue;}

                                bool original = true;

                                for (auto path : H8[vcurrent->storage_site().info()][getCone(vcurrent, vk)]) {
                                    if (original && (vk == path.first || vk == path.second) && (vcurrent == path.first || vcurrent == path.second)) {original = false;}
                                }

                                if (original) {

                                    if (nonAnchor) {
                                        spannerEdge edge = std::make_pair(vcurrent, vk);
                                        H8[vcurrent->storage_site().info()][getCone(vcurrent, vk)].push_back(edge);
                                        H8[vk->storage_site().info()][getCone(vk, w)].push_back(edge);
                                    }
                                    else if (!nonAnchor) {
                                        spannerEdge edge = std::make_pair(vcurrent, vk);
                                        H8[vcurrent->storage_site().info()][getCone(vcurrent, vk)].push_back(edge);
                                        H8[vk->storage_site().info()][getCone(vk, vcurrent)].push_back(edge);
                                    }
                                }
                            }
                        }

                        ++position;
                        continue;

                    }

                }

            }


            }

        }

    }

    inline void processSpanner(vector<spannerCones> &H8,
                            vector<anchorCones> &anchorEdges,
                            vector<yaoCones> &yaoEdges, vector<fanCones> &pointFans,
                            vector<numYaoEdges> &yaoEdgeCount, vector<Vertex_handle> &handles,
                            SDG2 &sdg) {

        size_t u_id = handles[0]->storage_site().info();
        size_t charge = 0;

        for (auto u : handles) {

            u_id = u->storage_site().info();
            Vertex_circulator Circ = sdg.incident_vertices(u);

            for (size_t cone = 0; cone < 4; cone++) {

                charge = H8[u_id][cone].size();

                if (charge == 2) {

                    auto vnext = yaoEdges[u_id][cone].first;

                    bool plusFirst = false;
                    bool minusFirst = false;
                    bool endOfChain = true;

                    // determine if you are in the middle of the duplicate edge chain
                    if (H8[vnext->storage_site().info()][(cone+1)%4].size() == 2 && yaoEdges[vnext->storage_site().info()][(cone+2)%4].first != u
                        && anchorEdges[u_id][cone].first != vnext && anchorEdges[vnext->storage_site().info()][(cone+2)%4].first != u) {

                        plusFirst = true;
                        endOfChain = false;

                    }



                    if (!plusFirst && H8[vnext->storage_site().info()][(cone-1)%4].size() == 2 && yaoEdges[vnext->storage_site().info()][(cone+2)%4].first != u
                        && anchorEdges[u_id][cone].first != vnext && anchorEdges[vnext->storage_site().info()][(cone+2)%4].first != u) {

                        minusFirst = true;
                        endOfChain = false;

                    }

                    bool inDuplicateChain = true;
                    size_t position = 1;
                    auto vprev = vnext;

                    // cannot progress further through the chain --> figure out which orientation to go backwards with
                    if (endOfChain) {
                        inDuplicateChain = false;

                        // test out plusFirst
                        auto v1 = pointFans[u_id][(cone+1)%4].first;

                        Circ = sdg.incident_vertices(u);
                        while (Circ != v1) {++Circ;}

                        size_t fanTotal = yaoEdges[u_id][(cone+1)%4].second;
                        size_t fanCount = 1;

                        while (fanCount < (fanTotal+1)) {

                            if (!plusFirst && yaoEdges[Circ->storage_site().info()][(cone-1)%4].first == u && yaoEdges[u_id][(cone+1)%4].first != Circ &&
                                anchorEdges[Circ->storage_site().info()][(cone-1)%4].first != u && anchorEdges[u_id][(cone+1)%4].first != Circ) {
                                    plusFirst = true;
                                    vprev = Circ;
                                }

                            ++fanCount;
                            ++Circ;

                        }

                        if (!plusFirst) {

                                v1 = pointFans[u_id][(cone-1)%4].first;
                                while (Circ != v1) {++Circ;}

                                fanTotal = yaoEdges[u_id][(cone-1)%4].second;
                                fanCount = 1;

                                while (fanCount < (fanTotal+1)) {

                                    if (!minusFirst && yaoEdges[Circ->storage_site().info()][(cone+1)%4].first == u && yaoEdges[u_id][(cone-1)%4].first != Circ &&
                                        anchorEdges[Circ->storage_site().info()][(cone+1%4)].first != u && anchorEdges[u_id][(cone-1)%4].first != Circ) {
                                            minusFirst = true;
                                            vprev = Circ;
                                        }

                                    ++fanCount;
                                    ++Circ;

                                }

                        }

                    }

                    // moving forward through the chain
                    while (inDuplicateChain) {

                        if (position % 2 == 1) {

                            if (plusFirst) {

                                vnext = yaoEdges[vnext->storage_site().info()][(cone+1)%4].first;

                                if (H8[vnext->storage_site().info()][cone].size() != 2) {
                                    inDuplicateChain = false;
                                    continue;
                                }

                                ++position;

                            }

                            else {

                                vnext = yaoEdges[vnext->storage_site().info()][(cone-1)%4].first;

                                if (H8[vnext->storage_site().info()][cone].size() != 2) {
                                    inDuplicateChain = false;
                                    continue;
                                }

                                ++position;

                            }

                        }

                        if (position % 2 == 0) {

                            if (plusFirst) {

                                vnext = yaoEdges[vnext->storage_site().info()][cone].first;

                                if (H8[vnext->storage_site().info()][(cone+1)%4].size() != 2) {
                                    inDuplicateChain = false;
                                    continue;
                                }

                                ++position;

                            }

                            else {

                                vnext = yaoEdges[vnext->storage_site().info()][cone].first;

                                if (H8[vnext->storage_site().info()][(cone-1)%4].size() != 2) {
                                    inDuplicateChain = false;
                                    continue;
                                }

                                ++position;

                            }

                        }

                    }

                    // the end of the chain has been reached --> vnext = wk+1
                    inDuplicateChain = true;

                    vprev = vnext;
                    size_t chainCount = 0;
                    size_t fanCount = 1;
                    size_t fanTotal = 0;
                    bool previousFound = false;
                    bool startFound = false;
                    size_t edgeCone = 0;
                    Vertex_circulator Circ = sdg.incident_vertices(vnext);

                    cout << endl << endl;
                    cout << u_id << "] " << cone << endl;
                    cout << vprev->storage_site().info() << endl;
                    if (minusFirst) {cout << "minusFirst" << endl;}
                    if (plusFirst) {cout << "plusFirst" << endl;}
                    cout << "position: " << position << endl << endl;

                    while (inDuplicateChain) {

                        vnext = vprev;
                        Circ = sdg.incident_vertices(vnext);

                        cout << "position: " << position << endl;

                        cout << vnext->storage_site().info() << " is vnext!" << endl;

                        previousFound = false;

                        if (position % 2 == 1) {

                            if (plusFirst) {

                                auto v1 = pointFans[vnext->storage_site().info()][(cone+2)%4].first;

                                while (Circ != v1) {++Circ;}

                                fanCount = 1;
                                fanTotal = yaoEdges[vnext->storage_site().info()][(cone+2)%4].second;

                                while (fanCount < (fanTotal + 1)) {

                                    if (!previousFound && yaoEdges[Circ->storage_site().info()][cone].first == vnext && yaoEdges[vnext->storage_site().info()][(cone+2)%4].first != Circ &&
                                        anchorEdges[Circ->storage_site().info()][cone].first != vnext && anchorEdges[vnext->storage_site().info()][(cone+2)%4].first != Circ) {


                                        for (auto edge : H8[vnext->storage_site().info()][(cone-1)%4]) {

                                            if ((Circ == edge.first || Circ == edge.second) && (vnext == edge.first || vnext == edge.second)) {

                                                previousFound = true;
                                                vprev = Circ;
                                                fanCount = fanTotal+1;

                                                if (H8[vprev->storage_site().info()][cone].size() != 2) {startFound = true;}

                                                cout << vprev->storage_site().info() << " is previous!" << endl << endl;

                                            }

                                        }

                                    }

                                    ++Circ;
                                    ++fanCount;

                                }



                            }

                            else {

                                auto v1 = pointFans[vnext->storage_site().info()][(cone+2)%4].first;

                                while (Circ != v1) {++Circ;}

                                fanCount = 1;
                                fanTotal = yaoEdges[vnext->storage_site().info()][(cone+2)%4].second;

                                while (fanCount < (fanTotal + 1)) {

                                    if (!previousFound && yaoEdges[Circ->storage_site().info()][cone].first == vnext && yaoEdges[vnext->storage_site().info()][(cone+2)%4].first != Circ
                                     && anchorEdges[Circ->storage_site().info()][cone].first != vnext && anchorEdges[vnext->storage_site().info()][(cone+2)%4].first != Circ) {

                                        for (auto edge : H8[vnext->storage_site().info()][(cone+1)%4]) {

                                            if ((Circ == edge.first || Circ == edge.second) && (vnext == edge.first || vnext == edge.second)) {

                                                previousFound = true;
                                                vprev = Circ;
                                                fanCount = fanTotal+1;

                                                if (H8[vprev->storage_site().info()][cone].size() != 2) {startFound = true;}

                                                cout << vprev->storage_site().info() << " is previous!" << endl << endl;

                                            }

                                        }

                                    }

                                    ++Circ;
                                    ++fanCount;

                                }

                            }

                        }

                        if (position % 2 == 0) {

                            if (plusFirst) {

                                auto v1 = pointFans[vnext->storage_site().info()][(cone+1)%4].first;

                                while (Circ != v1) {++Circ;}

                                fanCount = 1;
                                fanTotal = yaoEdges[vnext->storage_site().info()][(cone+1)%4].second;

                                while (fanCount < (fanTotal+1)) {

                                    if (!previousFound && yaoEdges[Circ->storage_site().info()][(cone-1)%4].first == vnext && yaoEdges[vnext->storage_site().info()][(cone+1)%4].first != Circ
                                     && anchorEdges[Circ->storage_site().info()][(cone-1)%4].first != vnext && anchorEdges[vnext->storage_site().info()][(cone+1)%4].first != Circ) {

                                        for (auto edge : H8[vnext->storage_site().info()][cone]) {

                                            if ((Circ == edge.first || Circ == edge.second) && (vnext == edge.first || vnext == edge.second)) {

                                                previousFound = true;
                                                vprev = Circ;
                                                fanCount = fanTotal+1;

                                                if (H8[vprev->storage_site().info()][(cone-1)%4].size() != 2) {startFound = true;}

                                                cout << vprev->storage_site().info() << " is previous!" << endl << endl;

                                            }

                                        }

                                    }

                                    ++Circ;
                                    ++fanCount;

                                }


                            }

                            else {

                                auto v1 = pointFans[vnext->storage_site().info()][(cone-1)%4].first;

                                while (Circ != v1) {++Circ;}

                                fanCount = 1;
                                fanTotal = yaoEdges[vnext->storage_site().info()][(cone-1)%4].second;

                                while (fanCount < (fanTotal+1)) {

                                    if (!previousFound && yaoEdges[Circ->storage_site().info()][(cone+1)%4].first == vnext && yaoEdges[vnext->storage_site().info()][(cone-1)%4].first != Circ
                                     && anchorEdges[Circ->storage_site().info()][(cone+1)%4].first != vnext && anchorEdges[vnext->storage_site().info()][(cone-1)%4].first != Circ) {

                                        for (auto edge : H8[vnext->storage_site().info()][cone]) {

                                            if ((Circ == edge.first || Circ == edge.second) && (vnext == edge.first || vnext == edge.second)) {

                                                previousFound = true;
                                                vprev = Circ;
                                                fanCount = fanTotal+1;

                                                if (H8[vprev->storage_site().info()][(cone+1)%4].size() != 2) {startFound = true;}

                                                cout << vprev->storage_site().info() << " is previous!" << endl << endl;

                                            }

                                        }

                                    }

                                    ++Circ;
                                    ++fanCount;

                                }

                            }

                        }

                        ++position;


                        if (previousFound) {

                            ++chainCount;

                            cout << endl << "----------------------------" << endl << endl;

                            cout << "chainCount: " << chainCount << endl;

                            if (chainCount % 2 == 0) {

                                size_t edgeCone = getCone(vprev, vnext);
                                bool removed = false;
                                spannerEdge currentEdge = std::make_pair(vprev, vnext);

                                cout << "currentEdge: <" << currentEdge.first->storage_site().info() << ", " << currentEdge.second->storage_site().info() << "> " << endl;
                                cout << "edgeCone: " << edgeCone << endl;
                                size_t i = 0;
                                // remove the edge from vprev
                                for (auto edge : H8[vprev->storage_site().info()][edgeCone]) {

                                    if (!removed && (currentEdge.first == edge.first || currentEdge.first == edge.second) && (currentEdge.second == edge.first || currentEdge.second == edge.second)) {
                                        H8[vprev->storage_site().info()][edgeCone].erase(H8[vprev->storage_site().info()][edgeCone].begin() + i);
                                        removed = true;
                                    }
                                    ++i;

                                }

                                i = 0;
                                removed = false;
                                // remove the edge from vnext

                                if (plusFirst) {

                                    for (auto edge : H8[vnext->storage_site().info()][(edgeCone+1)%4]) {

                                    if (!removed && (currentEdge.first == edge.first || currentEdge.first == edge.second) && (currentEdge.second == edge.first || currentEdge.second == edge.second)) {

                                        cout << "plusFirst] DOUBLE EDGE FOUND! --> (" << edge.first->storage_site().info() << ", " << edge.second->storage_site().info() << "> " << endl << endl;

                                        H8[vnext->storage_site().info()][(edgeCone+1)%4].erase(H8[vnext->storage_site().info()][(edgeCone+1)%4].begin()+i);
                                        removed = true;
                                    }
                                    ++i;

                                    }

                                }

                                else {

                                    for (auto edge : H8[vnext->storage_site().info()][(edgeCone-1)%4]) {

                                    if (!removed && (currentEdge.first == edge.first || currentEdge.first == edge.second) && (currentEdge.second == edge.first || currentEdge.second == edge.second)) {

                                        cout << "minusFirst] DOUBLE EDGE FOUND! --> (" << edge.first->storage_site().info() << ", " << edge.second->storage_site().info() << "> " << endl << endl;

                                        H8[vnext->storage_site().info()][(edgeCone-1)%4].erase(H8[vnext->storage_site().info()][(edgeCone-1)%4].begin()+i);
                                        removed = true;;
                                    }
                                    ++i;

                                    }

                                }



                            }

                            if (startFound) {inDuplicateChain = false;}

                        }

                        else {
                            inDuplicateChain = false;
                        }
                    }
                }
            }
        }
    } // end namespace BKPX 2015
}

template<typename RandomAccessIterator, typename OutputIterator>
void BKPX2015(RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false) {

    using namespace bkpx2015;

    // construct Linf Delaunay triangulation

  /*
    vector<Point_2> P =
    {
        {0,0},
        {2,1},
        {3,3},
        {1,4},
        {-1.5,2}
    };


    SDG2 sdg;
    Site_2 site = Site_2::construct_site_2(P.front());

    for( id_type id=0; id<P.size(); ++id )
    {
        site = Site_2::construct_site_2(P[id]);
        sdg.insert(site,id);
        cout << id << "] ";
        cout << site.point() << "\n";
    }
    cout << endl;

    assert( sdg.is_valid(true, 1) );

    cout<<endl<<endl;

    // Print the points and IDs
    for( auto it = sdg.finite_vertices_begin();
         it != sdg.finite_vertices_end(); ++it )
    {
        cout << it->storage_site().info() << "] ";
        cout << it->site().point() << endl;
    }
    cout << endl;
  */

   // ifstream ifs("data6.txt");
   // assert( ifs );

    vector<Point_2> P(pointsBegin, pointsEnd);

    SDG2 sdg;
    Site_2 site;
    id_type id = 0;

    for (Point_2 point : P) {
      site = Site_2::construct_site_2(point);
      sdg.insert(site, id);
      ++id;
    }

    assert( sdg.is_valid(true, 1) );
    cout << endl << endl;

    //  Print the points and IDs
    for( auto it = sdg.finite_vertices_begin();
         it != sdg.finite_vertices_end(); ++it )
    {
        cout << it->storage_site().info() << "] ";
        cout << it->site().point() << endl;
    }
    cout << endl;

    //N is the number of vertices in the delaunay triangulation.
    size_t n = sdg.number_of_vertices();
    if(n > SIZE_T_MAX - 1) return;

    // store the vertex handles
    vector<Vertex_handle> handles(n);

    for (auto v = sdg.finite_vertices_begin(); v != sdg.finite_vertices_end(); v++) {
        handles[v->storage_site().info()] = v;
    }

    //construct YaoEdges
    vector<yaoCones> yaoEdges(n);
    vector<fanCones> pointFans(n);
    vector<numYaoEdges> yaoEdgeCount(n);

    addYaoEdges(yaoEdges, pointFans, yaoEdgeCount, handles, sdg);

    cout << "yao edge terminals directed from u" << endl;

    for (auto u : handles) {

        size_t u_id = u->storage_site().info();

        cout << u_id << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

            if (yaoEdges[u_id][cone].second == 0) {
                cout << "x ";
            }
            else {
              //  size_t neighbor_id = yaoEdges[u_id][cone].first
                cout << (yaoEdges[u_id][cone].first)->storage_site().info() << " ";
            }

            if (cone == 3) {
                cout << endl;
            }

        }

    }

    cout << endl << endl;

    cout << "number of edges in cones" << endl;

    for (auto u : handles) {

        size_t u_id = u->storage_site().info();

        cout << u_id << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

              // num edges in that cone
            cout << yaoEdges[u_id][cone].second << " ";


            if (cone == 3) {
                cout << endl;
            }

        }

    }

    cout << endl << endl;

    cout << "number of YAO EDGES in cones" << endl;

    for (auto u : handles) {

        size_t u_id = u->storage_site().info();

        cout << u_id << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

              // num edges in that cone
            cout << yaoEdgeCount[u_id][cone] << " ";


            if (cone == 3) {
                cout << endl;
            }

        }

    }

    cout << endl << endl;

    cout << "fans of the cones" << endl;

    for (auto u : handles) {

        size_t u_id = u->storage_site().info();

        cout << u_id << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

            if (yaoEdges[u_id][cone].second == 0) {
                cout << "< x, x > ";
            }
            else {
                cout << "< " << pointFans[u_id][cone].first->storage_site().info() << ", " << pointFans[u_id][cone].second->storage_site().info() << " > ";
            }

            if (cone == 3) {
                cout << endl;
            }

        }

    }

    // identify the anchors
    vector<anchorCones> anchorEdges(n);

    determineAnchors(anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, sdg);


    cout << endl << "anchors" << endl;

    for (auto u : handles) {

        size_t u_id = u->storage_site().info();

        cout << u_id << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

            if (anchorEdges[u_id][cone].first == nullptr)
            {
                cout << "x ";
            }
            else
            {
                cout << (anchorEdges[u_id][cone].first)->storage_site().info() << " ";
            }
              // num edges in that cone

            if (cone == 3) {
                cout << endl;
            }

        }

    }

    cout << endl << endl << "strong anchors" << endl;

    for (auto u : handles) {

        size_t u_id = u->storage_site().info();

        cout << u_id << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

            if (anchorEdges[u_id][cone].second == Strong || anchorEdges[u_id][cone].second == StrongSelected)
            {
                cout << (anchorEdges[u_id][cone].first)->storage_site().info() << " ";
            }
            else
            {
                cout << "x ";
            }
              // num edges in that cone

            if (cone == 3) {
                cout << endl;
            }

        }

    }

    cout << endl << endl << "selected anchors" << endl;

    for (auto u : handles) {

        size_t u_id = u->storage_site().info();

        cout << u_id << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

            if (anchorEdges[u_id][cone].second == StrongSelected || anchorEdges[u_id][cone].second == WeakSelected)
            {
                cout << (anchorEdges[u_id][cone].first)->storage_site().info() << " ";
            }
            else
            {
                cout << "x ";
            }
              // num edges in that cone

            if (cone == 3) {
                cout << endl;
            }

        }

    }

    cout << endl << endl << "start of odd chain" << endl;

    for (auto u : handles) {

        size_t u_id = u->storage_site().info();

        cout << u_id << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

            if (anchorEdges[u_id][cone].second == StartOddChain)
            {
                cout << (anchorEdges[u_id][cone].first)->storage_site().info() << " ";
            }
            else
            {
                cout << "x ";
            }
              // num edges in that cone

            if (cone == 3) {
                cout << endl;
            }

        }

    }

    cout << endl << endl;

   // construct H8 --> degree 8 spanner
    vector<spannerCones> H8(n);

    degreeEightSpanner(H8, anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, sdg);

    cout << "degree 8 spanner edges" << endl;

    for (auto u : handles) {

        cout << u->storage_site().info() << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

            if (H8[u->storage_site().info()][cone].size() == 0) {
             //   cout << "(cone " << cone << " has no edge)";
            }
            else {
                for (auto edge : H8[u->storage_site().info()][cone]) {
                    cout << "(" << (edge.first)->storage_site().info() << "," << (edge.second)->storage_site().info() <<") ";
                }
            }
        }

        cout << endl;

    }

    cout << endl << endl;

    cout << "degree 8 spanner charges" << endl;

    for (auto u : handles) {

        cout << u->storage_site().info() << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

            cout << H8[u->storage_site().info()][cone].size() << " ";

        }

        cout << endl;

    }

    cout << endl << endl;


    // AS of here step 2 is complete. it likely has logical errors that will need to be assessed, but it builds so that's good :-)

  processSpanner(H8, anchorEdges, yaoEdges, pointFans, yaoEdgeCount, handles, sdg);

  cout << "degree 4 spanner edges" << endl;

  for (auto u : handles) {

        cout << u->storage_site().info() << "] ";

        for (size_t cone = 0; cone < 4; cone++) {

            if (H8[u->storage_site().info()][cone].size() == 0) {
             //   cout << "(cone " << cone << " has no edge)";
            }
            else {
                for (auto edge : H8[u->storage_site().info()][cone]) {
                    cout << "(" << (edge.first)->storage_site().info() << "," << (edge.second)->storage_site().info() <<") ";
                }
            }
        }

        cout << endl;

    }

}


}

#endif // GSNUNF_BKPX2015_H
