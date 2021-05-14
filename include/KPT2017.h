//Needs optimizing currently testing.
#ifndef GSNUNF_KPT2017_H
#define GSNUNF_KPT2017_H

//Base libraries.
#include <cmath>         // ceil, floor, isinf
#include <limits>
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

//Boost library
#include <boost/functional/hash.hpp> // size_t pair hash

//CGAL library
#include <CGAL/algorithm.h>
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Line_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

//Project library
#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "metrics.h"
#include "TDDelaunay.h"
#include "utilities.h"


namespace gsnunf {

    using namespace std;

    namespace kpt2017 {

    typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
    typedef Epick                     K;
    typedef K::Point_2                Point_2;
    typedef K::FT                     FT;

    typedef HalfThetaTriangulation<K> TD_Delaunay_2;
    typedef TD_Delaunay_2::Vertex_descriptor Vertex_descriptor;


    enum Color {
        BLUE, WHITE, RED = WHITE, GREEN = WHITE
    };

    //Cone angles.
    const double tan30 = TAN30;
    const FT cot30 = 1 / tan30;

    const FT alpha = PI/3;

    //Slopes of the cone boundary lines.
    const vector<double> bisectorSlopes{ INF, tan30, -1*tan30, INF, tan30, -1*tan30 };
    const vector<double> orthBisectorSlopes{ 0, -1*cot30, cot30, 0, -1*cot30, cot30 };

    //Finds the cone of p containing vertex q, for this algorithm all vertices have 6 cones (0-5) with an angle of (PI/3).
    inline size_t getSingleCone(const size_t p, const size_t q, const vector<Point_2> &h){
        const Point_2 refPoint( h.at(p).x() - tan30, h[p].y() + 1 );
        //Point refPoint(h[p]->point().x(), h[p] ->point().y() + 1);

        double theta = get_angle<K>(refPoint, h[p], h.at(q));

        size_t cone = (theta / alpha);

        return cone;
    }

    //Compute max of getCone(p,q) and (getCone(q,p)+3)%6, is used to make sure cones are calculated correctly.
    inline size_t getCone( const size_t p, const size_t q, const vector<Point_2> &h ) {
        if( h[p] < h[q] ) {
            return getSingleCone(p,q,h);
        } else {
            return (getSingleCone(q,p,h)+3)%6;
        }
    }

    //Finds the bisector length of a given edge.
    inline K::FT bisectorLength( const pair<size_t,size_t> &e, const vector<Point_2> &h ) {

        size_t cone = getCone(e.first, e.second, h);

        double xCord = h.at(e.first).x();
        double yCord = h[e.first].y() + 1;

        assert(cone<6);
        assert(e.first<h.size());

        xCord = h[e.first].x() - orthBisectorSlopes.at(cone);

        Point_2 bisectorPoint(xCord, yCord);

        K::Line_2 bisectorLine(h[e.first], bisectorPoint);

        Point_2 intersectionPoint = bisectorLine.projection(h[e.second]);

        double bisectorLen = distance(h[e.first], intersectionPoint);

        return bisectorLen;
    }
    struct bisectorLengthComp {
        bool operator()( const auto &lhs, const auto &rhs ) {
            return lhs < rhs;
        }
    };
//
//    /*
//      Step 3: Add incident edges consists of 2 sub steps.
//      (3.1) Starts with the empty set E_A.
//      (3.2) For each edge in L (sorted in non-decreasing order) Let i be the cone of p containing q if E_A has no edges with endpoint p in the
//            neighborhood of p in cone i and E_A has no edges with endpoint q in the neighborhood of q in cone i+3 then add edge (p,q) to E_A.
//    */
//    inline void addIncident( vector<pair<size_t,size_t>> &E_A,
//                             pointConeMap &AL_E_A,
//                             const double alpha,
//                             const vector<Vertex_handle> &h,
//                             const vector<pair<pair<size_t,size_t>,double>> &l ) {
//        //Loops through the entire set L.
//        for( auto e : l ) {
//
//            //Separates the edge (p,q) into the vertices p and q.
//            const size_t& p = e.first.first;
//            const size_t& q = e.first.second;
//
//            //Computes the cone of p containing q.
//            size_t p_cone = getCone(p, q, h),
//                   q_cone = getCone(q, p, h);
//
//            /*Evaluates the emptiness of a cone, given a vertex in the set E_A.
//              If a cone is empty then the set E_A does not contain an edge with the
//              given endpoint in the cone calculated above, and the status will be
//              set to true.*/
//            bool p_cone_empty = AL_E_A.find(make_pair(p, p_cone)) == AL_E_A.end();
//            bool q_cone_empty = AL_E_A.find(make_pair(q, q_cone)) == AL_E_A.end();
//
//            /*Checks that both cone neighborhood are empty, if these are both empty
//            then the condition for step 3 is met and (p,q) is added to E_A. (3.2)*/
//            if( p_cone_empty && q_cone_empty ){
//                E_A.push_back(e.first);
//
//                //Adds (p,q) to an adjacency list for future calculation.
//                AL_E_A.emplace(make_pair(p, p_cone), q);
//                AL_E_A.emplace(make_pair(q, q_cone), p);
//            }
//        }
//    }
//
//    inline void canonicalNeighborhood( vector<size_t>& canNeighbors,
//                                       const size_t& p,
//                                       const size_t& r,
//                                       const size_t cone,
//                                       const Delaunay &dt,
//                                       const vector<Vertex_handle> &h,
//                                       const edgeBisectorMap &b,
//                                       bool printLog = false ) {
//
//        pair<size_t,size_t> e = make_pair(p, r);
//
//        /*Vertex circulator oriented to r to find fist and last end vertex. Once r is the circulator is oriented to the first vertex in the cone,
//          that is in the canonical neighborhood. Once found all neighbors are added in clockwise order. For a vertex to be in the canonical
//          neighborhood it must have a bisector length greater than or equal to that of (p,r)
//        */
//        auto N_p = dt.incident_vertices(h[p]);
//
//        while(++N_p != h[r]);
//
//        while(!dt.is_infinite(++N_p) && getCone(p, N_p->info(), h) == cone && (b.at(make_pair(p, N_p->info())) > b.at(e)
//                || abs(b.at(make_pair(p, N_p->info())) - b.at(e)) < EPSILON));
//
//
//        while(!dt.is_infinite(--N_p) && getCone(p, N_p->info(), h) == cone && (b.at(make_pair(p, N_p->info())) > b.at(e)
//                  || abs(b.at(make_pair(p,N_p->info())) - b.at(e)) < EPSILON)){
//            canNeighbors.push_back(N_p->info());
//        }
//    }
//
//    /*
//      Step 4: Add canonical edges consists of 4 sub steps.
//      (4.1) Let r be an element of the canonical neighborhood of p in the cone 0 of p.
//      (4.2) Add inner edges if total neigborhood edges is 3 or more.
//      (4.3) If r is an end vertex and there is more than one edge in the neighborhood add the edge with endpoint r.
//      (4.4) Consider the first and last edge in the canoncial neighborhood. 3 criteria to add.
//        (4.4 a) If the edges are in cone 1 or 5 with respect to a and z add.
//        (4.4 b) If the edges are in cone 2 or 4 with respect to a and z and cone for has no edge with an end edge point in E_A add.
//        (4.4 c) Checks if end edges have a end point a or z in E_A and an edge different from one made with vertex b or y in cone 2 or 4 woth respect
//                to a and z if found the edge (b,c) or (w,y) is added.
//    */
//    inline void addCanonical( vector<pair<size_t,size_t>> &E_CAN,
//                              const size_t p,
//                              const size_t r,
//                              const double alpha,
//                              const Delaunay &dt,
//                              const vector<Vertex_handle> &h,
//                              const edgeBisectorMap &b,
//                              pointConeMap& AL_e_a,
//                              bool printLog=false ) {
//
//        //Creates an edge (p,r)
//        pair<size_t,size_t> e = make_pair(p, r);
//
//        //Computes the cone of p containing r.
//        size_t p_cone  = getCone(p, r, h);
//
//        //Set of the canonical neighborhood of p in the cone of p containing r. (This cone will be considered as cone 0)
//        vector<size_t> canNeighbors;
//
//        canonicalNeighborhood( canNeighbors, p, r, p_cone, dt, h, b );
//        assert(canNeighbors.size()>0);
//
//        //Number of edges in the neighborhood.
//        int canEdges = canNeighbors.size() - 1;
//
//        //Must be at least 1 edge.
//        if(canEdges > 1){
//            //Add inner edges if total neighborhood edges is 3 or more. (4.2)
//            for(int i = 1; i < canEdges-1; i++){
//                E_CAN.emplace_back(canNeighbors.at(i), canNeighbors.at(i + 1));
//                assert(dt.is_edge( h.at(canNeighbors.at(i)), h.at(canNeighbors.at(i + 1))));
//            }
//
//            //End edges in the canonical neighborhood.
//            const vector<pair<size_t,size_t>> canExtrema {
//                make_pair( canNeighbors.at(1), canNeighbors.front() ),
//                make_pair( canNeighbors.at(canEdges - 1), canNeighbors.back() )
//            };
//
//            //If r is an end vertex and there is more than one edge in the neighborhood add the edge with endpoint r. (4.3
//            for( auto e : canExtrema ) {
//                if( e.second == r && canEdges > 1) {
//                    E_CAN.push_back(e);
//                    assert(dt.is_edge( h.at(e.first), h.at(e.second)));
//                }
//            }
//
//            //First and last edges in the canonical neighborhood are considered and added by 3 criteria. (4.4)
//            vector<int> cone(6);
//            for( size_t i=0; i<6; ++i ) {
//                cone[i] = (p_cone+i)%6;
//            }
//
//            //If the edges are in cone 1 or 5 with respect to a and z add. (4.4 a)
//            for( auto i=0; i<canExtrema.size(); ++i ) {
//                const auto e = canExtrema[i];
//                const int z_cone = 1 + int(i==1)*4;
//                if( getCone(e.second, e.first, h) == cone[z_cone] ) {
//                    E_CAN.push_back(e);
//                    if(printLog) cout<<e.first<<"-"<<e.second<<"["<<i<<"]\\"<<cone[z_cone]<<"/,";
//
//                    assert( dt.is_edge( h.at(e.second), h.at(e.first) ) );
//                }
//            }
//
//            const vector<pointConeMap::iterator> endpointZ {
//                AL_e_a.find( make_pair(canExtrema[0].second, cone[2]) ),
//                AL_e_a.find( make_pair(canExtrema[1].second, cone[4]) )
//            };
//            const auto blank = AL_e_a.end(); //Iterator to end of map to check if an edge exists.
//
//            //If the edges are in cone 2 or 4 with respect to a and z and cone for has no edge with an end edge point in E_A add. (4.4 b)
//            for( auto i=0; i<canExtrema.size(); ++i ) {
//                const auto e = canExtrema[i];
//                const int z_cone = 2 + int(i==1)*2;
//                if( endpointZ[i] == blank && getCone(e.second, e.first, h) == cone[z_cone] ){
//                    E_CAN.push_back(e);
//
//                    assert( dt.is_edge( h.at(e.first), h.at(e.second) ) );
//                }
//            }
//
//            /*Checks if end edges have an end point a or z in E_A and an edge different from one made with vertex b or y in cone 2 or 4 woth respect
//              to a and z if found the edge (b,c) or (w,y) is added. (4.4 c)*/
//            // (a,b)
//            for( auto i=0; i<canExtrema.size(); ++i ) {
//                const auto e = canExtrema[i];
//                const int z_cone = 2 + int(i==1)*2;
//                if( getCone(e.second, e.first, h) == cone[z_cone]
//                 && endpointZ[i] != blank
//                 && endpointZ[i]->second != e.first ) {
//                    vector<size_t> zCanNeighbors;
//                    canonicalNeighborhood(
//                        zCanNeighbors, e.second, endpointZ[i]->second,
//                        cone[z_cone], dt, h, b
//                    );
//                    auto y = find( zCanNeighbors.begin(), zCanNeighbors.end(), e.first );
//                    auto w = y;
//                    w += int(y == zCanNeighbors.begin());
//                    w -= int(y == zCanNeighbors.end()-1);
//
//                    assert( y != zCanNeighbors.end() );
//                    assert( y != w );
//                    assert(dt.is_edge( h.at(*w), h.at(*y)));
//                    E_CAN.emplace_back(*w, *y);
//                }
//            }
//        }
//    }
template< class AnchorListMap, class Triangulation, class PointContainer >
void findAnchors( AnchorListMap &anchors, Triangulation &D, const PointContainer &P )
{
    anchors =
    {
        { BLUE, {} },
        { WHITE, {} }
    };

    // find and classify anchors
    for( auto vit=D.finite_vertices_begin();
         vit!=D.finite_vertices_end(); ++vit )
    {
        auto w = *vit;
        cout<<w<<") "<<D.point(w)<<endl;

        Vertex_descriptor local_anchors[3];
        K::FT local_minimum_bisectors[3] = {INF,INF,INF};

        // get the edges for which w is the target
        for( auto eit=D.negative_cone_edges_begin(w);
             eit!=D.negative_cone_edges_end(w); ++eit )
        {
            auto e = *eit;
            auto u = D.target(e),
                 v = D.source(e);
            size_t cone = getCone( u, v, P ),
                   flattened_cone = cone/2;
            auto e_bisector_length = bisectorLength(make_pair(u,v),P);
            if( e_bisector_length < local_minimum_bisectors[flattened_cone]) {
                local_anchors[flattened_cone] = v;
                local_minimum_bisectors[flattened_cone] = e_bisector_length;
            }
            cout<<"  "<<D.source(e)<<" is in cone "<<cone<<endl;
        }

        // set the anchor for each cone
        for( size_t i=0; i<3; ++i ) {
            if( local_minimum_bisectors[i] < INF ) {
                Color c = static_cast<Color>(int(i!=0));
                anchors[c].emplace_back( local_anchors[i], w );
                cout<<"  added "<< (c==BLUE?"blue":"white")<<" anchor "<<local_anchors[i]<<"-"<< w<<endl;
            }
        }
    }
}

template< class AnchorList, class PointContainer, class EdgeList >
void addWhiteAnchors(AnchorList &whiteAnchors, const PointContainer &P, EdgeList &A )
{
    std::sort( whiteAnchors.begin(), whiteAnchors.end(),
        [&] ( const auto &lhs, const auto &rhs )
        {
            return bisectorLength(lhs,P) < bisectorLength(rhs,P);
        }
    );

    // Add all white anchors whose adjacent white cone is empty
    struct PositiveAndNegativeWhiteConesForVertex
    {    // each vertex can have zero or one of both positive and negative white anchors
        bool positive = true;
        bool negative = true;
    };
    unordered_map<Vertex_descriptor, PositiveAndNegativeWhiteConesForVertex> AdjacentWhiteConeIsEmpty;

    for( auto whiteAnchor : whiteAnchors )
    {
        cout<<"WhiteAnchor:"<<whiteAnchor.first<<"-"<<whiteAnchor.second<<"\n";

        Vertex_descriptor higherVertex = whiteAnchor.first,
                          lowerVertex = whiteAnchor.second;

        if( P[higherVertex].y() < P[lowerVertex].y() )
        {
            swap(higherVertex,lowerVertex);
        }

        bool higherVertexIsKnown = contains( AdjacentWhiteConeIsEmpty, higherVertex );
        bool higherVertexAgrees = // vertex with higher y value needs its negative (second) cones checked
        {
            !higherVertexIsKnown
            || AdjacentWhiteConeIsEmpty[higherVertex].negative
        };
        bool lowerVertexIsKnown = contains( AdjacentWhiteConeIsEmpty, lowerVertex );
        bool lowerVertexAgrees = // other vertex needs its positive (first) cones checked
        {
            !lowerVertexIsKnown
            || AdjacentWhiteConeIsEmpty[lowerVertex].positive
        };

        if( higherVertexAgrees && lowerVertexAgrees )
        {
            cout<<"Both agree\n";
            A.insert(whiteAnchor);

            AdjacentWhiteConeIsEmpty.try_emplace(higherVertex, PositiveAndNegativeWhiteConesForVertex() );
            AdjacentWhiteConeIsEmpty[higherVertex].negative = false;

            AdjacentWhiteConeIsEmpty.try_emplace(lowerVertex, PositiveAndNegativeWhiteConesForVertex() );
            AdjacentWhiteConeIsEmpty[lowerVertex].positive = false;
        }
    }
}
template< class Triangulation, class EdgeList, class PointContainer, class AnchorList >
void addCanonicalEdgesInBlueCones( const Triangulation &D, const EdgeList &A, const PointContainer &P, const AnchorList &blueAnchors, EdgeList &S )
{
    for( auto blueAnchor : blueAnchors )
    {
        //auto source = blueAnchor.first;
        typedef typename Triangulation::Edge_descriptor Edge_descriptor;
        vector<Edge_descriptor> fan;
        auto target = blueAnchor.second;
        D.fan_of_cone(target,1,fan);

        auto eit = fan.begin();
        auto previousVertex = D.source(*eit);
        // get fan of blueAnchor's target vertex cone 1 (negative/blue)
        for( ++eit; eit!=fan.end(); ++eit )
        {
            auto thisVertex = D.source(*eit);
            auto canonicalEdge = make_pair( thisVertex, previousVertex );
            if( !contains(A, canonicalEdge)
             && !contains(A, reverse_pair(canonicalEdge)) )
            {
                cout<<"Canonical edge "<<canonicalEdge.first<<"-"<<canonicalEdge.second<<"\n";
                S.insert(canonicalEdge);
            }
            previousVertex = thisVertex;
        }
    }
}


} // namespace KPT2017


// Main algorithm.
template<typename RandomAccessIterator, typename OutputIterator>
void KPT2017(RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false)
{
    using namespace kpt2017;

    //Angle of the cones. Results in 6 cones for a given vertex.

    vector<Point_2> P( pointsBegin, pointsEnd );

    TD_Delaunay_2 D( P.begin(), P.end() );

    map<Color,vector<pair<size_t,size_t>>> Anchors;
    findAnchors( Anchors, D, P );

    // Step 1. add all blue anchors to A
    typedef set<pair<size_t,size_t>> EdgeList;
    EdgeList A( Anchors[BLUE].begin(),
                Anchors[BLUE].end() );

    // Step 2.
    addWhiteAnchors( Anchors[WHITE], P, A );

    // Step 3.
    EdgeList S(A);
    // Add to S every canonical edge in negative blue cones (cone 1) if the edge isn't in A
    addCanonicalEdgesInBlueCones(D,A,P,Anchors[BLUE],S);


    // Edge list is only needed for printing. Remove for production.
    vector<pair<Point_2,Point_2>> edgeList;

    // Send resultant graph to output iterator
    for(auto e : S) {
        // Edge list is only needed for printing. Remove for production.
        edgeList.emplace_back(P.at(e.first), P.at(e.second));

        *result = make_pair(P.at(e.first), P.at(e.second));
        ++result;
        *result = make_pair(P.at(e.second), P.at(e.first));
        ++result;
    }

    // START PRINTER NONSENSE
    if(printLog) {
        GraphPrinter printer(5);
        GraphPrinter::OptionsList options;

        options = {
            {"color", printer.inactiveEdgeColor},
            {"line width", to_string(printer.inactiveEdgeWidth)}
        };
        printer.drawEdgesOfHalfTheta(D, options);

        options = { // active edge options
            {"color", printer.activeEdgeColor},
            {"line width", to_string(printer.activeEdgeWidth)}
        };
        printer.drawEdges(edgeList.begin(), edgeList.end(), options);


        options = {
            {"vertex", make_optional(to_string(printer.vertexRadius))}, // vertex width
            {"color", make_optional(printer.backgroundColor)}, // text color
            {"fill", make_optional(printer.activeVertexColor)}, // vertex color
            {"line width", make_optional(to_string(0))} // vertex border (same color as text)
        };
        GraphPrinter::OptionsList borderOptions = {
            {"border", make_optional(to_string(printer.vertexRadius))}, // choose shape of vertex
            {"color", printer.activeEdgeColor}, // additional border color
            {"line width", to_string(printer.inactiveEdgeWidth)}, // additional border width
        };
        printer.drawVerticesWithInfo(D.points_begin(), D.points_end(), options, borderOptions);

        printer.print("KPT2017");
        cout << "\n";
    }
    // END PRINTER NONSENSE

} // function KPT2017

} // namespace gsnunf

#endif // GSNUNF_KPT2017_H
