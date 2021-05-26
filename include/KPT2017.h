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
inline size_t getSingleCone(const size_t p, const size_t q, const vector<Point_2> &h)
{
    const Point_2 refPoint( h.at(p).x() - tan30, h[p].y() + 1 );
    //Point refPoint(h[p]->point().x(), h[p] ->point().y() + 1);

    double theta = get_angle<K>(refPoint, h[p], h.at(q));

    size_t cone = (theta / alpha);

    return cone;
}

//Compute max of getCone(p,q) and (getCone(q,p)+3)%6, is used to make sure cones are calculated correctly.
inline size_t getCone( const size_t p, const size_t q, const vector<Point_2> &h )
{
    return p < q ?
        getSingleCone(p,q,h)
        : ( getSingleCone(q,p,h)+3 ) % 6;
}

//Compute max of getCone(p,q) and (getCone(q,p)+3)%6, is used to make sure cones are calculated correctly.
inline Color getColor( const size_t p, const size_t q, const vector<Point_2> &h )
{
    cout<<"Getting color of "<<p<<"-"<<q<<endl;
    return getCone(p,q,h) % 3 == 1 ? BLUE : WHITE;
}

//Finds the bisector length of a given edge.
inline K::FT bisectorLength( const pair<size_t,size_t> &e, const vector<Point_2> &h )
{
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

        // negative cones
        Vertex_descriptor local_anchors[3];
        K::FT local_minimum_bisectors[3] = {INF,INF,INF};

        // get the edges for which w is the target
        for( auto eit=D.negative_cone_edges_begin(w);
             eit!=D.negative_cone_edges_end(w); ++eit )
        {
            auto e = *eit;
            auto v = D.source(e);
            size_t cone = getCone( w, v, P ),
                   flattened_cone = cone/2;
            auto e_bisector_length = bisectorLength(make_pair(w,v),P);
            if( e_bisector_length < local_minimum_bisectors[flattened_cone] )
            {
                local_anchors[flattened_cone] = v;
                local_minimum_bisectors[flattened_cone] = e_bisector_length;
            }
            cout<<"  -"<<D.source(e)<<" is in cone "<<cone<<endl;
        }

        // set the anchor for each negative cone
        for( size_t i=0; i<3; ++i )
        {
            if( local_minimum_bisectors[i] < INF )
            {
                // i is equal to floor( cone of edge (v,w) / 2 ), so blue cone (1) is at i==0
                Color c = getColor( local_anchors[i], w, P );//same but less clear -> static_cast<Color>(int(i!=0));
                anchors[c].emplace_back( local_anchors[i], w );
                cout<<"  added "<< (c==BLUE?"blue":"white")<<" anchor "<<local_anchors[i]<<"-"<< w<<endl;
            }
        }

    }
}
template< class Set>
inline bool vertexIsEmpty(const Vertex_descriptor vertex, const bool isHigherVertex, const Set &AdjacentWhiteConeIsEmpty)
{
    bool vertexIsKnown = contains( AdjacentWhiteConeIsEmpty, vertex );
    return // vertex with higher y value needs its negative (second) cones checked
           (!vertexIsKnown)
        || ( isHigherVertex && AdjacentWhiteConeIsEmpty.at(vertex).negative)
        || (!isHigherVertex && AdjacentWhiteConeIsEmpty.at(vertex).positive);
}
template< class Set>
inline void fillVertex(const Vertex_descriptor vertex, const bool isHigherVertex, Set &AdjacentWhiteConeIsEmpty)
{
    typedef typename Set::mapped_type PositiveAndNegativeWhiteConesForVertex;
    AdjacentWhiteConeIsEmpty.try_emplace(vertex, PositiveAndNegativeWhiteConesForVertex() );
    AdjacentWhiteConeIsEmpty[vertex].negative = !isHigherVertex
                                                && AdjacentWhiteConeIsEmpty[vertex].negative;
    AdjacentWhiteConeIsEmpty[vertex].positive =  isHigherVertex
                                                && AdjacentWhiteConeIsEmpty[vertex].positive;
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

    typedef unordered_map<Vertex_descriptor,
                          PositiveAndNegativeWhiteConesForVertex> PolarVertexStatusMap;
    PolarVertexStatusMap AdjacentWhiteConeIsEmpty;

    for( auto whiteAnchor : whiteAnchors )
    {
        cout<<"WhiteAnchor:"<<whiteAnchor.first<<"-"<<whiteAnchor.second<<"\n";

        Vertex_descriptor source = whiteAnchor.first,
                          target = whiteAnchor.second;
        bool sourceIsHigher = P[target].y() < P[source].y();

        if( vertexIsEmpty(source, sourceIsHigher,AdjacentWhiteConeIsEmpty)
         && vertexIsEmpty(target,!sourceIsHigher,AdjacentWhiteConeIsEmpty) )
        {
            cout<<"Both agree\n";
            A.insert(whiteAnchor);

            fillVertex(source,  sourceIsHigher, AdjacentWhiteConeIsEmpty);
            fillVertex(target, !sourceIsHigher, AdjacentWhiteConeIsEmpty);
        }
    }
}

template< class Triangulation, class EdgeList, class PointContainer, class AnchorList, class AdjacencyList >
void addBlueCanonicalEdges( const Triangulation &D, const EdgeList &A, const PointContainer &P, const AnchorList &blueAnchors, AdjacencyList &S )
{
    for( auto blueAnchor : blueAnchors )
    {
        typedef typename Triangulation::Edge_descriptor Edge_descriptor;
        vector<Edge_descriptor> fan;
        //auto source = blueAnchor.first;
        auto target = blueAnchor.second;
        D.fan_of_cone(target,1,fan);

        auto eit = fan.begin();
        auto previousVertex = D.source(*eit);
        // get fan of blueAnchor's target vertex cone 1 (negative/blue)
        for( ++eit; eit!=fan.end(); ++eit )
        {
            auto thisVertex = D.source(*eit);

            pair<Vertex_descriptor,Vertex_descriptor> canonicalEdge;
            bool edgeExists;
            boost::tie( canonicalEdge, edgeExists )
                = D.eitherEdge( thisVertex, previousVertex );

            if(edgeExists && !contains(A, canonicalEdge) )
            {
                cout<<"Canonical edge "<<canonicalEdge.first<<"-"<<canonicalEdge.second<<"\n";
                S.try_emplace( canonicalEdge.second, typename AdjacencyList::mapped_type() );
                S[canonicalEdge.second].insert(canonicalEdge.first);
            }
            previousVertex = thisVertex;
        }
    }
}

template< class Triangulation, class EdgeList, class PointContainer, class AnchorList, class AdjacencyList >
void addWhiteCanonicalEdges( const Triangulation &D, const EdgeList &A, const PointContainer &P, const AnchorList &whiteAnchors, AdjacencyList &S )
{
    for( auto whiteAnchor : whiteAnchors )
    {
        typedef typename Triangulation::Edge_descriptor Edge_descriptor;

        auto source = whiteAnchor.first;
        auto target = whiteAnchor.second;
        size_t cone = D.getCone(target,source);
        vector<Edge_descriptor> fan;
        D.fan_of_cone(target,cone,fan);

        auto boundaryVertex = D.source(fan.back()),
             previousVertex = D.source(fan.front());
        int direction = 1;

        if( (cone+1) % 3 == 0 )
        { //white side of anchor is CW
            direction = -1;
            swap(previousVertex,boundaryVertex);
        }

        auto eit = fan.begin();
        while( D.source(*eit) != source )
            ++eit;

        previousVertex = D.source(*eit);

        for( eit+=direction; previousVertex!=boundaryVertex; eit+=direction )
        {
            auto thisVertex = D.source(*eit);
            cout<<"Considering white canonical edge "<<previousVertex<<"-"<<thisVertex<<endl;

            pair<Vertex_descriptor,Vertex_descriptor> canonicalEdge;
            bool edgeExists;
            boost::tie( canonicalEdge, edgeExists )
                = D.eitherEdge( thisVertex, previousVertex );

            // For white cones, add can. edge if its anchor is not in A
            if(edgeExists && !contains(A, canonicalEdge) )
            {
                cout<<"Canonical edge "<<canonicalEdge.first<<"-"<<canonicalEdge.second<<"\n";
                S.try_emplace( canonicalEdge.second, typename AdjacencyList::mapped_type() );
                S[canonicalEdge.second].insert(canonicalEdge.first);
            }
            previousVertex = thisVertex;
        }
    }
}

template<class AdjacencyList >
void createShortcut( const Vertex_descriptor &p, const Vertex_descriptor &q, const Vertex_descriptor &r, AdjacencyList &adj) {

    adj[q].erase(p);
    adj[q].erase(r);

    adj.try_emplace(r,typename AdjacencyList::mapped_type());
    adj[r].insert(p);
}

template<class Triangulation, class AdjacencyList>
void addBlueShortcuts( const Triangulation &D, AdjacencyList &S_not_A )
{
    for( auto v : S_not_A )
    {
        // check if there are two edges that share a target
        cout<<"  |"<<v.first<<"|="<<v.second.size()<<"\n";
        if( v.second.size() > 1 )
        {
            auto p = *v.second.begin(),
                 q = v.first,
                 r = *(--v.second.end());

            cout<< "    Blue shortcut edge "<<*v.second.begin()
                <<"-"<<v.first<<"-"<<*(--v.second.end())<< " --> "
                <<*v.second.begin()<<"-"<<*(--v.second.end())<<"\n";

            assert( D.edgeExists( make_pair(p,q)) && D.edgeExists( make_pair(r,q)) );
            createShortcut( p, q, r, S_not_A );
        }
    }
}



template< class Triangulation, class EdgeList, class PointContainer, class AnchorList, class AdjacencyList >
void addWhiteShortcuts( const Triangulation &D, const EdgeList &A, const PointContainer &P, const AnchorList &whiteAnchors, AdjacencyList &S )
{
    for( auto whiteAnchor : whiteAnchors )
    {
        typedef typename Triangulation::Edge_descriptor Edge_descriptor;

        auto source = whiteAnchor.first;
        auto target = whiteAnchor.second;
        cout<<"Looking for shortcuts of white anchor "<<source<<"-"<<target<<endl;
        size_t cone = D.getCone(target,source);
        vector<Edge_descriptor> fan;
        D.fan_of_cone(target,cone,fan);

        auto boundaryVertex = D.source(fan.back()),
             previousVertex = D.source(fan.front());
        int direction = 1;

        if( (cone+1) % 3 == 0 ) { //white side of anchor is CW
            direction = -1;
            swap(previousVertex,boundaryVertex);
        }
        //cout<<"direction="<<direction<<endl;

        if( source != boundaryVertex )
        {
            auto i = fan.begin(); // use i here to be consistent with paper
            while( D.source(*i) != source )
                ++i;

            previousVertex = D.source(*i);

            //cout<<"index in fan="<<(i-fan.begin())<<endl;

            i+=direction;

            for( auto i_vertex = D.source(*i);
                 i_vertex != boundaryVertex;
                 i_vertex = D.source(*i) )
            {
                //cout<<"i_vertex="<<i_vertex<<endl;

                if( getColor(previousVertex,i_vertex,P) == WHITE )
                {
                    //cout<<"color is white"<<endl;
                    i+=direction;
                    previousVertex = i_vertex;
                }
                else if( i_vertex != boundaryVertex )
                {
                    //cout<<"color is blue"<<endl;
                    auto j = i;
                    j+=direction;

                    auto j_vertex = D.source(*j),
                         j_prev = i_vertex,
                         max_j = j;

                    auto j_angle = get_angle<K>(P[target], P[i_vertex], P[j_vertex]),
                         maxAngle = j_angle;

                    for( j=j; j_prev!=boundaryVertex; j+=direction )
                    {    // j is in a white cone AND all vertices from i+1 to boundaryVertex are on the same side of i,j
                        j_vertex = D.source(*j);
                        j_angle = get_angle<K>(P[target], P[i_vertex], P[j_vertex]);

                        if( getColor(j_prev,j_vertex,P) == WHITE
                         && maxAngle < j_angle )
                        {
                            maxAngle = j_angle;
                            max_j = j;
                        }
                        j_prev = j_vertex;
                    }
                    i = max_j;
                    j_vertex = D.source(*max_j);
                    j_prev = D.source(*(max_j-=direction));

                    cout<<"White shortcut edge "<<j_vertex<<"-"<<i_vertex<<"\n";
                    S.try_emplace( i_vertex, typename AdjacencyList::mapped_type() );
                    S[i_vertex].insert(D.source(*max_j));
                    if(S.find(j_prev) != S.end())
                    {
                        S[j_prev].erase(j_vertex);
                    }
                    previousVertex = j_prev;
                }
            }
        }
    }
}

} // namespace KPT2017


// Main algorithm.
template<typename RandomAccessIterator, typename OutputIterator>
void KPT2017(RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false)
{
    using namespace kpt2017;

    vector<Point_2> P( pointsBegin, pointsEnd );

    TD_Delaunay_2 D( P.begin(), P.end() );
    size_t n = D.number_of_vertices();

    map<Color,vector<pair<size_t,size_t>>> Anchors;
    findAnchors( Anchors, D, P );

    auto AnchorComp = [&P]( const pair<size_t,size_t> &lhs, const pair<size_t,size_t> &rhs )
    {
        return   lhs.second  < rhs.second
            || ( lhs.second == rhs.second
                && getCone(lhs.second, lhs.first,P) < getCone(rhs.second,rhs.first,P) );
    };
    // Step 1. add all blue anchors to A
    set<pair<size_t,size_t>,decltype(AnchorComp)> A(
        Anchors[BLUE].begin(),
        Anchors[BLUE].end(),
        AnchorComp );

    // Step 2.
    addWhiteAnchors( Anchors[WHITE], P, A );

    // A should be degree 4 or less
    vector<size_t> degree(n,0);

    for( auto a : A ) {
        ++degree[a.first];
        ++degree[a.second];
    }
    assert( *max_element(degree.begin(), degree.end()) <= 4 );

    // Step 3.
    typedef map<size_t, set<size_t>> AdjacencyList;
    AdjacencyList S_not_A;
    // Add to S every canonical edge in negative blue cones (cone 1) if the edge isn't in A
    addBlueCanonicalEdges(D,A,P,Anchors[BLUE],S_not_A);

    //Step 4.
    addBlueShortcuts(D,S_not_A);

    // should be degree 4 or less
    for( auto s : S_not_A ) {
        for( auto a : s.second ) {
            ++degree[a];
            ++degree[s.first];
        }
    }
    assert( *max_element(degree.begin(), degree.end()) <= 4 );

    //Step 5.
    addWhiteCanonicalEdges(D,A,P,Anchors[WHITE],S_not_A);

    //Step 6.
    addWhiteShortcuts(D,A,P,Anchors[WHITE],S_not_A);

    // Since we didn't set S to A in step 3, we need to combine them now
    set<pair<size_t,size_t>> S( A.begin(), A.end() );

    cout<<"Adding S_not_A to S...\n";
    for( auto v : S_not_A )
    {
        if( v.second.size() > 0)
        {
            cout<<*v.second.begin()<<"-"<< v.first<<endl;
            S.emplace( *v.second.begin(), v.first );
        }
    }

    // should be degree 4 or less
    fill( degree.begin(), degree.end(), 0 );
    for( auto s : S ) {
        ++degree[s.first];
        ++degree[s.second];
    }
    assert( *max_element(degree.begin(), degree.end()) <= 4 );

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
        GraphPrinter printer(9);
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
