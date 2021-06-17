//Needs optimizing currently testing.
#ifndef GSNUNF_BGHP2010_H
#define GSNUNF_BGHP2010_H

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

namespace bghp2010 {

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
typedef Epick                     K;
typedef K::Point_2                Point_2;
typedef K::FT                     FT;

typedef HalfThetaTriangulation<K> TD_Delaunay_2;
typedef TD_Delaunay_2::Vertex_descriptor Vertex_descriptor;

enum EdgeLabel { // per cone
    FIRST = 1,
    CLOSEST = 0,
    LAST = -1
};

enum ConePolarity {
    POSITIVE = 0, // even cones
    NEGATIVE = 1  // odd cones
};

inline ConePolarity getConePolarity( const size_t p, const size_t q, const vector<Point_2> &h )
{
    return ConePolarity( getCone(p,q,h) % 2 );
}

template< class Container >
inline double get_canonical_angle( const size_t p, const size_t q, const size_t r, const Container &P )
{
    return CGAL::min( get_angle(p,q,r,P), get_angle(r,q,p,P) );
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

inline size_t iPlusOne(const size_t i)
{
    return ((i+1) % 6);
}

inline size_t iLessOne(const size_t i)
{
    return ((i-1+6) % 6);
}

template< class KeyEdgesMap, class Triangulation, class PointContainer >
void findKeyEdges( KeyEdgesMap &KeyEdges, Triangulation &D, const PointContainer &P )
{
    typedef typename Triangulation::Edge_descriptor Edge_descriptor;

    // find closest in each cone
    for( auto vit=D.finite_vertices_begin();
         vit!=D.finite_vertices_end(); ++vit )
    {
        auto w = *vit;

        for( size_t cone=1; cone<6; cone+=2 )
        {
            vector<Edge_descriptor> fan;
            D.fan_of_cone(w,cone,fan);

            // Find closest
            Vertex_descriptor closestKnown;
            K::FT closestKnownBisector = INF;

            for( auto e : fan )
            {
                auto v = D.source(e);
                auto length = bisectorLength( make_pair(w,v), P );
                if( length < closestKnownBisector )
                {
                    closestKnown = v;
                    closestKnownBisector = length;
                }
            }
            size_t flattenedCone = cone/2;
            if( closestKnownBisector < INF )
            {
                KeyEdges[CLOSEST][w][flattenedCone] = closestKnown;
            }

            // Find first and last
            if( fan.size() > 0 )
            {
                KeyEdges[FIRST][w][flattenedCone] = D.source( fan.back() );
                KeyEdges [LAST][w][flattenedCone] = D.source( fan.front() );
            }
        }
    }
}
template< class Triangulation, class PointContainer, class KeyEdgesContainer >
bool is_i_relevant( const Vertex_descriptor child,
                    const Vertex_descriptor parent,
                    const Vertex_descriptor grandparent,
                    size_t i,
                    Triangulation &D,
                    const PointContainer &P,
                    KeyEdgesContainer &KeyEdges )
{
    if( child == SIZE_T_MAX || grandparent == SIZE_T_MAX )
        return false;

    size_t child_cone = getCone(child,grandparent,P);

    if( child_cone != i )
        return false;

    return ( child == KeyEdges[FIRST][parent][iPlusOne(i)/2] && child != KeyEdges[CLOSEST][parent][iPlusOne(i)/2] )
        || ( child == KeyEdges [LAST][parent][iLessOne(i)/2] && child != KeyEdges[CLOSEST][parent][iLessOne(i)/2] );
}

template< class Triangulation, class PointContainer, class EdgeList, class KeyEdgesContainer >
bool is_i_distant( const Vertex_descriptor w,
                   const size_t i,
                   const size_t charge,
                   Triangulation &D,
                   const PointContainer &P,
                   const EdgeList &E,
                   KeyEdgesContainer &KeyEdges )
{
    if( charge < 2 )
        return false;

    Vertex_descriptor u = D.parent(w,i),
        first = KeyEdges[FIRST][w][iPlusOne(i)/2],
        last = KeyEdges [LAST][w][iLessOne(i)/2];

    //cout<< "w:"<<w<<" u:"<<u<<" i:"<<i<< " first:"<<first<<" last:"<<last<<endl;

    //return false;
    return !contains( E, make_pair(w,u) )
        && is_i_relevant( first, w, u, iPlusOne(i), D, P, KeyEdges )
        && is_i_relevant( last,  w, u, iLessOne(i), D, P, KeyEdges );
}

template< class AdjacencyList, class PointSet, class ChargeList, class EdgeList>
void addClosestInNegativeCones(const AdjacencyList &ClosestEdges,
                               const PointSet &P,
                               ChargeList &Charges,
                               EdgeList &E,
                               bool printLog)
{
    const size_t n = P.size();
    // Add closest in each negative cone for each vertex
    if(printLog) cout<<"\nClosest\n";
    for( size_t u=0; u<n; ++u ) {
        for( auto v : ClosestEdges[u] ) {
            if( v != SIZE_T_MAX )
            {
                size_t cone = getCone(u,v,P);
                auto pair1 = make_pair(u,cone),
                     pair2 = make_pair(u,(cone+3)%6);
                Charges.try_emplace( pair1, 0 );
                Charges[pair1]++;
                Charges.try_emplace( pair2, 0 );
                Charges[pair2]++;

                E.emplace( u,v );
                if( printLog ) cout<<v<<" "<<u<<"\n";
            }
        }
    }
}

template< class KeyEdgeList, class PointSet, class Triangulation, class ChargeList, class EdgeList>
void add_i_relevantNeighbors(  KeyEdgeList &KeyEdges,
                               const PointSet &P,
                               const Triangulation &D,
                               ChargeList &Charges,
                               EdgeList &E,
                               bool printLog)
{
    const size_t n = P.size();

    // Add first and last in each negative cone if it is (i+1)-relevant
    if(printLog) cout<<"\nFirst and Last\n";
    for( size_t u=0; u<n; ++u ) {
        // get edges from positive cones
        for( auto it= D.positive_cone_edges_begin(u);
                  it!=D.positive_cone_edges_end(u);
                  ++it )
        {
            auto e = *it;
            auto v = D.target(e);
            size_t posCone = getCone(u,v,P);
            for( int j=-1; j<=1; ++j )
            {
                size_t i = (posCone + 6 + j) % 6;
                //size_t iLessOne = (i - 1 + 6) % 6;
                EdgeLabel w_label = static_cast<EdgeLabel>(j); // FIRST or LAST
                auto w = KeyEdges[w_label][u][i/2];

                if( is_i_relevant( w, u, v, posCone, D, P, KeyEdges ) ) {

                    size_t cone = getCone(u,v,P);
                    auto pair2 = make_pair(v,cone);
                    Charges.try_emplace( pair2, 0 );
                    Charges[pair2]++;

                    E.emplace( u,w );
                    if( printLog ) cout<<w<<" "<<u<<"\n";
                }
            }
        }
    }
}

template< class KeyEdgeList, class PointSet, class Triangulation, class ChargeList, class EdgeList>
void handle_i_distantCharge2s(  KeyEdgeList &KeyEdges,
                               const PointSet &P,
                               const Triangulation &D,
                               ChargeList &Charges,
                               EdgeList &E,
                               bool printLog)
{
    const size_t n = P.size();

    for( auto it=Charges.begin();
         it!=Charges.end(); ++it )
    {
        auto cone = *it;
        Vertex_descriptor w = cone.first.first;
        size_t i = cone.first.second;
        size_t& charge = cone.second;
        if( is_i_distant( w, i, charge, D, P, E, KeyEdges ) )
        {
            auto next = KeyEdges[FIRST][w][iPlusOne(i)/2],
                 prev = KeyEdges[LAST] [w][iLessOne(i)/2];
            E.emplace( next, prev );

            // remove the edge to w that is farthest from w's grandparent in cone i
            Vertex_descriptor u = D.parent(w,i);
            auto thetaNext = get_canonical_angle(w,u,next,P),
                 thetaPrev = get_canonical_angle(w,u,prev,P);
            auto remove = make_pair( w, thetaNext > thetaPrev ? next : prev );
            assert( contains(E,remove) );
            E.erase(remove);

            // update charges
            auto nextCharge = make_pair(next,iPlusOne(i)),
                 prevCharge = make_pair(prev,iLessOne(i));
            ++(Charges[nextCharge]);
            ++(Charges[prevCharge]);
            --charge;
        }
    }
}

template< class KeyEdgeList, class PointSet, class Triangulation, class ChargeList, class EdgeList>
void handleOtherCharge2s(  KeyEdgeList &KeyEdges,
                               const PointSet &P,
                               const Triangulation &D,
                               ChargeList &Charges,
                               EdgeList &E,
                               bool printLog)
{
    const size_t n = P.size();

    for( auto it=Charges.begin();
         it!=Charges.end(); ++it )
    {
        auto cone = *it;
        Vertex_descriptor w = cone.first.first;
        size_t i = cone.first.second;
        size_t& charge = cone.second;
        if( charge == 2
         && Charges[make_pair(w,iLessOne(i))] == 1
         && Charges[make_pair(w,iPlusOne(i))] == 1 )
        {
            auto w_parent = D.parent(w,i);
            if( w_parent == SIZE_T_MAX )
                continue;
            auto next = KeyEdges[FIRST][w][iPlusOne(i)/2],
                 prev = KeyEdges[LAST] [w][iLessOne(i)/2];

            auto remove = make_pair( w, w == KeyEdges[LAST][w_parent][i] ? prev : next );

            //assert( contains(E,remove) );

            if (printLog) cout<<"Edge "<<remove.first<<"-"<<remove.second<<" ";
            if( printLog && !contains(E,remove) ) cout<<"not ";
            if( printLog) cout<<"found\n";

            E.erase( remove );

            // update charge
            --charge;
        }
    }
}

} // namespace bghp2010


// Main algorithm.
template<typename RandomAccessIterator, typename OutputIterator>
void BGHP2010(RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false)
{
    using namespace bghp2010;

    // Step 1
    vector<Point_2> P( pointsBegin, pointsEnd );
    TD_Delaunay_2 D( P.begin(), P.end() );

    {
        Timer tim;
        size_t n = D.number_of_vertices();

        map<EdgeLabel, vector<vector<size_t>>> KeyEdges = {
            { CLOSEST, vector<vector<size_t>>(n, vector<size_t>(3, SIZE_T_MAX)) },
            { FIRST, vector<vector<size_t>>(n, vector<size_t>(3, SIZE_T_MAX)) },
            { LAST, vector<vector<size_t>>(n, vector<size_t>(3, SIZE_T_MAX)) }
        };

        findKeyEdges( KeyEdges, D, P );

        set<pair<size_t,size_t>> E;
        map<pair<Vertex_descriptor,size_t>, size_t> Charges;

        // Step 2
        addClosestInNegativeCones( KeyEdges[CLOSEST], P, Charges, E, printLog );
        add_i_relevantNeighbors( KeyEdges, P, D, Charges, E, printLog );

        // Sanity check after step 2...
        //      each negative (odd) cone should have at most 1 edge,
        //      each positive (even) cone should have at most 2 edges

        if(printLog)cout<<"\nCharges after step 2\n";
        for(auto charge: Charges)
        {
            size_t cone = charge.first.second;
            if(printLog)cout<< charge.first.first<< "  "<<cone<<"   "<<charge.second<<"   "<<((cone%2==1&&charge.second <= 1)||(cone%2==0 && charge.second<=2)? "OK":"FAIL")<<"\n";
        }

        // Step 3
        handle_i_distantCharge2s( KeyEdges, P, D, Charges, E, printLog );

        if(printLog)cout<<"\nCharges after step 3\n";
        for(auto charge: Charges)
        {
            size_t cone = charge.first.second;
            if(printLog)cout<< charge.first.first<< "  "<<cone<<"   "<<charge.second<<"   "<<((cone%2==1&&charge.second <= 1)||(cone%2==0 && charge.second<=2)? "OK":"FAIL")<<"\n";
        }
        cout<<endl;

        // Step 4
        handleOtherCharge2s( KeyEdges, P, D, Charges, E, printLog );

        if(printLog)cout<<"\nCharges after step 4\n";
        for(auto charge: Charges)
        {
            size_t cone = charge.first.second;
            if(printLog)cout<< charge.first.first<< "  "<<cone<<"   "<<charge.second<<"   "<<((cone%2==1&&charge.second <= 1)||(cone%2==0 && charge.second<=2)? "OK":"FAIL")<<"\n";
        }

        // Send resultant graph to output iterator
        for(auto e : E)
        {
            *result = e;
            ++result;
//            *result = reverse_pair(e);
//            ++result;
        }





        // START PRINTER NONSENSE
        if(printLog)
        {
            vector<pair<Point_2,Point_2>> edgeList;

            for(auto e : E)
            {
                edgeList.emplace_back(P.at(e.first), P.at(e.second));
            }

            GraphPrinter printer(0.01);
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

            printer.print("BGHP2010");
            cout << "\n";
        }
        // END PRINTER NONSENSE
    }
} // function BGHP2010

} // namespace gsnunf

#endif // GSNUNF_BGHP2010_H
