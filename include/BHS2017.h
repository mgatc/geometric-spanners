#ifndef GSNUNF_BHS2017_H
#define GSNUNF_BHS2017_H

#include <cmath>         // ceil, floor, isinf
#include <limits>
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

#include <boost/functional/hash.hpp> // size_t pair hash

#include <CGAL/algorithm.h> //
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Line_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

//#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "metrics.h"
#include "utilities.h"


namespace gsnunf {

using namespace std;

namespace bhs2017 {

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K>      Vb;
typedef CGAL::Triangulation_face_base_2<K>                          Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef CGAL::Aff_transformation_2<K>                               Transformation;
typedef Delaunay::Vertex_handle                                     Vertex_handle;
typedef Delaunay::Vertex_circulator                                 Vertex_circulator;
typedef CGAL::Vector_2<K>                                           Vector_2;
typedef CGAL::Line_2<K>                                             Line;
typedef Delaunay::Point                                             Point;
typedef Delaunay::Finite_vertices_iterator                          Finite_vertices_iterator;
typedef Delaunay::Finite_edges_iterator                             Finite_edges_iterator;

typedef pair<size_t,size_t>                                         size_tPair;
typedef boost::hash<size_tPair>                                     size_tPairHash;
typedef unordered_map<size_tPair,bool,size_tPairHash>               size_tPairMap;
typedef unordered_map<pair<size_t,size_t>, double, pointPairHash, edgeEquality> edgeBisectorMap;


//Function to find the cone of p containg vertex q, for this algorithm all vertices have 6 cones (0-5) with an angle of (PI/3).
inline size_t getCone(const vector<Vertex_handle>& H, const size_t p, const size_t q, const double alpha){

    double tan30 = tan(PI/6);

    Point refPoint(H[p]->point().x() - tan30, H[p]->point().y() + 1);

    double theta = get_angle<bhs2017::K>(refPoint, H[p]->point(), H[q]->point());

    size_t cone = ( theta / alpha );

    return cone;
}


//Step 2. Fins the bisector length of all edges in the delaunay triangulation.
inline K::FT bisectorLength( const vector<Vertex_handle>& H, const pair<size_t,size_t>& e, double alpha) {


    double tan30 = tan(PI/6);
    double cot30 = 1/tan30;

    vector<double> bisectorSlopes = {INF, tan30, -1*tan30, INF, tan30, -1*tan30};
    vector<double> orthBisectorSlopes{0, -1*cot30, cot30, 0, -1*cot30, cot30};

    size_t cone = getCone(H, e.first, e.second, alpha);

    double xCord = H[e.first]->point().x();
    double yCord = H[e.first]->point().y() + 1;

    if(cone % 3 != 0){
        xCord = (bisectorSlopes[cone] * H[e.first]->point().x() + 1) / bisectorSlopes[cone];
    }

    Point bisectorPoint(xCord, yCord);

    Line bisectorLine(H[e.first]->point(), bisectorPoint);

    Point intersectionPoint = bisectorLine.projection(H[e.second]->point());

    double bisectorLen = distance( H[e.first]->point(), intersectionPoint );

    return bisectorLen;
}

//Step 4 add cannonical edges. Consists of 4 sub steps.
inline void addCannonical(size_t p, size_t r, vector<pair<size_t, size_t>> &E_CAN, const vector<pair<size_t, size_t>> &E_A, const edgeBisectorMap &B,
 const vector<bitset<6>> &coneStatus, const Delaunay &DT, const vector<Vertex_handle> &handles, const size_t n, const double alpha){

    pair<size_t, size_t> e = make_pair(p,r);

    size_t cone  = getCone(handles, p, r, alpha);

    /*Vertex circulator oriented to r to find fist and last end vertex. Once found all neighbors are added
      starting from first end to last end vertex into a vector.
    */
    auto N_p = DT.incident_vertices(handles[p]);

    while(++N_p != handles[r]);

    while(!DT.is_infinite(++N_p) && getCone(handles, p, N_p->info(), alpha) == cone &&
        ( B.at(make_pair(p,N_p->info())) > B.at(e) || abs( B.at(make_pair(p,N_p->info())) - B.at(e)) < EPSILON ) );

    vector<size_t> canNeighbors;

    while(!DT.is_infinite(--N_p) && getCone(handles, p, N_p->info(), alpha) == cone &&
        ( B.at(make_pair(p,N_p->info())) > B.at(e) || abs( B.at(make_pair(p,N_p->info())) - B.at(e)) < EPSILON ) ){

        canNeighbors.push_back(N_p->info());
    }

//        cout << "CAN NEIGHBORS\n";
//        for( auto e: canNeighbors){
//            cout << e << "\n";
//        }

    //Add inner edges if total neigborhood edges is 3 or more. (4.2)
    for(int i=1; i<int(canNeighbors.size())-2; i++){
        E_CAN.emplace_back(canNeighbors.at(i),canNeighbors.at(i+1));
    }

    //If r is an end edge add the edge with endpoint r. (4.3)
    int canEdges = canNeighbors.size()-1;

    if(canNeighbors.front() == r && canEdges > 1){
        E_CAN.emplace_back(r,canNeighbors.at(1));
    }
    if(canNeighbors.back() == r && canEdges > 1){
        E_CAN.emplace_back(canNeighbors.back(), canNeighbors.at(canEdges-1));
    }

    //First and last edges in the cannonical neighborhood are condidered and added by 3 criteria. (4.4)
    if(canEdges > 0){

        pair<size_t,size_t> canFirst = make_pair(canNeighbors.front(), canNeighbors[1]);
        pair<size_t, size_t> canLast = make_pair(canNeighbors[canEdges-1], canNeighbors.back());

        //If the edges are in cone 1 or 5 with respect to a and z add. (4.4 a)
        if(getCone(handles, canFirst.first, canFirst.second, alpha) == (cone + 1) %6){
            E_CAN.push_back(canFirst);
        }
        if(getCone(handles, canLast.second, canLast.first, alpha) == (cone + 5) %6){
            E_CAN.push_back(canLast);
        }

        //If the edges are in cone 2 or 4 with respect to a and z and cone for has no edge with an end edge point in E_A add. (4.4 b)
        if(getCone(handles, canFirst.first, canFirst.second, alpha) == (cone + 2) % 6 && !coneStatus.at(canFirst.first)[(cone + 2) %6]){
                E_CAN.push_back(canFirst);
            }
        if(getCone(handles, canLast.second, canLast.first, alpha) == (cone + 4) % 6 && !coneStatus.at(canLast.second)[(cone + 4) %6]){
                E_CAN.push_back(canLast);
            }

        //Checks if end edges have a end point a or z in E_A and an edge different from the end edge in cone 2 or 4 woth respect to a and z. (4.4 c)
        if(getCone(handles, canFirst.first, canFirst.second, alpha) == (cone + 2) % 6){

        for(auto e : E_A){
                if(e.first == canFirst.first && getCone(handles, canFirst.first, e.second, alpha) == (cone + 2) % 6){
                    if(e.second != canFirst.second){
                        E_CAN.emplace_back(make_pair(canFirst.second, e.second));
                    }
                }

                if(e.second == canFirst.first && getCone(handles, canFirst.first, e.first, alpha) == (cone + 2) % 6){
                    if(e.first != canFirst.second){
                        E_CAN.emplace_back(make_pair(canFirst.second, e.first));
                    }
                }
            }

//            auto N_a = DT.incident_vertices(handles[canFirst.first]);
//
//            while(--N_a != handles[canFirst.second]);
//
//            if(!DT.is_infinite(--N_a) && getCone(handles, canFirst.first, N_a->info(), alpha) == (cone + 2) % 6){
//                E_CAN.emplace_back(make_pair(canFirst.second, N_a->info()));
//                cout << "Added (B,C)";
//            }
        }

        if(getCone(handles, canLast.second, canLast.first, alpha) == (cone + 4) % 6){

            for(auto e : E_A){
                if(e.first == canLast.second && getCone(handles, canLast.second, e.second, alpha) == (cone + 4) % 6){
                    if(e.second != canLast.first){
                        E_CAN.emplace_back(make_pair(e.second, canLast.first));
                    }
                }

                if(e.second == canLast.second && getCone(handles, canLast.second, e.first, alpha) == (cone + 4) % 6){
                    if(e.first != canLast.first){
                        E_CAN.emplace_back(make_pair(e.first, canFirst.first));
                    }
                }
            }

//            auto N_z = DT.incident_vertices(handles[canLast.second]);
//
//            while(--N_z != handles[canLast.first]);
//
//            if(!DT.is_infinite(++N_z) && getCone(handles, canLast.second, N_z->info(), alpha) == (cone + 4) % 6){
//                E_CAN.emplace_back(make_pair(N_z->info(), canLast.first));
//                cout << "Added (W,Y)";
//            }
        }
    }

}


} // namespace BHS2017

template< typename RandomAccessIterator, typename OutputIterator >
void BHS2017( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false ) {
    using namespace bhs2017;

    const double alpha = PI / 3;

    //if(printLog) cout<<"alpha:"<<alpha<<",";

    // Construct Delaunay triangulation
    bhs2017::Delaunay DT( pointsBegin, pointsEnd );
    size_t n = DT.number_of_vertices();
    if( n > SIZE_T_MAX - 1 ) return;

    vector<bhs2017::Vertex_handle> handles(n);

    // Add IDs
    size_t i=0;
    for( auto v=DT.finite_vertices_begin(); v!=DT.finite_vertices_end(); ++v ) {
        v->info() = i;
        handles[i] = v;
        ++i;
    }

    // Put edges in a vector, then sort on weight
    vector<pair<size_t,size_t>> L;

    for( auto e=DT.finite_edges_begin(); e!=DT.finite_edges_end(); ++e ) {
        L.emplace_back( e->first->vertex( (e->second+1)%3 )->info(),
                        e->first->vertex( (e->second+2)%3 )->info() );
    }

    edgeBisectorMap B(L.size());

    for(auto e : L){
        B.emplace(e, bisectorLength(handles, e, alpha));
    }

    sort( L.begin(), L.end(), [&]( const auto& lhs, const auto& rhs ) {
        return B.at(lhs) < B.at(rhs);
    });

    //Add incident set E_A.
    vector<pair<size_t,size_t> > E_A;

    vector<bitset<6>> coneStatus(n);

    //Loop through all edges ordered in L.
    for(auto e : L){

        size_t p = e.first;
        size_t q = e.second;

        //Check if cone i of p containing q has no edges in E_A with end point q in the same cone neighborhood.
        size_t p_cone = getCone(handles, p, q, alpha);

        bool p_cone_status = coneStatus.at(p)[p_cone];
        bool q_cone_status = coneStatus.at(q)[(p_cone + 3) % 6];

        if(!p_cone_status && !q_cone_status){
            E_A.push_back(e);
            coneStatus.at(p)[p_cone] = true;
            coneStatus.at(q)[(p_cone + 3) % 6] = true;
        }
    }

    //Add cannonical E_CAN
    vector<pair<size_t, size_t>> E_CAN;

    for(auto  e : E_A){

    addCannonical(e.first, e.second, E_CAN, E_A, B, coneStatus, DT, handles, n, alpha);
    addCannonical(e.second, e.first, E_CAN, E_A, B, coneStatus, DT, handles, n, alpha);

    }

//    cout << "E_A\n";
//    for(auto e: E_A){
//        cout << e.first << " " << e.second << "\n";
//    }
//
//    cout << "\nE_CAN\n";
//    for(auto e: E_CAN){
//        cout << e.first << " " << e.second << "\n";
//    }

//    cout << "\nUnique\n";

    //Union of sets E_A and E_CAN for final edge set.
    // Combine E and E_star, remove duplicates
    E_A.insert( E_A.end(), E_CAN.begin(), E_CAN.end() );
    sort( E_A.begin(), E_A.end(), []( const auto& l, const auto& r) {
        return (min(l.first, l.second) < min(r.first, r.second)
            || (min(l.first, l.second) == min(r.first, r.second) && max(l.first, l.second) < max(r.first, r.second)));
    } );
    E_A.erase( unique(E_A.begin(), E_A.end(), []( const auto& l, const auto& r) {
        return ( l.first == r.first && l.second == r.second )
            || ( l.first == r.second && l.second == r.first );
    }), E_A.end() );

//    cout << "\nE\n";
//    for(auto e: E_A){
//        cout << e.first << " " << e.second << "\n";
//    }

    // Edge list is only needed for printing. Remove for production.
    vector< pair<Point,Point> > edgeList;

    // Send resultant graph to output iterator
        for( auto e : E_A ) {
            // Edge list is only needed for printing. Remove for production.
            edgeList.emplace_back( handles.at(e.first)->point(), handles.at(e.second)->point() );

            *result = make_pair( handles.at(e.first)->point(), handles.at(e.second)->point() );
            ++result;
            *result = make_pair( handles.at(e.second)->point(), handles.at(e.first)->point() );
            ++result;
        }

        // START PRINTER NONSENSE
        if( printLog ) {
            GraphPrinter printer(1);
            GraphPrinter::OptionsList options;

            options = {
                { "color", printer.inactiveEdgeColor },
                { "line width", to_string(printer.inactiveEdgeWidth) }
            };
            printer.drawEdges( DT, options );

            options = { // active edge options
                { "color", printer.activeEdgeColor },
                { "line width", to_string(printer.activeEdgeWidth) }
            };
            printer.drawEdges( edgeList.begin(), edgeList.end(), options );


            options = {
                { "vertex", make_optional( to_string(printer.vertexRadius) ) }, // vertex width
                { "color", make_optional( printer.backgroundColor ) }, // text color
                { "fill", make_optional( printer.activeVertexColor ) }, // vertex color
                { "line width", make_optional( to_string(0) ) } // vertex border (same color as text)
            };
            GraphPrinter::OptionsList borderOptions = {
                { "border", make_optional( to_string(printer.vertexRadius) ) }, // choose shape of vertex
                { "color", printer.activeEdgeColor }, // additional border color
                { "line width", to_string(printer.inactiveEdgeWidth) }, // additional border width
            };
            printer.drawVerticesWithInfo( DT, options, borderOptions );

            printer.print( "BHS2017" );
            cout<<"\n";
        }
        // END PRINTER NONSENSE

} // function BHS2017

} // namespace gsnunf

#endif // GSNUNF_BHS2017_H
