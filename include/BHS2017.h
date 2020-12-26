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

    unordered_map<pair<size_t,size_t>,double,pointPairHash> B(L.size());

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

        size_t p = e.first;
        size_t r = e.second;

        size_t cone  = getCone(handles, p, r, alpha);

        /*Vertex circulator oriented to r to find fist and last end vertex. Once found all neighbors are added
          starting from first end to last end vertex into a vector.
        */
        auto N_p = DT.incident_vertices(handles[p]);

        while(++N_p != handles[r]);

        while(!DT.is_infinite(++N_p) && getCone(handles, p, N_p->info(), alpha) == cone && !( B.find({p, N_p->info()}) == B.end() && B.find({N_p->info(), p}) == B.end()) &&
            ( B.at(make_pair(p,N_p->info())) > B.at(e) || abs( B.at(make_pair(p,N_p->info())) - B.at(e)) < EPSILON ) );

        vector<size_t> canNeighbors;

        while(!DT.is_infinite(--N_p) && getCone(handles, p, N_p->info(), alpha) == cone && !( B.find({p, N_p->info()}) == B.end() && B.find({N_p->info(), p}) == B.end()) &&
            ( B.at(make_pair(p,N_p->info())) > B.at(e) || abs( B.at(make_pair(p,N_p->info())) - B.at(e)) < EPSILON ) ){

            canNeighbors.push_back(N_p->info());
        }

        //Add inner edges if total neigborhood edges is 3 or more. (4.2)
        for(int i=1; i<int(canNeighbors.size())-2; i++){
            E_CAN.emplace_back(canNeighbors.at(i),canNeighbors.at(i+1));
        }

        //If r is an end edge add the edge with endpoint r. (4.3)
        int canSize = canNeighbors.size()-1;

        if(canNeighbors.at(0) == r && canSize > 0){
            E_CAN.emplace_back(r,canNeighbors.at(1));
        }
        if(canNeighbors.at(canSize)== r && canSize > 0){
            E_CAN.emplace_back(canNeighbors.at(canSize), canNeighbors.at(canSize-1));
        }

    }

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
