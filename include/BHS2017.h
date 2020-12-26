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

bool selectEdge( const Delaunay& T, size_tPairMap &E, const Vertex_handle i, const Vertex_handle j, const size_t n, bool printLog = false ) {
    assert( T.is_edge( i, j ) );
    //if( printLog ) cout<<"add:("<<i->info()<<","<<j->info()<<") ";

    auto existing = E.begin();
    bool inserted = false;
    tie(existing,inserted) = E.try_emplace( makeNormalizedPair( i->info(), j->info() ), false );
    if(!inserted) existing->second = true;

    return inserted;
}

inline K::FT bisectorLength( const vector<Vertex_handle>& H, const pair<size_t,size_t>& e, const double alpha) {

    double tan30 = tan(PI/6);
    double cot30 = 1/tan30;

    vector<double> bisectorSlopes = {INF, tan30, -1*tan30, INF, tan30, -1*tan30};
    vector<double> orthBisectorSlopes{0, -1*cot30, cot30, 0, -1*cot30, cot30};

    Point refPoint(H[e.first]->point().x() - tan30, H[e.first]->point().y() + 1);

    double theta = get_angle<bhs2017::K>(refPoint, H[e.first]->point(), H[e.second]->point());

    size_t cone = ( theta / alpha );

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

void addIncident(const vector<Vertex_handle>& H, const vector<pair<size_t, size_t>> L, const double alpha, vector<pair<size_t, size_t>> &E_A){

    double tan30 = tan(PI/6);

    for(auto i = L.begin(); i != L.end(); ++i){

        Point refPoint(H[i->first]->point().x() - tan30, H[i->first]->point().y() + 1);

        double theta = get_angle<bhs2017::K>(refPoint, H[i->first]->point(), H[i->second]->point());

        size_t cone = ( theta / alpha );

        cout << cone << "\n";
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
    vector<pair<size_t,size_t> > L;

    for( auto e=DT.finite_edges_begin(); e!=DT.finite_edges_end(); ++e ) {
        L.emplace_back( e->first->vertex( (e->second+1)%3 )->info(),
                        e->first->vertex( (e->second+2)%3 )->info() );
    }

    sort( L.begin(), L.end(), [&]( const auto& lhs, const auto& rhs ) {
        return bisectorLength( handles, lhs, alpha ) < bisectorLength( handles, rhs, alpha );
    });

    vector<pair<size_t,size_t> > E_A;

    addIncident(handles, L, alpha, E_A);



    // START PRINTER NONSENSE
    if( printLog ) {
        GraphPrinter printer(1);
        GraphPrinter::OptionsList options;

        options = {
            { "color", printer.inactiveEdgeColor },
            { "line width", to_string(printer.inactiveEdgeWidth) }
        };
        printer.drawEdges( DT, options );

//        options = { // active edge options
//            { "color", printer.activeEdgeColor },
//            { "line width", to_string(printer.activeEdgeWidth) }
//        };
//        printer.drawEdges( edgeList.begin(), edgeList.end(), options );


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
