#ifndef GSNUNF_BCC2012_H
#define GSNUNF_BCC2012_H

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
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

//#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "metrics.h"
#include "utilities.h"


namespace gsnunf {

using namespace std;

namespace bcc2012 {

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K>      Vb;
typedef CGAL::Triangulation_face_base_2<K>                          Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef CGAL::Aff_transformation_2<K>                               Transformation;
typedef Delaunay::Vertex_handle                                     Vertex_handle;
typedef Delaunay::Vertex_circulator                                 Vertex_circulator;
typedef CGAL::Vector_2<K>                                           Vector_2;
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

inline K::FT edgeLength( const vector<Vertex_handle>& H, const pair<size_t,size_t>& e ) {
    return distance( H[e.first]->point(), H[e.second]->point() );
}

} // namespace bcc2012

template< typename RandomAccessIterator, typename OutputIterator >
void BCC2012_7( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false ) {
    using namespace bcc2012;

    const double alpha = PI / 4;

    //if(printLog) cout<<"alpha:"<<alpha<<",";

    // Construct Delaunay triangulation
    bcc2012::Delaunay DT( pointsBegin, pointsEnd );
    size_t n = DT.number_of_vertices();
    if( n > SIZE_T_MAX - 1 ) return;

    vector<bcc2012::Vertex_handle> handles(n);

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
        return edgeLength( handles, lhs ) < edgeLength( handles, rhs );
    });
    vector<size_t> closest( n, SIZE_T_MAX );

    for( auto e : L ) {
        // Set closest
        if( closest[e.first] == SIZE_T_MAX )
            closest[e.first] = e.second;
        if( closest[e.second] == SIZE_T_MAX )
            closest[e.second] = e.first;

    }

    for(size_t i=0;i<n;i++)
        cout<<i<<": "<<closest[i]<<"\n";















    // Done. Send edges from G_prime with value == true (selected by both endpoints) to output.

    // Edge list is only needed for printing. Remove for production.
//    vector< pair<Point,Point> > edgeList;
//    edgeList.reserve( G_prime.size() );

    // Send resultant graph to output iterator
//    for( auto e : G_prime ) {
//        if( e.second ) { // e.second holds the bool value of whether both vertices of an edge selected the edge
//            // Edge list is only needed for printing. Remove for production.
//            //edgeList.emplace_back( handles.at(e.first.first)->point(), handles.at(e.first.second)->point() );
//
//            *result = make_pair( handles.at(e.first.first)->point(), handles.at(e.first.second)->point() );
//            ++result;
//            *result = make_pair( handles.at(e.first.second)->point(), handles.at(e.first.first)->point() );
//            ++result;
//        }
//    }


    //
    //
    // START PRINTER NONSENSE
    //
    //

//    if( printLog ) {
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

        printer.print( "bcc2012" );
        cout<<"\n";
//    }

    //
    //
    // END PRINTER NONSENSE
    //
    //

} // function BCC2012_7

} // namespace gsnunf

#endif // GSNUNF_BCC2012_H

