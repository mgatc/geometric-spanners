#ifndef GSNUNF_BCC2012_H
#define GSNUNF_BCC2012_H

#include <bitset>
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

//typedef pair<size_t,size_t>                                         size_tPair;
//typedef boost::hash<size_tPair>                                     size_tPairHash;
//typedef unordered_map<size_tPair,bool,size_tPairHash>               size_tPairMap;
//
//bool selectEdge( const Delaunay& T, size_tPairMap &E, const Vertex_handle i, const Vertex_handle j, const size_t n, bool printLog = false ) {
//    assert( T.is_edge( i, j ) );
//    //if( printLog ) cout<<"add:("<<i->info()<<","<<j->info()<<") ";
//
//    auto existing = E.begin();
//    bool inserted = false;
//    tie(existing,inserted) = E.try_emplace( makeNormalizedPair( i->info(), j->info() ), false );
//    if(!inserted) existing->second = true;
//
//    return inserted;
//}

inline K::FT edgeLength( const vector<Vertex_handle>& H, const pair<size_t,size_t>& e ) {
    return distance( H[e.first]->point(), H[e.second]->point() );
}

template< typename T >
inline void fillCone( vector<T>& filled, vector<T>& partial, const size_t p, const size_t cone, bool printLog = false ) {
    const size_t NUM_CONES = filled.front().size();
    // This cone may have been partially covered by edges falling on the
    // boundaries at cone_p and cone_p+1. Because these cones have been
    // filled, if they are also partial, we need to mark the other partially
    // filled cone as filled and remove the partial status. This could
    // happen more than once, so do it in a loop
    filled.at(p)[cone] = true;

    if(printLog) cout<<"    fillCone["<<p<<"]["<<cone<<"]\n";

    size_t nextCone = (cone+1)%NUM_CONES;
    while( partial.at(p)[nextCone] ) {
        if(printLog) cout<<"    partial at "<<nextCone<<",fillCone["<<p<<"]["<<nextCone<<"]\n";
        filled.at(p)[nextCone] = true;
        partial.at(p)[nextCone] = false;
        nextCone = (nextCone+1)%NUM_CONES;
    }
    size_t prevCone = cone;
    while( partial.at(p)[prevCone] ) {
        if(printLog) cout<<"    partial at "<<prevCone;
        prevCone = (prevCone-1+NUM_CONES)%NUM_CONES;
        if(printLog) cout<<",fillCone["<<p<<"]["<<prevCone<<"]\n";
        filled.at(p)[prevCone] = true;
        partial.at(p)[prevCone] = false;
    }
}

} // namespace bcc2012

template< typename RandomAccessIterator, typename OutputIterator >
void BCC2012_7( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false ) {
    using namespace bcc2012;

    const size_t NUM_CONES = 8;
    const double alpha = 2*PI / NUM_CONES;

    if(printLog) cout<<"\nnumCones:"<<NUM_CONES<<"\n";
    if(printLog) cout<<"alpha:"<<alpha<<"\n";

    // Construct Delaunay triangulation
    bcc2012::Delaunay DT( pointsBegin, pointsEnd );
    size_t n = DT.number_of_vertices();
    if( n > SIZE_T_MAX - 1 ) return;
    if(printLog) cout<<"n:"<<n<<"\n";

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

    /*  Store the ID of the closest vertex in closest[vertex_id] = closest_id.

        Track cones that are filled and cones that are "partially" filled.
        Partially filled cones are cones whose clockwise-most boundary has an edge,
        but neither the indicated cone nor cone (cone+1)%8 have been filled.
    */
    vector<size_t> closest( n, SIZE_T_MAX );
    vector<bitset<NUM_CONES>> filled(n);
    vector<bitset<NUM_CONES>> partial(n);
    vector<pair<size_t,size_t>> E; // output edge list
    vector<pair<size_t,size_t>> E_star; // edges added from "Wedge"

    for( auto pq : L ) {
        size_t p = pq.first,
               q = pq.second;
        if(printLog) cout<<"p-q:"<<p<<" - "<<q<<"\n";

        if( filled.at(p).count() == NUM_CONES || filled.at(q).count() == NUM_CONES )
            continue;

        if(printLog) cout<<"  p_filled:"<<filled.at(p)<<"\n";
        if(printLog) cout<<"  p_partial:"<<partial.at(p)<<"\n";
        if(printLog) cout<<"  q_filled:"<<filled.at(q)<<"\n";
        if(printLog) cout<<"  q_partial:"<<partial.at(q)<<"\n";



        // Politely ask p if it wants an edge to q
        if( closest.at(p) == SIZE_T_MAX ) // Set closest
            closest.at(p) = q;

        // Determine angle of pq compared to closest of p
        double theta_p = get_angle<bcc2012::K>(
            handles.at(closest.at(p))->point(),
            handles.at(p)->point(),
            handles.at(q)->point()
        );
        size_t cone_p = size_t( theta_p / alpha ),
               cone_pPrev = (cone_p-1+NUM_CONES)%NUM_CONES;
        bool qOnBoundary = theta_p - cone_p*alpha < EPSILON,
             pGivenConeFilled = filled.at(p)[cone_p],
             pPrevConeFilled = filled.at(p)[cone_pPrev],
             pAbides = (!qOnBoundary && !pGivenConeFilled)
                    || ( qOnBoundary &&(!pGivenConeFilled || !pPrevConeFilled) );

        if(printLog) cout<<"  cone_p:"<<cone_p<<"\n";
        if(printLog && qOnBoundary) cout<<"  qOnBoundary\n";
        if(printLog && pAbides ) cout<<"  p abides!\n";



        // Politely ask q if it wants an edge to p
        if( closest.at(q) == SIZE_T_MAX ) // Set closest
            closest.at(q) = p;

        // Determine angle of pq compared to closest of q
        double theta_q = get_angle<bcc2012::K>(
            handles.at(closest.at(q))->point(),
            handles.at(q)->point(),
            handles.at(p)->point()
        );
        size_t cone_q = size_t( theta_q / alpha ),
               cone_qPrev = (cone_q-1+NUM_CONES)%NUM_CONES;
        bool pOnBoundary = theta_q - cone_q*alpha < EPSILON,
             qGivenConeFilled = filled.at(q)[cone_q],
             qPrevConeFilled = filled.at(q)[cone_qPrev],
             qAbides = (!pOnBoundary && !qGivenConeFilled)
                    || ( pOnBoundary &&(!qGivenConeFilled || !qPrevConeFilled) );

        if(printLog) cout<<"  cone_q:"<<cone_q<<"\n";
        if(printLog && pOnBoundary) cout<<"  pOnBoundary\n";
        if(printLog && qAbides ) cout<<"  q abides!\n";



        // Only continue if p and q both consent to add the edge
        // Add the edge then do some bookkeeping to track filled and partially filled cones
        if( pAbides && qAbides ) {
            E.emplace_back(p,q);
            // Bookkeeping for p
            if( qOnBoundary ) {
                // check the fill status of the given cone and the (cone-1+NUM_CONES)%NUM_CONES cone
                // if both are empty, mark cone as partially filled
                if( !pGivenConeFilled && !pPrevConeFilled ) {
                    partial.at(p)[cone_p] = true;
                    if(printLog) cout<<"  onBoundary: Both cones empty! Partial added"<<cone_p<<"\n";
                }
                // if one is empty, mark the empty as filled
                else if( pGivenConeFilled ^ pPrevConeFilled ) { // one is filled
                    if(printLog) cout<<"  onBoundary: One cone empty!\n";
                    // Mark both as filled
                    fillCone( filled, partial, p, cone_p, printLog );
                    fillCone( filled, partial, p, cone_pPrev, printLog );
                }
            } else {
                fillCone( filled, partial, p, cone_p, printLog );
            }
            // Bookkeeping for q
            if( pOnBoundary ) {
                if( !qGivenConeFilled && !qPrevConeFilled ) {
                    partial.at(q)[cone_q] = true;
                    if(printLog) cout<<"  onBoundary: Both cones empty! Partial added"<<cone_q<<"\n";
                }
                // if one is empty, mark the empty as filled
                else if( qGivenConeFilled ^ qPrevConeFilled ) { // one is filled
                    if(printLog) cout<<"  onBoundary: One cone empty!\n";
                    // Mark both as filled
                    fillCone( filled, partial, q, cone_q, printLog );
                    fillCone( filled, partial, q, cone_qPrev, printLog );
                }
            } else {
                fillCone( filled, partial, q, cone_q, printLog );
            }
            // Call Wedge

        }
    }




    // Edge list is only needed for printing. Remove for production.
    vector< pair<Point,Point> > edgeList;
    edgeList.reserve( E.size() );

    // Send resultant graph to output iterator
    for( auto e : E ) {
        // Edge list is only needed for printing. Remove for production.
        edgeList.emplace_back( handles.at(e.first)->point(), handles.at(e.second)->point() );

        *result = make_pair( handles.at(e.first)->point(), handles.at(e.second)->point() );
        ++result;
        *result = make_pair( handles.at(e.second)->point(), handles.at(e.first)->point() );
        ++result;
    }


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

