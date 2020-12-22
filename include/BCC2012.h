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

inline size_t getCone( const vector<Vertex_handle>& handles, const vector<size_t>& closest, const size_t& p, const size_t& q, const double& alpha ) {
    return size_t( get_angle<bcc2012::K>(
            handles.at(closest.at(p))->point(),
            handles.at(p)->point(),
            handles.at(q)->point() )
        / alpha
    );
}



} // namespace bcc2012

template< typename RandomAccessIterator, typename OutputIterator >
void BCC2012_7( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false ) {
    using namespace bcc2012;

    const size_t NUM_CONES = 8;
    const double alpha = 2*PI / NUM_CONES;
    const double PI_OVER_TWO = PI / 2;

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

    vector<size_t> closest( n, SIZE_T_MAX );
    vector<bitset<NUM_CONES>> filled(n);
    vector<pair<size_t,size_t>> E; // output edge list
    vector<pair<size_t,size_t>> E_star; // edges added from "Wedge"

    for( auto pq : L ) {
        size_t p = pq.first,
               q = pq.second;
        if(printLog) cout<<"p-q:"<<p<<" - "<<q<<"\n";
        if(printLog) cout<<"p-q:"<<handles.at(p)->point()<<" - "<<handles.at(q)->point()<<"\n";

        if(printLog) cout<<"  p_filled:"<<filled.at(p)<<"\n";
        if(printLog) cout<<"  q_filled:"<<filled.at(q)<<"\n";


        if( filled.at(p).count() == NUM_CONES || filled.at(q).count() == NUM_CONES )
            continue;

        // Politely ask p if it wants an edge to q
        if( closest.at(p) == SIZE_T_MAX ) // Set closest
            closest.at(p) = q;

        size_t cone_p = getCone( handles, closest, p, q, alpha ),
               cone_pPrev = (cone_p-1+NUM_CONES)%NUM_CONES;

        double theta_p = get_angle<bcc2012::K>(
            handles.at(closest.at(p))->point(),
            handles.at(p)->point(),
            handles.at(q)->point()
        );

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

        size_t cone_q = getCone( handles, closest, q, p, alpha ),
               cone_qPrev = (cone_q-1+NUM_CONES)%NUM_CONES;

        double theta_q = get_angle<bcc2012::K>(
            handles.at(closest.at(q))->point(),
            handles.at(q)->point(),
            handles.at(p)->point()
        );

        bool pOnBoundary = theta_q - cone_q*alpha < EPSILON,
             qGivenConeFilled = filled.at(q)[cone_q],
             qPrevConeFilled = filled.at(q)[cone_qPrev],
             qAbides = (!pOnBoundary && !qGivenConeFilled)
                    || ( pOnBoundary &&(!qGivenConeFilled || !qPrevConeFilled) );

        if(printLog) cout<<"  cone_q:"<<cone_q<<"\n";
        if(printLog && pOnBoundary) cout<<"  pOnBoundary\n";
        if(printLog && qAbides ) cout<<"  q abides!\n";



        // Only continue if p and q both consent to add the edge
        if( pAbides && qAbides ) {
            E.emplace_back(p,q);

            // Wedge on each cone of pq and qp
            // There will be at least one for each, but there could
            // be two cones for one or both pq and qp if the edge
            // falls on the boundary of a cone and the cone is not already filled
            vector<vector<size_t>> wedge;
            // Bookkeeping for p
            if( qOnBoundary ) {
                if(!filled.at(p)[cone_pPrev])
                    wedge.push_back({p,q,cone_pPrev});
                filled.at(p)[cone_pPrev] = true;
            }
            if(!filled.at(p)[cone_p])
                wedge.push_back({p,q,cone_p});
            filled.at(p)[cone_p] = true;

            // Bookkeeping for q
            if( pOnBoundary ) {
                if(!filled.at(q)[cone_qPrev])
                    wedge.push_back({q,p,cone_qPrev});
                filled.at(q)[cone_qPrev] = true;
            }
            if(!filled.at(q)[cone_q])
                wedge.push_back({q,p,cone_q});
            filled.at(q)[cone_q] = true;


            // Wedge on p, q
            for( auto params : wedge ) {
                // p is params.at(0)
                // q is params.at(1)
                // cone is params.at(2)

                // q_m[i] holds the circulator for q_{m-i}
                vector<Vertex_circulator> q_m(3);
                // find q
                q_m[2] = DT.incident_vertices( handles.at(params.at(0)) );
                while( ++q_m[2] != handles.at(params.at(1)) ); // point to q
                const auto q_i = q_m[2];

                // Process edges from q_j to q_i
                // set and increment q_m and q_{m-1} so the sequence will be correct at the start of the loop
                for( size_t i=0; i<q_m.size()-1; ++i ) {
                    q_m[i] = q_i;
                    ++q_m[i];
                }
                // Get the first and last vertex in cone_p, called q_j and q_k
                // Rotate q_j CCW until we leave the cone
                while( !DT.is_infinite(++q_m[0])
                && getCone(handles,closest,params.at(0),q_m[0]->info(),alpha) == params.at(2) ) {
                    if( q_m[2] != q_i
                    ||( q_m[2] == q_i && get_angle<K>(handles.at(params.at(0))->point(), q_i->point(), q_m[1]->point()) > PI_OVER_TWO ) ) {
                        E_star.emplace_back( q_m[1]->info(), q_m[2]->info() );
                        if(printLog) cout<<"  wedge["<<params.at(0)<<"][cone"<<params.at(2)<<"]: "<<q_m[1]->info()<< " - " <<q_m[2]->info()<<"\n";
                    }
                    q_m[2] = q_m[1];
                    q_m[1] = q_m[0];
                };

                // Process edges from q_k to q_i
                for( size_t i=0; i<q_m.size()-1; ++i ) {
                    q_m[i] = q_i;
                    --q_m[i];
                }
                q_m[2] = q_i;
                while( !DT.is_infinite(--q_m[0])
                && getCone(handles,closest,params.at(0),q_m[0]->info(),alpha) == params.at(2) ) {
                    if( q_m[2] != q_i
                    ||( q_m[2] == q_i && get_angle<K>(q_m[1]->point(), q_i->point(), handles.at(params.at(0))->point()) > PI_OVER_TWO ) ) {
                        E_star.emplace_back( q_m[1]->info(), q_m[2]->info() );
                        if(printLog)cout<<"  wedge["<<params.at(0)<<"][cone"<<params.at(2)<<"]: "<<q_m[1]->info()<< " - " <<q_m[2]->info()<<"\n";
                    }
                    q_m[2] = q_m[1];
                    q_m[1] = q_m[0];
                };
            }
        }
    }

    // Combine E and E_star, remove duplicates
    E.insert( E.end(), E_star.begin(), E_star.end() );
    sort( E.begin(), E.end() );
    E.erase( unique(E.begin(), E.end(), []( const auto& l, const auto& r) {
        return ( l.first == r.first && l.second == r.second )
            || ( l.first == r.second && l.second == r.first );
    }), E.end() );


    // Edge list is only needed for printing. Remove for production.
    vector< pair<Point,Point> > edgeList;
    edgeList.reserve( E.size()+E_star.size() );

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

    if( printLog && n <= 200 ) {
        GraphPrinter printer(0.006);
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

        printer.print( "bcc2012a" );
        cout<<"\n";
    }

    //
    //
    // END PRINTER NONSENSE
    //
    //

} // function BCC2012_7

} // namespace gsnunf

#endif // GSNUNF_BCC2012_H

