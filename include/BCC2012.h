#ifndef GSNUNF_BCC2012_H
#define GSNUNF_BCC2012_H

#include <bitset>
#include <cmath>         // ceil, floor, isinf
#include <functional>
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

inline size_t getPreviousCone( const size_t& cone, const size_t& numCones ) {
    return (cone-1+numCones)%numCones;
}

inline bool vertexAgreesOnEdge( const vector<Vertex_handle>& handles,
                                vector<size_t>& closest,
                                const vector<bitset<8>>& filled,
                                const size_t p,
                                const size_t q,
                                const double alpha,
                                const size_t numCones,
                                size_t& cone,
                                size_t& conePrev,
                                bool& qOnBoundary ) {

    if( closest.at(p) == SIZE_T_MAX ) // First, make sure the closest vertex is set
        closest.at(p) = q;

    cone = getCone( handles, closest, p, q, alpha );
    conePrev = getPreviousCone( cone, numCones );

    double theta = get_angle<bcc2012::K>(
        handles.at(closest.at(p))->point(),
        handles.at(p)->point(),
        handles.at(q)->point()
    );

    qOnBoundary = theta - cone*alpha < EPSILON;

    bool pGivenConeFilled = filled.at(p)[cone],
         pPrevConeFilled = filled.at(p)[conePrev];

    return (!qOnBoundary && !pGivenConeFilled)
        || ( qOnBoundary &&(!pGivenConeFilled || !pPrevConeFilled) );
}

inline void updateVertexConeStatus( vector<bitset<8>>& filled, vector<vector<size_t>>& wedge,
                                    const size_t p, const size_t q,
                                    const bool qOnBoundary, const size_t cone, const size_t conePrev ) {
    if( qOnBoundary ) {
        if(!filled.at(p)[conePrev])
            wedge.push_back({p,q,conePrev});
        filled.at(p)[conePrev] = true;
    }
    if(!filled.at(p)[cone])
        wedge.push_back({p,q,cone});
    filled.at(p)[cone] = true;
}

inline void wedge( const Delaunay& DT, const vector<Vertex_handle>& handles, const vector<size_t>& closest, const vector<size_t>& params,
                   vector<pair<size_t,size_t>>& addToE_star, const Vertex_circulator& q_i, const double alpha ) {
    // params
    // p is params.at(0)
    // q is params.at(1)
    // cone is params.at(2)

    // q_m[i] holds the circulator for q_{m-i}
    vector<Vertex_circulator> q_m(3);

    // Setup function objects for wedge
    vector<std::function<Vertex_circulator(Vertex_circulator&)>> step;
    step.push_back( [] ( Vertex_circulator& c ) { return ++c; } );
    step.push_back( [] ( Vertex_circulator& c ) { return --c; } );

    for( size_t i=0; i<step.size(); i++ ) {
        fill( q_m.begin(), q_m.end(), q_i );

        // set and increment q_m and q_{m-1} so the sequence will be correct at the start of the loop
        for( size_t j=0; j<q_m.size()-1; ++j )
            step[i](q_m[j]);

        size_t a = params.at(0);
        //     b = q_i
        size_t c = params.at(0);

        if(i==0)
            c = q_m[1]->info();
        else
            a = q_m[1]->info();

        // Get the first and last vertex in cone_p, called q_j and q_k
        // Rotate q_j CCW until we leave the cone
        while( !DT.is_infinite(q_m[0])
        && getCone(handles,closest,params.at(0),q_m[0]->info(),alpha) == params.at(2)
        && !DT.is_infinite(step[i](q_m[0]))
        && getCone(handles,closest,params.at(0),q_m[0]->info(),alpha) == params.at(2) ) {
            if( q_m[2] != q_i
            ||( q_m[2] == q_i && get_angle<K>(handles.at(a)->point(), q_i->point(), handles.at(c)->point()) > PI_OVER_TWO ) ) {
                addToE_star.emplace_back( q_m[1]->info(), q_m[2]->info() );
            }
            q_m[2] = q_m[1];
            q_m[1] = q_m[0];
        };
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

    vector<size_t> closest( n, SIZE_T_MAX ); // id of closest vertex for orienting cones
    vector<bitset<NUM_CONES>> filled(n); // status of each cone
    vector<pair<size_t,size_t>> E; // output edge list
    vector<pair<size_t,size_t>> E_star; // edges added from "Wedge"

    for( auto pq : L ) {
        size_t p = pq.first,
               q = pq.second;

        if(printLog) cout<<"p-q:"<<p<<" - "<<q<<"\n";
        if(printLog) cout<<"p-q:"<<handles.at(p)->point()<<" - "<<handles.at(q)->point()<<"\n";
        if(printLog) cout<<"  p_filled:"<<filled.at(p)<<"\n";
        if(printLog) cout<<"  q_filled:"<<filled.at(q)<<"\n";

        // If either p or q's cone is filled, don't even bother
        if( filled.at(p).count() == NUM_CONES || filled.at(q).count() == NUM_CONES )
            continue;

        // Politely ask p if it wants an edge to q
        size_t cone_p = 0,
           cone_pPrev = 0;
        bool qOnBoundary = false;
        bool pAbides = vertexAgreesOnEdge( handles, closest, filled, p, q, alpha, NUM_CONES,
                                           cone_p, cone_pPrev, qOnBoundary );

        if(printLog) cout<<"  cone_p:"<<cone_p<<"\n";
        if(printLog && qOnBoundary) cout<<"  qOnBoundary\n";
        if(printLog && pAbides ) cout<<"  p abides!\n";


        // Politely ask q if it wants an edge to p
        size_t cone_q = 0,
           cone_qPrev = 0;
        bool pOnBoundary = false;
        bool qAbides = vertexAgreesOnEdge( handles, closest, filled, q, p, alpha, NUM_CONES,
                                           cone_q, cone_qPrev, pOnBoundary );

        if(printLog) cout<<"  cone_q:"<<cone_q<<"\n";
        if(printLog && pOnBoundary) cout<<"  pOnBoundary\n";
        if(printLog && qAbides ) cout<<"  q abides!\n";

        // Only continue if p and q both consent to add the edge
        if( pAbides && qAbides ) {
            E.emplace_back(p,q); // Place the edge

            // Wedge on each cone of pq and qp
            // There will be at least one for each, but there could
            // be two cones for one or both pq and qp if the edge
            // falls on the boundary of a cone and the cone is not already filled
            vector<vector<size_t>> W; // holds the parameters for each call to wedge

            // Bookkeeping for p
            updateVertexConeStatus( filled, W, p, q, qOnBoundary, cone_p, cone_pPrev );

            // Bookkeeping for q
            updateVertexConeStatus( filled, W, q, p, pOnBoundary, cone_q, cone_qPrev );

            vector<pair<size_t,size_t>> addToE_star;

            // Wedge on p, q
            for( auto params : W ) {
                // find q
                auto q_z = DT.incident_vertices( handles.at(params.at(0)) );
                while( ++q_z != handles.at(params.at(1)) ); // point to q
                const auto q_i(q_z);

                wedge( DT, handles, closest, params, addToE_star, q_i, alpha );
            }
            E_star.insert( E_star.end(), addToE_star.begin(), addToE_star.end() );
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

