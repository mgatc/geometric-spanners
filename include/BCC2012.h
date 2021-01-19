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


enum Q_primePosition  { between_j_i=0, between_i_k=1, not_set=2 };

inline K::FT edgeLength( const vector<Vertex_handle>& H, const pair<size_t,size_t>& e ) {
    return distance( H[e.first]->point(), H[e.second]->point() );
}

inline size_t getCone( const vector<Vertex_handle>& handles,
                       const vector<size_t>& closest,
                       const size_t p,
                       const size_t q,
                       const double alpha ) {
    return size_t( get_angle<bcc2012::K>(
            handles.at(closest.at(p))->point(),
            handles.at(p)->point(),
            handles.at(q)->point() )
        / alpha
    );
}

inline size_t getPreviousCone( const size_t cone, const size_t numCones ) {
    return (cone-1+numCones)%numCones;
}

template< size_t DEGREE, size_t NUM_CONES = DEGREE+1  >
inline bool vertexAgreesOnEdge( const vector<Vertex_handle>& handles,
                                vector<size_t>& closest,
                                const vector<bitset<NUM_CONES>>& filled,
                                const size_t p,
                                const size_t q,
                                size_t& cone,
                                size_t& conePrev,
                                bool& qOnBoundary ) {

    if( closest.at(p) == SIZE_T_MAX ) { // First, make sure the closest vertex is set
        closest.at(p) = q;
        qOnBoundary = true; // the only vertex that falls on a boundary should be the closest
    }
    const double ALPHA = 2*PI/NUM_CONES;
    cone = getCone( handles, closest, p, q, ALPHA );
    conePrev = getPreviousCone( cone, NUM_CONES );

    double theta = get_angle<bcc2012::K>(
        handles.at(closest.at(p))->point(),
        handles.at(p)->point(),
        handles.at(q)->point()
    );

    bool pGivenConeFilled = filled.at(p)[cone],
         pPrevConeFilled = filled.at(p)[conePrev];

    return (!qOnBoundary && !pGivenConeFilled)
        || ( qOnBoundary &&(!pGivenConeFilled || !pPrevConeFilled) );
}

template< size_t DEGREE, size_t NUM_CONES = DEGREE+1 >
inline void updateVertexConeStatus( vector<bitset<NUM_CONES>>& filled,
                                    vector<vector<size_t>>& wedge,
                                    const size_t p,
                                    const size_t q,
                                    const bool qOnBoundary,
                                    const size_t cone,
                                    const size_t conePrev ) {
    if( qOnBoundary ) {
        if(!filled.at(p)[conePrev])
            wedge.push_back({p,q,conePrev});
        filled.at(p)[conePrev] = true;
    }
    if(!filled.at(p)[cone])
        wedge.push_back({p,q,cone});
    filled.at(p)[cone] = true;
}

inline double get_min_angle( const Vertex_handle& p,
                             const Vertex_handle& q,
                             const Vertex_handle& r ) {
    return CGAL::min(
        get_angle<K>( p->point(), q->point(), r->point() ),
        get_angle<K>( r->point(), q->point(), p->point() )
    );
}

/*
 *  Performs algorithm "wedge" with the exception of line 1. The
 *  for loop must be accomplished outside of the function.
 *  Degree 6 and 7 wedge algorithms are implemented as templates.
 *  This is the primary template, but will never be used. Instead,
 *  the specialized templates are used, defined below.
 */
template<size_t DEGREE, size_t NUM_CONES=DEGREE+1>
inline void wedge( const Delaunay& DT,
                   const vector<Vertex_handle>& handles,
                   const vector<size_t>& closest,
                   const vector<size_t>& params,
                   vector<pair<size_t,size_t>>& addToE_star,
                   const Vertex_circulator& q_i,
                   const bool printLog ) {
    assert(DEGREE==6||DEGREE==7);
}

template<>
inline void wedge<7>( const Delaunay& DT,
                      const vector<Vertex_handle>& handles,
                      const vector<size_t>& closest,
                      const vector<size_t>& params,
                      vector<pair<size_t,size_t>>& addToE_star,
                      const Vertex_circulator& q_i,
                      const bool printLog ) {

    const double ALPHA = 2*PI/8;

    const size_t& p = params.at(0);
    const size_t& q = params.at(1);
    const size_t& cone = params.at(2);

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

        size_t a = p;
        //     b = q_i
        size_t c = p;

        if(i==0)
            c = q_m[1]->info();
        else
            a = q_m[1]->info();

        // Get the first and last vertex in cone_p, called q_j and q_k
        // Rotate q_j CCW until we leave the cone
        while( !DT.is_infinite(q_m[0])
        && getCone( handles, closest, p, q_m[0]->info(), ALPHA ) == cone
        && !DT.is_infinite(step[i](q_m[0]))
        && getCone( handles, closest, p, q_m[0]->info(), ALPHA ) == cone ) {
            if( q_m[2] != q_i
            ||( q_m[2] == q_i && get_angle<K>(handles.at(a)->point(), q_i->point(), handles.at(c)->point()) > PI_OVER_TWO ) ) {
                addToE_star.emplace_back( q_m[1]->info(), q_m[2]->info() );
            }
            step[i](q_m[2]);
            step[i](q_m[1]);
        };
    }
}

template<>
inline void wedge<6>( const Delaunay& DT,
                      const vector<Vertex_handle>& handles,
                      const vector<size_t>& closest,
                      const vector<size_t>& params,
                      vector<pair<size_t,size_t>>& addToE_star,
                      const Vertex_circulator& q_i,
                      const bool printLog ) {

    const double ALPHA = 2*PI/7;

    const size_t& p = params.at(0);
    const size_t& q = params.at(1);
    const size_t& cone = params.at(2);

    // Line 2: Build Q
    auto N_p = DT.incident_vertices(handles[p]);

    while( ++N_p != handles[q] ); // point to q aka q_i

    while( !DT.is_infinite(++N_p) // move CCW until we leave the cone
        && ( getCone(handles, closest, p, N_p->info(), ALPHA) == cone
            || N_p->info() == q ) );

    vector<size_t> Q; // the ordered neighbors in the current cone

    while( !DT.is_infinite(--N_p) // move CW until we leave the cone, adding each to Q
        && ( getCone(handles, closest, p, N_p->info(), ALPHA) == cone
            || N_p->info() == q ) ) { // handles the case where q is on a boundary
        Q.push_back(N_p->info());
    }
    assert( !Q.empty() );

    // Line 3: Build Q'
    unordered_set<size_t> Q_prime; // select elements from Q (line 3)

    size_t j = 0, // the first point in the ordered list of neighbors is q_j
           i = 0, // the opposite end of the edge
           k = Q.size()-1; // the last point in the ordered list is q_k

    // find the index of q_i in Q
    for( int n=0; n<Q.size(); ++n )
        if( Q[n] == q )
            i = n;


    Q_primePosition Q_primePos = not_set;

    for( int n=j+1; n<k; ++n ) {
        if( n != i // q_i
        && get_angle<K>( handles.at(Q.at(n+1))->point(),
                         handles.at(Q.at(n))->point(),
                         handles.at(Q.at(n-1))->point() ) < SIX_PI_OVER_SEVEN ) {
            Q_prime.insert(Q[n]);
            Q_primePos = static_cast<Q_primePosition>(n>i); // if we have passed i, n>i will be true == 1 == between_i_k
        }
        if(printLog)cout<<get_angle<K>( handles.at(Q.at(n+1))->point(),
                             handles.at(Q.at(n))->point(),
                             handles.at(Q.at(n-1))->point() )<<" <>? "<<SIX_PI_OVER_SEVEN<<"\n";
    }
    if(printLog) {
        cout<<"p:"<<p<<" q:"<<q<<"\n";
        cout<<"Q:";
        for( auto v : Q ) cout<<v<<" ";

        cout<<"\n";
        cout<<"Q':";
        for( auto v : Q_prime ) cout<<v<<" ";

        cout<<"\n";

        cout<<"j:"<<j<<" i:"<<i<<" k:"<<k<<endl;
        cout<<"j+1:"<<int(j+1)<<" i-2:"<<int(i-2)<<" k:"<<k<<endl;
    }

    // TODO:I don't trust the types here... i and k are size_t. Converting to int isn't good...
    // Line 4: Add select edges
    for( int n=j+1; n<int(i)-2; ++n )
        if( !contains( Q_prime, Q.at(n)) && !contains( Q_prime, Q.at(n+1)) )
            addToE_star.emplace_back( Q.at(n), Q.at(n+1) );

    if(printLog)cout<<"add to E_star:";
    if(printLog)for( auto e : addToE_star ) cout<<e.first<<" "<<e.second<<" - ";

    for( int n=i+1; n<int(k)-2; ++n )
        if( !contains( Q_prime, Q.at(n)) && !contains( Q_prime, Q.at(n+1)) )
            addToE_star.emplace_back( Q.at(n), Q.at(n+1) );

    if(printLog)cout<<"add to E_star:";
    if(printLog)for( auto e : addToE_star ) cout<<e.first<<" "<<e.second<<" - ";

    size_t f = i, // will hold the index in Q of the first point in Q_prime
           a = i;

    // Line 5:
    switch( Q_primePos ) {
    case between_i_k:
        // Line 6-7
        if( i != j
         && i-1 != j
         && get_angle<K>(handles.at(p)->point(), handles.at(Q.at(i))->point(), handles.at(Q.at(i-1))->point()) > FOUR_PI_OVER_SEVEN )
            addToE_star.emplace_back( Q.at(i), Q.at(i-1) );

        if(printLog&&i!=j&&i-1!=j)cout<<get_angle<K>(handles.at(p)->point(), handles.at(Q.at(i))->point(), handles.at(Q.at(i-1))->point())<<"\n";

        while( ++f < Q.size() && !contains( Q_prime, Q.at(f) ) );

        a = f - int(f==Q.size())*Q.size(); // the first point after f not in Q_prime, avoid vector overflow
        while( ++a < Q.size() && contains(Q_prime, Q.at(a) ) );

        if(printLog)cout<<"f:"<<f<<" a:"<<a<<"\n";

        if( f == i+1 ) {
            // TODO: Need to check for equality in this one too
            if( a != k && get_angle<K>( handles.at(Q.at(i+1))->point(),
                                         handles.at(Q.at(i))->point(),
                                         handles.at(p)->point() ) < FOUR_PI_OVER_SEVEN ) {
                addToE_star.emplace_back(Q.at(f), Q.at(a));
            }
            if(printLog)cout<<get_angle<K>( handles.at(Q.at(i+1))->point(),
                                         handles.at(Q.at(i))->point(),
                                         handles.at(p)->point() )<<"\n";
            if( i != j && i !=k && f+1 != k && get_angle<K>( handles.at(Q.at(i+1))->point(),
                                           handles.at(Q.at(i))->point(),
                                           handles.at(p)->point() ) > FOUR_PI_OVER_SEVEN ) {
                addToE_star.emplace_back(Q.at(i), Q.at(f+1));
            }
        } else {
            size_t l = Q.size(); // the last point in Q_prime
            while( --l > 0 && !contains( Q_prime, Q.at(l) ) );
            if( l > 0 ) {
                size_t b = l; // set to max n that is less than l in Q but not Q_prime
                while( --b > 0 && contains( Q_prime, Q.at(b) ) );

                if(printLog)cout<<"l:"<<l<<" b:"<<b<<endl;
                if( l == k-1 ) {
                    addToE_star.emplace_back( Q.at(l), Q.at(b) );
                } else {
                    addToE_star.emplace_back( Q.at(b), Q.at(l+1) );
                    if( contains( Q_prime, Q.at(l-1) ) )
                        addToE_star.emplace_back( Q.at(l), Q.at(l-1) );
                }
            }
        }
        break;
    case between_j_i:
        // Line 6-7
        if( i != k
         && i+1 != k
         && get_angle<K>(handles.at(Q.at(i+1))->point(), handles.at(Q.at(i))->point(), handles.at(p)->point()) > FOUR_PI_OVER_SEVEN )
        {
            addToE_star.emplace_back( Q.at(i), Q.at(i+1) );
        }

        if(printLog)if(i!=k&&i+1!=k)cout<<get_angle<K>(handles.at(Q.at(i+1))->point(), handles.at(Q.at(i))->point(), handles.at(p)->point())<<"\n";

        while( --f > 0 && !contains( Q_prime, Q.at(f) ) );

        a = f + int(f==0)*Q.size(); // the first point after f not in Q_prime, protect against unsigned underflow
        while( --a > 0 && contains(Q_prime, Q.at(a) ) );

        if(printLog)cout<<"f:"<<f<<" a:"<<a<<"\n";

        if( f == i-1 ) {
            if(printLog)cout<<get_angle<K>( handles.at(p)->point(),
                                         handles.at(Q.at(i))->point(),
                                         handles.at(Q.at(i-1))->point() )<<"\n";
            // TODO: Need to check for equality in this one too
            if( a != j && get_angle<K>( handles.at(p)->point(),
                                         handles.at(Q.at(i))->point(),
                                         handles.at(Q.at(i-1))->point() ) < FOUR_PI_OVER_SEVEN ) {
                addToE_star.emplace_back(Q.at(f), Q.at(a));
            }
            if( i != j && i !=k && f-1 != j && get_angle<K>( handles.at(p)->point(),
                                           handles.at(Q.at(i))->point(),
                                           handles.at(Q.at(i-1))->point() ) > FOUR_PI_OVER_SEVEN ) {
                addToE_star.emplace_back(Q.at(i), Q.at(f-1));
            }
        } else {
            size_t l = 0; // the last point in Q_prime
            while( ++l < Q.size() && !contains( Q_prime, Q.at(l) ) );
            if( l > 0 ) {
                size_t b = l; // set to max n that is less than l in Q but not Q_prime
                while( ++b < Q.size() && contains( Q_prime, Q.at(b) ) );

//                cout<<"l:"<<l<<" b:"<<b<<endl;
                if( l == j+1 ) {
                    addToE_star.emplace_back( Q.at(l), Q.at(b) );
                } else {
                    addToE_star.emplace_back( Q.at(b), Q.at(l-1) );
                    if( contains( Q_prime, Q.at(l+1) ) )
                        addToE_star.emplace_back( Q.at(l), Q.at(l+1) );
                }
            }
        }
        break;
    }

    if(printLog)cout<<"add to E_star:";
    if(printLog)for( auto e : addToE_star ) cout<<e.first<<" "<<e.second<<" - ";

    if(printLog)cout<<"\n\n";
}


} // namespace bcc2012

template< size_t DEGREE = 7, size_t NUM_CONES = DEGREE+1,
          typename RandomAccessIterator, typename OutputIterator >
void BCC2012( RandomAccessIterator pointsBegin,
              RandomAccessIterator pointsEnd,
              OutputIterator result,
              bool printLog = false ) {
    using namespace bcc2012;

    assert( DEGREE == 7 || DEGREE == 6 );
    const double ALPHA = 2*PI/NUM_CONES;

//    if(printLog) cout<<"\nnumCones:"<<NUM_CONES<<"\n";
//    if(printLog) cout<<"ALPHA:"<<ALPHA<<"\n";

    // Construct Delaunay triangulation
    bcc2012::Delaunay DT( pointsBegin, pointsEnd );
    size_t n = DT.number_of_vertices();
    if( n > SIZE_T_MAX - 1 ) return;
//    if(printLog) cout<<"n:"<<n<<"\n";

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

//        if(printLog) cout<<"p-q:"<<p<<" - "<<q<<"\n";
//        if(printLog) cout<<"p-q:"<<handles.at(p)->point()<<" - "<<handles.at(q)->point()<<"\n";
//        if(printLog) cout<<"  p_filled:"<<filled.at(p)<<"\n";
//        if(printLog) cout<<"  q_filled:"<<filled.at(q)<<"\n";

        // If either p or q's cone is filled, don't even bother
        if( filled.at(p).count() == NUM_CONES || filled.at(q).count() == NUM_CONES )
            continue;

        // Politely ask p if it wants an edge to q
        size_t cone_p = 0,
           cone_pPrev = 0;
        bool qOnBoundary = false;
        bool pAbides = vertexAgreesOnEdge<DEGREE>( handles, closest, filled, p, q,
                                           cone_p, cone_pPrev, qOnBoundary );

//        if(printLog) cout<<"  cone_p:"<<cone_p<<"\n";
//        if(printLog && qOnBoundary) cout<<"  qOnBoundary\n";
//        if(printLog && pAbides ) cout<<"  p abides!\n";


        // Politely ask q if it wants an edge to p
        size_t cone_q = 0,
           cone_qPrev = 0;
        bool pOnBoundary = false;
        bool qAbides = vertexAgreesOnEdge<DEGREE>( handles, closest, filled, q, p,
                                           cone_q, cone_qPrev, pOnBoundary );

//        if(printLog) cout<<"  cone_q:"<<cone_q<<"\n";
//        if(printLog && pOnBoundary) cout<<"  pOnBoundary\n";
//        if(printLog && qAbides ) cout<<"  q abides!\n";

        // Only continue if p and q both consent to add the edge
        if( pAbides && qAbides ) {
            E.emplace_back(p,q); // Place the edge

            // Wedge on each cone of pq and qp
            // There will be at least one for each, but there could
            // be two cones for one or both pq and qp if the edge
            // falls on the boundary of a cone and the cone is not already filled
            vector<vector<size_t>> W; // holds the parameters for each call to wedge

            // Bookkeeping for p
            updateVertexConeStatus<DEGREE>( filled, W, p, q, qOnBoundary, cone_p, cone_pPrev );

            // Bookkeeping for q
            updateVertexConeStatus<DEGREE>( filled, W, q, p, pOnBoundary, cone_q, cone_qPrev );

            vector<pair<size_t,size_t>> addToE_star;

            // Wedge on p, q
            for( auto params : W ) {
                // find q
                auto q_z = DT.incident_vertices( handles.at(params.at(0)) );
                while( ++q_z != handles.at(params.at(1)) ); // point to q
                const auto q_i(q_z);

                wedge<DEGREE>( DT, handles, closest, params, addToE_star, q_i, printLog );
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

    if( printLog && n <= 200 ) {
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
    }

    //
    //
    // END PRINTER NONSENSE
    //
    //

} // function BCC2012



} // namespace gsnunf

#endif // GSNUNF_BCC2012_H

