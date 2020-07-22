#ifndef GSNUNF_POLYGONSPANNER_H
#define GSNUNF_POLYGONSPANNER_H

#include <algorithm>
#include <queue>

#include "DelaunayGraph.h"
#include "GeometricSpannerPrinter.h"
#include "SpanningGraph.h"
#include "StretchFactor.h"



namespace gsnunf {

namespace polygon_spanner {

    enum VertexStatus { unknown, known, complete };

    template< typename T >
    void add_cross_edges( T& G, const typename T::template VertexMap<VertexStatus>& status,
                         const typename T::Vertex_handle &p, const typename T::Vertex_handle &q, const typename T::Vertex_handle &r ){
        typename T::Vertex_circulator N = G._DT.incident_vertices(q);

        while( --N != p ); // loop until N points to p
        do {
            typename T::Vertex_handle v_1 = N->handle();
            typename T::Vertex_handle v_2 = (--N)->handle();
            G.addToEventQueue( v_1, 1 ); // focus1 on v_1
            G.addToEventQueue( v_2, 2 );// focus2 on v_2

            // get status of v_1 and v_2, if they are complete, do not add any edges
            if( status.at(v_1) != complete && !G._DT.is_infinite(v_1)
              && status.at(v_2) != complete && !G._DT.is_infinite(v_2) ) {
                G.add_edge( v_1, v_2 ); // add edge between N and CW N
                G.addToEventQueue( {v_1, v_2}, true );
            }
        } while( N != r );
    }

    template< typename T >
    void add_forward_edges( T& G, const typename T::template VertexMap<VertexStatus>& status,
                           const typename T::Vertex_handle &p, const typename T::Vertex_handle &q, const typename T::Vertex_handle &r ) {

        using Vertex_handle = typename T::Vertex_handle;
        using Vertex_circulator = typename T::Vertex_circulator;
        using Vector_2 = typename T::Vector_2;

        double alpha = G.get_angle( p,q,r );
        short subangles = rint( ceil( 2.0 * alpha / PI ) );
        double beta = alpha / subangles;

        std::vector<Vertex_handle> add( subangles, G._DT.infinite_vertex() ); // initialize add to infinite vertex

        Vertex_circulator N = G._DT.incident_vertices(q);
        while( --N != p ); // loop until N points to p

        double theta;
        short i;
        double deg = 180/PI;

        //cout<<p->point()<<" "<<q->point()<<" "<<r->point()<<" alpha:"<<alpha*deg<<" subangles:"<<subangles<<" beta:"<<beta*deg<<"\n";

        while( --N != r ) {
            G.addToEventQueue( N, 1 );// focus1 on N
            if( !G._DT.is_infinite(N) && status.at(N) != complete ) { // N is not infinite or complete
                theta = G.get_angle( p, q, N );
                i = int(theta/beta);
                assert( i <= subangles-1 );
                //cout<<"  theta:"<<theta*deg<<" i:"<<i<<"\n";

                if( G._DT.is_infinite( add[i] )
                  || Vector_2( N->point(), q->point() ).squared_length() < Vector_2( add[i]->point(), q->point() ).squared_length() )
                    add[i] = N->handle();   // if the saved vertex is infinite or longer than the current one, update
            }
        } //while( --N != r );

        for( Vertex_handle v : add )
            if( !G._DT.is_infinite(v) ) {
                if( G._E.find( v )->second.size() < 5 )   // neighbors of v_1 in PolygonSpanner edges
                    G.add_edge(q, v);
                else {

//                    std::cout<<v->point()<<" N_PS:"<<G._E.find( v )->second.size()<<std::endl;
//
//                    GeometricSpannerPrinter printer( .25f );
//                    printer.draw( G, "PolygonSpanner" );
                }

                G.addToEventQueue( {q, v}, true ); // add edge
                //cout<<"  adding edge:"<<q->point()<<" "<<v->point()<<"\n";
            }
    }

    template< typename T >
    void add_polygon_spanner_edges( T& G, const typename T::template VertexMap<VertexStatus>& status,
                                   const typename T::Vertex_handle &p, const typename T::Vertex_handle &q, const typename T::Vertex_handle &r ) {
        assert( !G._DT.is_infinite(p) && !G._DT.is_infinite(q) && !G._DT.is_infinite(r) );

        add_forward_edges( G, status, p, q, r );
        //add_cross_edges( G, status, p, q, r );
    }

    template< typename T >
    void find_s_1( const T& G, typename T::Vertex_circulator& C, const typename T::VertexSet& N ) {
        typename T::Vertex_circulator done(C);
        // NOTE: use --C/C-- for CW rotation and ++C/C++ for CCW rotation
        //cout<< "C is pointed at "<<C->point()<<"to start\n";
        // Loop until C is in N.
        while( !contains( N, --C ) && C != done );
        //cout<< "C is pointed at "<<C->point()<<" after the first loop\n";
        done = C;
        /* Loop until the last C was an infinite vertex or this
         * vertex is the starting vertex. If C is infinite, it
         * will be decremented once more by the postfix operator,
         * landing us on the first vertex CW of the infinite, a
         * valid s_1. If the starting vertex is reached instead
         * of the infinite, then it is a valid s_1.
         */
        while( !G._DT.is_infinite(C--) && C != done );
        //cout<< "C is pointed at "<<C->point()<<" after the second loop\n";
    }

} // namespace polygon_spanner



template< typename T >
void PolygonSpanner( DelaunayGraph<T>& P ) {

    using Vertex_handle = typename T::Vertex_handle;
    using Vertex_circulator = typename T::Vertex_circulator;
    using VertexHash = typename DelaunayGraph<T>::VertexHash;
    using VertexSet = typename DelaunayGraph<T>::VertexSet;

    using namespace polygon_spanner;

    // Create a vertex status map and fill it
    typename DelaunayGraph<T>::template VertexMap<VertexStatus> status;
    for( auto it = P._DT.finite_vertices_begin(); it != P._DT.finite_vertices_end(); ++it )
        status.emplace( it, VertexStatus::unknown );

    VertexHash on_outer_face; // property maps
    std::queue<Vertex_handle> level; // BFS queue
    DelaunayGraph PS( P );          // resultant graph object

    Vertex_circulator v_convex_hull = P._DT.incident_vertices( P._DT.infinite_vertex() ), // create a circulator of the convex hull
                      done( v_convex_hull );

    do { // cycle through convex hull to set vertex info
        on_outer_face.insert( v_convex_hull );  // incident_chords (called in next loop) relies on accurate on_outer_face values
        P.addToEventQueue( v_convex_hull, 0 );// focus0 on v_convex_hull
    } while( ++v_convex_hull != done );

    /* Process v_1

       Processing v_1 is slightly different than processing all other vertices.
       Theoretically, v_1 can be chosen arbitrarily, but because we have only
       symbolically duplicated vertices and edges for TransformPolygon, we must
       ensure v_1 is not a split vertex. A split vertex is any vertex on the
       convex hull with degree > 2, or any interior vertex with degree > 1.
       These two conditions can be represented in one statement by casting a
       boolean representation of "on convex hull" to int.
     */

    auto v_1 = P._DT.finite_vertices_begin();
    VertexSet& N_P = P._E.find( v_1 )->second;  // neighbors of v_1 in SpanningGraph edges
    do {
        N_P = P._E.find( v_1 )->second;
    } while( N_P.size() > 1 + int( contains( on_outer_face, v_1++ ) ) );

    Vertex_handle v_i = v_1; // choose v_1
    P.addToEventQueue( v_i, 0 );

    // Create and orient C so it points to s_1
    Vertex_circulator C = PS._DT.incident_vertices( v_i ); // neighbors of v_1 in DT
    VertexSet& N_PS = PS._E.find( v_i )->second;  // neighbors of v_1 in PolygonSpanner edges
    find_s_1( PS, C, N_P ); // verified

    done = C;                      // remember where we started
    Vertex_handle last = PS._DT.infinite_vertex(); // initialize last as infinite

    do { // investigate neighbors

        P.addToEventQueue( C, 1 ); // focus1 on C
        if( contains( N_P, C ) ) {        // if we found a vertex, try to add edges for partition
            //cout<<"found neighbor of "<<v_i->point()<<":  "<<C->point()<<"\n";
            if( last != PS._DT.infinite_vertex() ) {        // not a valid partition if last vertex was infinite
                add_polygon_spanner_edges( PS, status, last, v_i, C );
            }
            last = C;                                   // update last found vertex
        }
        // BFS Housekeeping
        level.push(C);  // queue C
        status.at(C) = known; // mark C as explored
    } while( --C != done && !PS._DT.is_infinite(C) ); // keep going until we reach done or infinite

    if( C == done && last != PS._DT.infinite_vertex() ) {                               // if we reached done, we still need to add the last partition
        add_polygon_spanner_edges( PS, status, last, v_i, C );
    }
    status.at(v_i) = complete;

    //cout<<"v_1:"<<v_i->point()<<endl;

    // Process v_i, i>1 //

    Vertex_handle s_1, s_m;
    Vertex_circulator s_j, s_k;

    do { // loop through level queue

        v_i = level.front();
        P.addToEventQueue( v_i, 0 ); // focus0 on v_i

        C = PS._DT.incident_vertices( v_i );
        N_PS = PS._E.find( v_i )->second;    // neighbors of v_1 in PolygonSpanner edges
        N_P = P._E.find( v_i )->second;  // neighbors of v_1 in Polygon edges

        if( N_PS.size() > 5 ) {
            std::cout<<v_i->point()<<" N_PS:"<<N_PS.size()<<" N_P:"<<N_P.size()<<std::endl;

//            GeometricSpannerPrinter printer( .25f );
//            printer.draw( PS, "PolygonSpanner" );
        }
        assert( N_PS.size() <= 5 ); // Lemma 3.3

        find_s_1( PS, C, N_P );

        Vertex_circulator done(C);
        Vertex_circulator s_1;

        do { // loop through neighbors in DT

            P.addToEventQueue( C, 1 ); // focus1 on C
            if( contains( N_P, C ) ) {    // if we found a vertex in P, try to add edges for partition
                if( CGAL::circulator_size( s_1 ) > 0 ) {        // not a valid partition if last vertex's neighbors
                    s_j = s_1;
                    s_k = C;
                    s_m = C->handle();

                    // loop until s_j is in N_PS or we went all the way to s_m
                    while( !( contains( N_PS, --s_j ) || s_j == s_m ) );
                    while( !contains( N_PS, ++s_k ) && s_k != s_j ); //edge(--s_k, v_i) is not in E_P

                    cout<<"s1, vi, sj:"<<PS.get_angle(s_1,v_i,s_j)<<"\n";
                    cout<<"sk, vi, sm:"<<PS.get_angle(s_k,v_i,s_m)<<"\n";

                    add_polygon_spanner_edges( PS, status, s_1, v_i, s_j );
                    add_polygon_spanner_edges( PS, status, s_k, v_i, s_m );

                    N_PS = PS._E.find( v_i )->second;     // update active edge list for this loop
                }
                last = C;                              // update last found vertex
            }
        } while( --C != done && !PS._DT.is_infinite(C) ); // keep going until we reach done or infinite

        if( C == done ) {
            P.addToEventQueue( C, 1 ); // focus1 on C
            s_1 = last->handle();
            s_j = last;
            s_k = C;
            s_m = C->handle();

            while( !contains( N_PS, --s_j ) && s_j->handle() != s_1 ); // while edge(++s_j, v_i) is not in E_P
            while( !contains( N_PS, ++s_k ) && s_k->handle() != s_m ); // while edge(--s_k, v_i) is not in E_P

            add_polygon_spanner_edges( PS, status, s_1, v_i, s_j );
            add_polygon_spanner_edges( PS, status, s_k, v_i, s_m );
        }

        done = C;

        // BFS housekeeping
        level.pop();
        do { // loop through neighbors in DT
            if( !PS._DT.is_infinite(C) && status.at(C) == unknown ) { // If C is NOT known, queue it and add to known
                level.push(C);
                status.at(C) = known;
                P.addToEventQueue( C, 1 ); // focus1 on C
            }

        } while( --C != done ); // keep going until we reach done or infinite
    } while( !level.empty() ); // level is not empty
    P.addToEventQueue( P._DT.infinite_vertex(), 0 ); // focus0 on infinite

    std::swap( P._E, PS._E );

    // Lemma 3.4
    // ((PI+1)*(2*PI/(3*cos(PI/6)))) = 10.01602416
    //assert( StretchFactor(P) <= (PI+1) ); // fails
    //assert( StretchFactor(P) <= ((PI+1)*(2*PI/(3*cos(PI/6)))) ); // fails

    // Test degree assumption given after lemma 3.4
    for( auto it=P._E.begin(); it!=P._E.end(); ++it ) {
        assert( it->second.size() <= 12 );
    }

} // PolygonSpanner( SpanningGraph &P )

} // namespace gsnunf

#endif // GSNUNF_POLYGONSPANNER_H
