#ifndef GSNUNF_POLYGONSPANNER_H
#define GSNUNF_POLYGONSPANNER_H

#include <algorithm>
#include <queue>

#include "DelaunayGraph.h"
#include "GeometricSpannerPrinter.h"
#include "SpanningGraph.h"
//#include "StretchFactor.h"



namespace gsnunf {

namespace polygon_spanner {

    template< typename T >
    void add_cross_edges( T& G, const typename T::VertexHash& known,
                         const typename T::Vertex_handle &p, const typename T::Vertex_handle &q, const typename T::Vertex_handle &r ){
        typename T::Vertex_circulator N = G._DT.incident_vertices(q);

        while( --N != p ); // loop until N points to p
        do {
            typename T::Vertex_handle v_1 = N->handle();
            typename T::Vertex_handle v_2 = (--N)->handle();
            G.addToEventQueue( v_1, 1 ); // focus1 on v_1
            G.addToEventQueue( v_2, 2 );// focus2 on v_2
            if(  !contains( known, v_1 ) && !G._DT.is_infinite(v_1)
              && !contains( known, v_2 ) && !G._DT.is_infinite(v_2) ) { // N is not known or infinite
                G.add_edge( v_1, v_2 ); // add edge between N and CW N
                G.addToEventQueue( {v_1, v_2}, true );
            }
        } while( N != r );
    }

    template< typename T >
    void add_forward_edges( T& G, const typename T::VertexHash& known,
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

        while( --N != r ) {
            G.addToEventQueue( N, 1 );// focus1 on N
            if( !contains( known, N ) && !G._DT.is_infinite(N) ) { // N is not known or infinite
                theta = G.get_angle( p, q, N );
                i = int(theta/beta);
                assert( i <= subangles-1 );

                if( G._DT.is_infinite( add[i] ) || Vector_2( N->point(), q->point() ).squared_length() < Vector_2( add[i]->point(), q->point() ).squared_length() )
                    add[i] = N->handle();   // if the saved vertex is infinite or longer than the current one, update
            }
        } //while( --N != r );

        for( Vertex_handle v : add )
            if( !G._DT.is_infinite(v) ) {
                G.add_edge(q, v);
                G.addToEventQueue( {q, v}, true ); // add edge
            }
    }

    template< typename T >
    void add_polygon_spanner_edges( T& G, const typename T::VertexHash& known,
                                   const typename T::Vertex_handle &p, const typename T::Vertex_handle &q, const typename T::Vertex_handle &r ) {
        assert( !G._DT.is_infinite(p) && !G._DT.is_infinite(q) && !G._DT.is_infinite(r) );

        //add_forward_edges( G, known, p, q, r );
        //add_cross_edges( G, known, p, q, r );
    }

    template< typename T >
    void find_s_1( const T& G, typename T::Vertex_circulator& C, const typename T::VertexSet& N ) {
        typename T::Vertex_circulator done(C);

        while( !contains( N, --C ) && C != done ); // loop while C is not in N

        done = C;

        while( !( --C == done || G._DT.is_infinite( C->handle() ) ) ); // loop until reaching s_1 again or an infinite vertex

        if( G._DT.is_infinite( C->handle() ) ) // if we stopped on an infinite vertex, step CW
            done = --C;
    }

} // namespace polygon_spanner



template< typename T >
void PolygonSpanner( DelaunayGraph<T>& SG ) {

    using Vertex_handle = typename T::Vertex_handle;
    using Vertex_circulator = typename T::Vertex_circulator;
    using VertexHash = typename DelaunayGraph<T>::VertexHash;
    using VertexSet = typename DelaunayGraph<T>::VertexSet;

    using namespace polygon_spanner;

    VertexHash known, on_outer_face; // property maps
    std::queue<Vertex_handle> level; // BFS queue
    DelaunayGraph<T> PS( SG );          // resultant graph object

    Vertex_circulator v_convex_hull = SG._DT.incident_vertices( SG._DT.infinite_vertex() ), // create a circulator of the convex hull
                      done( v_convex_hull );

    do { // cycle through convex hull to set vertex info
        on_outer_face.insert( v_convex_hull );  // incident_chords (called in next loop) relies on accurate on_outer_face values
        SG.addToEventQueue( v_convex_hull, 0 );// focus0 on v_convex_hull
    } while( ++v_convex_hull != done );

    // Process v_1 //

    auto v_1 = SG._DT.finite_vertices_begin(); // v_1 can be chosen arbitrarily
    //advance(v_1, 6); // used for testing

    Vertex_handle v_i = v_1; // choose v_1
    SG.addToEventQueue( v_convex_hull, 0 );// focus0 on v_convex_hull

    // Create and orient C so it points to s_1
    Vertex_circulator C = PS._DT.incident_vertices( v_i ); // neighbors of v_1 in DT
    VertexSet& N_PS = PS._E.find( v_i )->second;  // neighbors of v_1 in PolygonSpanner edges
    VertexSet& N_SG = SG._E.find( v_i )->second;  // neighbors of v_1 in SpanningGraph edges
    find_s_1( PS, C, N_SG );

    done = C;                      // remember where we started
    Vertex_handle last = PS._DT.infinite_vertex(); // initialize last as infinite

    do { // investigate neighbors
        SG.addToEventQueue( C, 1 ); // focus1 on C
        if( contains( N_SG, C ) ) {        // if we found an vertex, try to add edges for partition
            if( last != PS._DT.infinite_vertex() ) {        // not a valid partition if last vertex was infinite
                add_polygon_spanner_edges( PS, known, last, v_i, C );
            }
            last = C;                                   // update last found vertex
        }
        level.push( C->handle() );  // queue C
        known.insert( C->handle() );// mark C as explored
    } while( --C != done && !PS._DT.is_infinite(C) ); // keep going until we reach done or infinite

    if( C == done ) {                               // if we reached done, we still need to add the last partition
        add_polygon_spanner_edges( PS, known, last, v_i, C );
    }

    known.insert( v_i );

    // Process v_i, i>1 //

    Vertex_handle s_1, s_m;
    Vertex_circulator s_j, s_k;

    do { // loop through _level queue
        v_i = level.front();
        level.pop();
        SG.addToEventQueue( v_i, 0 ); // focus0 on v_i

        C = PS._DT.incident_vertices( v_i );
        N_PS = PS._E.find( v_i )->second;    // neighbors of v_1 in PolygonSpanner edges
        N_SG = SG._E.find( v_i )->second;  // neighbors of v_1 in SpanningGraph edges

        // Lemma 3.3
        //assert( N_PS.size() <= 5 );
        std::cout<<v_i->point()<<" N_PS:"<<N_PS.size()<<" N_SG:"<<N_SG.size()<<std::endl;

        find_s_1( PS, C, N_SG );

        Vertex_circulator done(C);
        Vertex_circulator last;

        do { // loop through neighbors in DT
            SG.addToEventQueue( C, 1 ); // focus1 on C
            if( contains( N_SG, C ) ) {    // if we found a vertex in SG, try to add edges for partition
                if( last != nullptr ) {        // not a valid partition if last vertex was null
                    s_1 = last->handle();
                    s_j = last;
                    s_k = C;
                    s_m = C->handle();

                    while( !contains( N_PS, --s_j ) && s_j->handle() != s_1 ); //edge(++s_j, v_i) is not in E_P
                    while( !contains( N_PS, ++s_k ) && s_k->handle() != s_m ); //edge(--s_k, v_i) is not in E_P

                    add_polygon_spanner_edges( PS, known, s_1, v_i, s_j );
                    add_polygon_spanner_edges( PS, known, s_k, v_i, s_m );

                    N_PS = PS._E.find( v_i )->second;     // update active edge list for this loop
                }
                last = C;                              // update last found vertex
            }
        } while( --C != done && !PS._DT.is_infinite(C) ); // keep going until we reach done or infinite

        if( C == done ) {
            SG.addToEventQueue( C, 1 ); // focus1 on C
            s_1 = last->handle();
            s_j = last;
            s_k = C;
            s_m = C->handle();

            while( !contains( N_PS, --s_j ) && s_j->handle() != s_1 ); // while edge(++s_j, v_i) is not in E_P
            while( !contains( N_PS, ++s_k ) && s_k->handle() != s_m ); // while edge(--s_k, v_i) is not in E_P

            add_polygon_spanner_edges( PS, known, s_1, v_i, s_j );
            add_polygon_spanner_edges( PS, known, s_k, v_i, s_m );
        }

        done = C;

        // BFS housekeeping
        do { // loop through neighbors in DT
            if( !contains( known, C ) && !PS._DT.is_infinite(C) ) { // If C is NOT known, queue it and add to known
                level.push(C);
                known.insert(C);
                SG.addToEventQueue( C, 1 ); // focus1 on C
            }
        } while( --C != done ); // keep going until we reach done or infinite

    } while( !level.empty() ); // level is not empty

    SG.addToEventQueue( SG._DT.infinite_vertex(), 0 ); // focus0 on infinite

    std::swap( SG._E, PS._E );

    // Lemma 3.4
    // ((PI+1)*(2*PI/(3*cos(PI/6)))) = 10.01602416
    //assert( StretchFactor(SG) <= (PI+1) ); // fails
    //assert( StretchFactor(SG) <= ((PI+1)*(2*PI/(3*cos(PI/6)))) ); // fails

    // Test degree assumption given after lemma 3.4
    for( auto it=SG._E.begin(); it!=SG._E.end(); ++it ) {
        assert( it->second.size() <= 36 );
    }

} // PolygonSpanner( SpanningGraph &SG )

} // namespace gsnunf

#endif // GSNUNF_POLYGONSPANNER_H
