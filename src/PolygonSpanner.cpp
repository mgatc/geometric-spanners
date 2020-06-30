#include "PolygonSpanner.h"

#include <list>
#include <unordered_set>
#include <math.h>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <cassert>



namespace gsnunf {

PolygonSpanner::PolygonSpanner( SpanningGraph &SG ) : DelaunayGraph( SG ) {

    // Process v_1
    auto v_1 = _DT->finite_vertices_begin(); // used for testing
    //advance(v_1, 6);

    Vertex_handle v_i = v_1; // choose v_1

    // Create and orient C so it points to s_1
    Vertex_circulator C = _DT->incident_vertices( v_i ); // neighbors of v_1 in DT
    Incident_vertices N_PS = _E.find( v_i )->second;     // neighbors of v_1 in PolygonSpanner edges
    Incident_vertices N_SG = SG._E.find( v_i )->second;  // neighbors of v_1 in SpanningGraph edges

    find_s_1_in_circulator( C, N_SG );

    Vertex_circulator done(C);                      // remember where we started
    Vertex_handle last = _DT->infinite_vertex();    // initialize last as an invalid, but allowed value

    do {
        if( N_SG.find( C->handle() ) != N_SG.end() ) {        // if we found an vertex, try to add edges for partition
            if( last != _DT->infinite_vertex() ) {        // not a valid partition if last vertex was infinite
                add_polygon_spanner_edges( last, v_i, C );
            }
            last = C;                                   // update last found vertex
        }
        _level.push( C->handle() );  // queue C
        _known.insert( C->handle() );// mark C as explored
    } while( --C != done && !_DT->is_infinite(C) ); // keep going until we reach done or infinite

    if( C == done ) {                               // if we reached done, we still need to add the last partition
        add_polygon_spanner_edges( last, v_i, C );
    }

    _known.insert( v_i );

    Vertex_handle s_1, s_m;
    Vertex_circulator s_j, s_k;

    // process v_i, i>1
    do { // loop through _level queue
        v_i = _level.front();
        _level.pop();

        C = _DT->incident_vertices( v_i );
        N_PS = _E.find( v_i )->second;     // neighbors of v_1 in PolygonSpanner edges
        N_SG = SG._E.find( v_i )->second;  // neighbors of v_1 in SpanningGraph edges

        assert( N_PS.size() <= 5 ); // Lemma 3.3

        find_s_1_in_circulator( C, N_SG );

        Vertex_circulator done(C);
        Vertex_circulator last;

        do { // loop through neighbors in DT
            if( N_SG.find( C->handle() ) != N_SG.end() ) {    // if we found an vertex in SG, try to add edges for partition
                if( last != nullptr ) {        // not a valid partition if last vertex was null
                    s_1 = last->handle();
                    s_j = last;
                    s_k = C;
                    s_m = C->handle();

                    while( N_PS.find( (--s_j)->handle() ) == N_PS.end() && s_j->handle() != s_1 ); //edge(++s_j, v_i) is not in E_P
                    while( N_PS.find( (++s_k)->handle() ) == N_PS.end() && s_k->handle() != s_m ); //edge(--s_k, v_i) is not in E_P

                    add_polygon_spanner_edges( s_1, v_i, s_j );
                    add_polygon_spanner_edges( s_k, v_i, s_m );

                    N_PS = _E.find( v_i )->second;     // update active edge list for this loop
                }
                last = C;                              // update last found vertex
            }
        } while( --C != done && !_DT->is_infinite(C) ); // keep going until we reach done or infinite

        if( C == done ) {                               // if we reached done, we still need to add the last partition
            s_1 = last->handle();
            s_j = last;
            s_k = C;
            s_m = C->handle();

            while( N_PS.find( (--s_j)->handle() ) == N_PS.end() && s_j->handle() != s_1 ); // while edge(++s_j, v_i) is not in E_P
            while( N_PS.find( (++s_k)->handle() ) == N_PS.end() && s_k->handle() != s_m );    // while edge(--s_k, v_i) is not in E_P

            add_polygon_spanner_edges( s_1, v_i, s_j );
            add_polygon_spanner_edges( s_k, v_i, s_m );
        }

        done = C;

        // BFS housekeeping
        do { // loop through neighbors in DT
            if( _known.find(C) == _known.end() && !_DT->is_infinite(C) ) { // If C is NOT known, queue it and add to known
                _level.push(C);
                _known.insert(C);
            }
        } while( --C != done ); // keep going until we reach done or infinite

    } while( !_level.empty() ); // level is not empty

} // PolygonSpanner::PolygonSpanner( SpanningGraph &SG )



void PolygonSpanner::add_cross_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
    Vertex_circulator N = _DT->incident_vertices(q);

    while( --N != p ); // loop until N points to p
    do {
        Vertex_handle v_1 = N->handle();
        Vertex_handle v_2 = (--N)->handle();
        if(  _known.find(v_1) == _known.end() && !_DT->is_infinite(v_1)
          && _known.find(v_2) == _known.end() && !_DT->is_infinite(v_2) ) { // N is not known or infinite
            add_edge( v_1, v_2 ); // add edge between N and CW N
        }
    } while( N != r );

}

void PolygonSpanner::add_forward_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
    double alpha = get_angle( p,q,r );
    short subangles = rint( ceil( 2.0 * alpha / PI ) );
    double beta = alpha / subangles;
    //cout<<"alpha:"<<alpha<< " subangles:"<<subangles<<" beta:"<<beta<<endl;

    Vertex_handle add[subangles];
    fill( add, add+subangles, _DT->infinite_vertex() ); // initialize add to infinite vertex

    Vertex_circulator N = _DT->incident_vertices(q);
    while( --N != p ); // loop until N points to p

    double theta;
    short i;

    while( --N != r ) {
        if( _known.find(N) == _known.end() && !_DT->is_infinite(N) ) { // N is not known or infinite
            theta = get_angle( p, q, N );
            i = int(theta/beta);
            //cout<< "  subangles: "<<subangles<< "    i: "<<i<< "    theta:"<<theta<<endl;
            assert( i <= subangles-1 );

            if( _DT->is_infinite( add[i] ) || Vector2D( N->point(), q->point() ).squared_length() < Vector2D( add[i]->point(), q->point() ).squared_length() )
                add[i] = N->handle();   // if the saved vertex is infinite or longer than the current one, update
        }
    } //while( --N != r );

    for( Vertex_handle v : add )
        if( !_DT->is_infinite(v) )
            add_edge(q, v);
}

void PolygonSpanner::add_polygon_spanner_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
    assert( !_DT->is_infinite(p) && !_DT->is_infinite(q) && !_DT->is_infinite(r) );

    add_forward_edges( p, q, r );
    add_cross_edges( p, q, r );
}

void PolygonSpanner::find_s_1_in_circulator( Vertex_circulator& C, const Incident_vertices& N ) {
    while( N.find( (--C)->handle() ) == N.end() ); // loop while N is not in E

    Vertex_circulator start(C);

    while( !( --C == start || _DT->is_infinite( C->handle() ) ) ); // loop until reaching s_1 again or an infinite vertex

    if( _DT->is_infinite( C->handle() ) ) // if we stopped on an infinite vertex, step CW
        start = --C;
}

} // namespace gsnunf
