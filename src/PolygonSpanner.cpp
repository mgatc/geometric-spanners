#include "PolygonSpanner.h"

#include <list>
#include <unordered_set>
#include <math.h>
#include <cmath>

#include <CGAL/Vector_2.h>



namespace gsnunf {

PolygonSpanner::PolygonSpanner( SpanningGraph &SG ) : DelaunayGraph( SG ) {

    // Process v_1
    const Vertex_handle v_1 = (_DT->finite_vertices_begin()); // choose v_1
    _known.insert( v_1 );

    // Create and orient C so it points to s_1
    Vertex_circulator C = _DT->incident_vertices( v_1 ); // neighbors of v_1 in DT
    Incident_vertices N = _E.find( v_1 )->second;    // neighbors of v_1 in E

    while( N.find( (--C)->handle() ) == N.end() ); // loop while N is not in E

    Vertex_circulator start = C;

    while( !( --C == start || _DT->is_infinite( C->handle() ) ) ); // loop until reaching s_1 again or an infinite vertex

    if( _DT->is_infinite( C->handle() ) ) // if we stopped on an infinite vertex, step CW
        start = --C;

    Vertex_circulator done(start);
    Vertex_handle last = _DT->infinite_vertex();

    do {
        _level.push( C->handle() );  // queue C
        _known.insert( C->handle() );// mark C as explored
        if( N.find( C->handle() ) != N.end() ) {        // if we found an vertex, try to add edges for partition
            if( last != _DT->infinite_vertex() ) {        // not a valid partition if last vertex was infinite
                add_polygon_spanner_edges( last, v_1, C );
            }
            last = C;                                   // update last found vertex
        }
    } while( --C != done && !_DT->is_infinite(C) ); // keep going until we reach done or infinite

    if( C == done ) {                                 // if we reached done, we still need to add the last partition
        add_polygon_spanner_edges( last, v_1, C );
    }

    Vertex_handle v_i, s_1, s_m;
    Vertex_circulator s_j, s_k;

    // process v_i, i>1
    do {
        v_i = _level.front();
        _level.pop();

        C = _DT->incident_vertices( v_i );
        N = _E.find( v_i )->second;

        while( N.find( (--C)->handle() ) == N.end() ); // loop while C is not in N

        Vertex_circulator start = C;

        while( !( --C == start || _DT->is_infinite( C->handle() ) ) ); // loop until reaching s_1 again or an infinite vertex

        if( _DT->is_infinite( C->handle() ) ) // if we stopped on an infinite vertex, step CW
            start = --C;

        Vertex_circulator done(start);
        Vertex_handle last = _DT->infinite_vertex();

        // found s_1

        do {
            if( _known.find(C) == _known.end() ) {
                _level.push(C);
                _known.insert(C);
            }
            if( N.find( C->handle() ) != N.end() ) {        // if we found an vertex, try to add edges for partition
                if( last != _DT->infinite_vertex() ) {        // not a valid partition if last vertex was infinite
                    s_1 = last;
                    s_j = s_1;
                    s_k = C;
                    ++s_k;
                    s_m = C;

                    while( _E.find( (--s_j)->handle() ) == _E.end() ); //edge(++s_j, v_i) is not in E_P
                    while( _E.find( (++s_k)->handle() ) == _E.end() ); //edge(--s_k, v_i) is not in E_P

                    add_polygon_spanner_edges( s_1, v_i, s_j );
                    add_polygon_spanner_edges( s_k, v_i, s_m );
                }
                last = C;                                   // update last found vertex
            }
        } while( --C != done && !_DT->is_infinite(C) ); // keep going until we reach done or infinite

        if( C == done ) {                                 // if we reached done, we still need to add the last partition
            s_1 = done;
            s_j = done;
            s_k = C;
            ++s_k;
            s_m = C;

            while( _E.find( (--s_j)->handle() ) == _E.end() ); //edge(++s_j, v_i) is not in E_P
            while( _E.find( (++s_k)->handle() ) == _E.end() ); //edge(--s_k, v_i) is not in E_P

            add_polygon_spanner_edges( s_1, v_i, s_j );
            add_polygon_spanner_edges( s_k, v_i, s_m );
        }

    } while( !_level.empty() ); // level is not empty

} // PolygonSpanner::PolygonSpanner( SpanningGraph &SG )



void PolygonSpanner::add_cross_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
    Vertex_circulator N = _DT->incident_vertices(q);
    Vertex_handle v_prev = p;

    while( --N != p ); // loop until N points to p
    do {
        add_edge( N, --N ); // add edge between N and CW N
    } while( N != r );

}

void PolygonSpanner::add_forward_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
    double theta = get_angle( p,q,r );
    short num_subangles = ceil( theta / PI / 2 );
    double subangle = theta / num_subangles;

    Vertex_handle add[num_subangles];
    Vertex_circulator N = _DT->incident_vertices(q);

    while( --N != p );

    double alpha;
    short i;

    do {
        if( _known.find(N) == _known.end() && !_DT->is_infinite(N) ) { // N is not removed or infinite
            alpha = get_angle( N, q, r );
            i = int(alpha/theta);
            cout<< "sub: "<<num_subangles<< "  i: "<<i<<endl;
            Vector2D a( N->point(),q->point() );
            Vector2D b( add[i]->point(),q->point() );

            if( a.squared_length() < b.squared_length() )
                add[i] = N->handle();
        }
    } while( --N != q );

    for( Vertex_handle v : add )
        add_edge(q, v);
}

void PolygonSpanner::add_polygon_spanner_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
    cout<< "v_1:  " << q->point() <<endl;
    cout<< " s_1: " << p->point() <<endl;
    cout<< " s_m: " << r->point() <<endl;

    //add_forward_edges( p, q, r );
    add_cross_edges( p, q, r );
}

double PolygonSpanner::get_angle( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
    Vector2D pq(p -> point(), q -> point()); //vector pq <- q-p // geometric vector, not std::vector
    Vector2D qr(q -> point(), r -> point()); //    vector qr <- r-q
    return acos( ((pq.x() * qr.x()) + pq.y() * qr.y()) / (sqrt(pq.squared_length())*sqrt(qr.squared_length())));
}

} // namespace gsnunf

