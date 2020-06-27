#include "PolygonSpanner.h"

#include <list>



namespace gsnunf {

PolygonSpanner::PolygonSpanner( const SpanningGraph &SG ) : Graph( SG ) {
//    level <- empty queue
//    known <- empty set (use unordered_set)
//    E_P <- copy( SG._E )

    // Process v_1
    const Vertex_handle v_1 = (SG._DT.finite_vertices_begin()); // choose v_1
//    known.insert( v_1 );
    cout << "v_1: " << v_1->point() << endl;

    // Create and orient C so it points to s_1
    Vertex_circulator C = SG._DT.incident_vertices( v_1 ); // neighbors of v_1 in DT
    Incident_vertices N = SG._E.find( v_1 )->second;    // neighbors of v_1 in E

    while( N.find( (--C)->handle() ) == N.end() ) cout << " rej1: " << C->point()<<endl; // loop while N is not in E

    Vertex_circulator start = C;
    cout<< " start: "<< start->point()<<endl; // prospective s_1

    while( !( --C == start || SG._DT.is_infinite( C->handle() ) ) ) cout << " rej2: "<<C->point()<<endl; // loop until reaching s_1 again or an infinite vertex

    if( SG._DT.is_infinite( C->handle() ) ) // if we stopped on an infinite vertex, step CW
        start = --C;

    Vertex_circulator done(start);
    Vertex_handle last = SG._DT.infinite_vertex();

    do {
      //level.insert( C->handle() );
      //known.insert( C->handle() );
        if( N.find( C->handle() ) != N.end() ) {        // if we found an vertex, try to add edges for partition
            if( last != SG._DT.infinite_vertex() ) {        // not a valid partition if last vertex was infinite
                add_polygon_spanner_edges( last, v_1, C );
            }
            last = C;                                   // update last found vertex
        }
    } while( --C != done && !SG._DT.is_infinite(C) ); // keep going until we reach done or infinite

    if( C == done ) {                                 // if we reached done, we still need to add the last partition
        add_polygon_spanner_edges( last, v_1, C );
    }

//    // process v_i, i>1
//    do
    //    v_i <- level.pop()
    //    N <- incident_vertices( v_n )
    //
    //    orient N so it points to s_1
    //    circulator done <- N
    //
    //    do
    //      if N is not known
    //        level.queue(N)
    //        known.add(N)
    //    while ++N != done AND N is not infinite
    //
    //    s_1 <- done, s_j <- done, s_k <- --N s_m <- N
    //
    //    while edge(++s_j, v_i) is not in E_P
    //    while edge(--s_k, v_i) is not in E_P

    //    add_polygon_spanner_edges( s_1, v_i, s_j )
    //    add_polygon_spanner_edges( s_k, v_i, s_m )

//    while level is not empty

} // PolygonSpanner::PolygonSpanner( SpanningGraph &SG )



void PolygonSpanner::add_cross_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
//    Vertex_circulator N = incident_vertices(q)
//    Vertex_handle v_prev = p
//    while --N != p
//    while v_prev <- N AND --N != r
//    add edge( v_prev, N )
//    add_edge(v_prev, N)
}

void PolygonSpanner::add_forward_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
//    theta <- GetAngle( p,q,r )
//    num_subangles <- ceiling( theta / 90 ) // may need radians instead of degrees (pi/2)
//    subangle <- theta / num_subangles
//    Vertex_handle add[num_subangles]
//    Vertex_circulator N = incident_vertices(q)
//
//    while --N != p
//
//    while --N != p
//    if N is not removed or infinite
//      alpha <- GetAngle( N, q, r )
//      i <- int(alpha/theta)
//      if vector(q,N) < vector(q,add[i])
//        add[i] = N
//    foreach v_n in add
//    add edge(q, v_n)
}

void PolygonSpanner::add_polygon_spanner_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {

    cout<< " s_1: " << p->point() <<endl;
    cout<< " s_m: " << r->point() <<endl;

    add_forward_edges( p, q, r );
    add_cross_edges( p, q, r );
}

double PolygonSpanner::get_angle( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
//    vector pq <- q-p // geometric vector, not std::vector
//    vector qr <- r-q
//    return arccos( dotproduct(pq, qr)/(length(pq)*length(qr)))
}

} // namespace gsnunf

