#include "PolygonSpanner.h"

#include <list>



namespace gsnunf {

PolygonSpanner::PolygonSpanner( const SpanningGraph &SG, VisitsAllowedTable visitsAllowed ) : Graph( SG ) {
//    queue <- empty queue
//    E_P <- copy( SG._E )
//
//    // process v_1
//    choose v_1
//    v_1.is_removed = true
//    circulator N <- incident_vertices( v_1 )
//    orient N so it points to s_1
//    circulator done <- N
//
//    do
//    queue(N)
//    while ++N != done AND N is not infinite
//
//    s_1 <- done, s_m <- --N
//
//    add_polygon_spanner_edges( s_1, v_1, s_m )
//
//    // process v_i, i>1
//    do
//    v_i <- queue.pop()
//    v_i.is_removed = true
//    N <- incident_vertices( v_n )
//
//    orient N so it points to s_1
//    circulator done <- N
//
//    do
//      if N is not removed
//        queue(N)
//    while ++N != done AND N is not infinite
//
//    s_1 <- done, s_j <- done, s_k <- --N s_m <- N
//
//    while edge(++s_j, v_i) is not in E_P
//    while edge(--s_k, v_i) is not in E_P

//    add_polygon_spanner_edges( s_1, v_i, s_j )
//    add_polygon_spanner_edges( s_k, v_i, s_m )

//    while queue is not emp

} // PolygonSpanner::PolygonSpanner( SpanningGraph &SG )



void PolygonSpanner::add_cross_delaunay_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
//    Vertex_circulator C = incident_vertices(q)
//    Vertex_handle v_prev = p
//    while --C != p
//    while v_prev <- C AND --C != r
//    add edge( v_prev, C )
//    add_edge(v_prev, C)
}

void PolygonSpanner::add_incident_delaunay_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
//    theta <- GetAngle( p,q,r )
//    num_subangles <- ceiling( theta / 90 ) // may need radians instead of degrees (pi/2)
//    subangle <- theta / num_subangles
//    Vertex_handle add[num_subangles]
//    Vertex_circulator C = incident_vertices(q)
//
//    while --C != p
//
//    while --C != p
//    if C is not removed or infinite
//      alpha <- GetAngle( C, q, r )
//      i <- int(alpha/theta)
//      if vector(q,C) < vector(q,add[i])
//        add[i] = C
//    foreach v_n in add
//    add edge(q, v_n)
}

void PolygonSpanner::add_polygon_spanner_edges( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
    add_incident_delaunay_edges( p, q, r );
    add_cross_delaunay_edges( p, q, r );
}

double PolygonSpanner::get_angle( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
//    vector pq <- q-p // geometric vector, not std::vector
//    vector qr <- r-q
//    return arccos( dotproduct(pq, qr)/(length(pq)*length(qr)))
}

} // namespace gsnunf

