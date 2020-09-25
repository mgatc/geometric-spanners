#ifndef GSNUNF_UTILITIES_H
#define GSNUNF_UTILITIES_H

#include <cmath>

#include <CGAL/Vector_2.h>

namespace gsnunf {

using namespace std;

const double PI = M_PI;
const double EPSILON = 0.000001;

template< class T >
bool contains( const T& V, const typename T::key_type& v ) {
    return V.find(v) != V.end();
}

/* If V contains v, remove v.
 * If V does not contain v, add it.
 * Return new value.
 */
template< class T >
bool toggle( T& V, const typename T::key_type& v ) {
    bool found = contains( V, v );
    if( found ) V.erase( v );
    else V.insert( v );
    return !found;
}

template< typename T >
inline std::pair<T,T> makeNormalizedPair( const T& i, const T& j ) {
    return make_pair(
        CGAL::min(i,j),
        CGAL::max(i,j)
    );
}

template< class DT >
double get_angle( const DT& T, const typename DT::Vertex_handle &p, const typename DT::Vertex_handle &q, const typename DT::Vertex_handle &r ) {
    assert( !T.is_infinite(p) );
    assert( !T.is_infinite(q) );
    assert( !T.is_infinite(r) );

    CGAL::Vector_2 pq( p->point(), q->point() );
    CGAL::Vector_2 rq( r->point(), q->point() );

    double result = atan2( rq.y(), rq.x() ) - atan2( pq.y(), pq.x() );

    // atan() returns a value between -PI and PI. From zero ("up"), CCW rotation is negative and CW is positive.
    // Our zero is also "up," but we only want positive values between 0 and 2*PI:

    result *= -1; // First, invert the result. This will associate CW rotation with positive values.
    if( result < 0 ) // Then, if the result is less than 0 (or epsilon for floats) add 2*PI.
        result += 2*PI;
    result = CGAL::min( result, 2*PI );
    //cout<<"angle("<<p->info()<<","<<q->info()<<","<<r->info()<<")="<<result<<" ";
    return result;
}

template< typename TriWithInfo, typename OutputIterator >
void getVertexInfo( const TriWithInfo& Triangulation, OutputIterator res ) {
    for( auto v=Triangulation.finite_vertices_begin(); v!=Triangulation.finite_vertices_end(); ++v ) {
        *res = to_string( v->info() );
       ++res;
    }
}



} // namespace gsnunf

#endif // GSNUNF_UTILITIES_H
