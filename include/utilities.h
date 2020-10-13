#ifndef GSNUNF_UTILITIES_H
#define GSNUNF_UTILITIES_H

#include <chrono>
#include <cmath>
#include <iostream>
#include <string>

#include <CGAL/Kernel/global_functions.h>
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

template< typename Point >
inline double distance( Point p, Point q ) {
    return CGAL::sqrt( CGAL::squared_distance(p,q) );
}

template< class K >
double get_angle( const typename K::Point_2& p, const typename K::Point_2& q, const typename K::Point_2& r ) {
    CGAL::Vector_2<K> pq( p, q );
    CGAL::Vector_2<K> rq( r, q );

    double result = atan2( pq.y(), pq.x() ) - atan2( rq.y(), rq.x() );

    // atan() returns a value between -PI and PI. From zero ("up"), CCW rotation is negative and CW is positive.
    // Our zero is also "up," but we only want positive values between 0 and 2*PI:

    result = fmod( result+2*PI, 2*PI );
    //result = CGAL::min( result, 2*PI );
    //cout<<"angle("<<p<<","<<q<<","<<r<<")="<<result*180/PI<<" ";
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
