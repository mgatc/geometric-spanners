#ifndef GSNUNF_UTILITIES_H
#define GSNUNF_UTILITIES_H

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Vector_2.h>



namespace gsnunf {

using namespace std;

const double PI = M_PI;
const double EPSILON = 0.000001;
const double INF = std::numeric_limits<double>::infinity();
const double MAX = std::numeric_limits<double>::max();

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

template< class Point, class OutputIterator >
void readPointsFromFile( OutputIterator out, const string outputFileName ) {
    //typedef typename OutputIterator::value_type Point;
    ifstream in(outputFileName);
    if (in.is_open()) {
        double x,y;
        Point p;
        while ( in >> p ) {
            *out = p;
            ++out;
        }
        in.close();
    }
}

template< class Point, class OutputIterator >
string generateRandomPoints( size_t n, double size, OutputIterator pointsOut ) {
    typedef CGAL::Creator_uniform_2<double,Point> Creator;

    auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( size );
    auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   size );
    auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( size );
    auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( size );


    auto g1s = CGAL::Random_points_in_square_2<Point,Creator>( size/4 );
    auto g2s = CGAL::Random_points_in_disc_2<Point,Creator>(   size/4 );
    auto g3s = CGAL::Random_points_on_square_2<Point,Creator>( size/4 );
    auto g4s = CGAL::Random_points_on_circle_2<Point,Creator>( size/4 );

    vector<Point> points;
    points.reserve(n);

    std::copy_n( g1, n*2/9, back_inserter(points) );
    std::copy_n( g2, n/9, back_inserter(points) );
    std::copy_n( g3, n*2/18, back_inserter(points) );
    std::copy_n( g4, n/18, back_inserter(points) );

    std::copy_n( g1s, n/9, back_inserter(points) );
    std::copy_n( g2s, n*2/9, back_inserter(points) );
    std::copy_n( g3s, n/18, back_inserter(points) );
    std::copy_n( g4s, n*2/18, back_inserter(points) );

    //points.emplace_back(0,0);

    // copy points to output iterator
    for( Point p : points )
        *(pointsOut++) = p;

    // copy points to file
    ofstream out;
    string fName;
    fName = "data-" + to_string(n) + "_" + to_string(size) + "x" + to_string(size) + ".txt";
    out.open( fName, ios::trunc );
    for( Point p : points )
        out << p << endl;

    out.close();

    return fName;
}

} // namespace gsnunf

#endif // GSNUNF_UTILITIES_H
