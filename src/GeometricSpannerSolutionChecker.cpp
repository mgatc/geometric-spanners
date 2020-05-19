#include "GeometricSpannerSolutionChecker.h"
#include "CgalComponents.h"

namespace GeometricSpanner {

typedef CGAL::Delaunay_triangulation_2<K> Dt;
typedef Dt::Vertex_handle Vertex_handle;

GeometricSpannerSolutionChecker::GeometricSpannerSolutionChecker( double epsilon ) :
    epsilon( epsilon )
    { }

bool GeometricSpannerSolutionChecker::check( const vector<Point> &P, const list<Point> &C ) {
    bool pass = true;
    // Make new Delaunay Triangulation, insert all disk centers
    Dt T( C.begin(), C.end() );

    for( Point p : P ) {
        // get nearest neighbor to p from DT
        Point nearest = T.nearest_vertex( p )->point();
        if( CGAL::squared_distance( p, nearest ) > 1+epsilon ) {
            pass = false;
            std::cout << "\nSOLUTION ERROR: " << p << " not covered." << std::endl;
        }
    }
    return pass;
}

} // GeometricSpanner
