#include <list>
#include <memory>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // CGAL
#include <CGAL/Delaunay_triangulation_2.h>                      // Triangulations
#include <CGAL/point_generators_2.h>                            // Random point generation, testing

#include "DelaunayGraph.h"
#include "PlanarSpanner.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation_2;
typedef CGAL::Creator_uniform_2<double,Point> Creator;

using namespace gsnunf;
using namespace std;



int main() {

    const double epsilon = 2;

    // RANDOM POINT SET
//    list<Point> points;
//    const double width = 25;
//    const int n = 500;
//    CGAL::Random_points_in_square_2<Point,Creator> g( width/2 );
//    std::copy_n( g, n, std::back_inserter(points) );

//    for( Point p : points )
//        cout<<p<<endl;


    // POINT SET FROM PAPER, PAGE 253
    list<Point> points = {
        Point( -1, 0.1 ),
        Point( -0.9, 3 ),
        Point( -2, 6 ),
        Point( -7, 3.1 ),
        Point( -6, -0.1 ),
        Point( -9, -0.2 ),
        Point( -7.7, -1 ),
        Point( -6.1, -1.5 ),
        Point( -10, -4 ),
        Point( -4, -3 ),
        Point( -1.5, -6 ),
        Point( 1, -9 ),
        Point( 4, -4 ),
        Point( 4.1, 0 ),
        Point( 3.9, 5.9 ),
        Point( 5, 3 ),
        Point( 5, -2 ),
        Point( 9, 1 )
    };

    // TESTING POINT SET
//    list<Point> points = {
//        Point(0,0),
//        Point(0,3),
//        Point(5,0),
//        Point(5,3),
//        Point(7,5),
//        Point(7,-2),
//        Point(-4,4),
//        Point(-2,1),
//        Point(7,-9),
//        Point(-11,-3),
//        Point(10,0),
//        Point(1,-10),
//        Point(5,2),
//        Point(8,8)
//    };

    Delaunay_triangulation_2 DT( points.begin(), points.end() );

    DelaunayGraph<Delaunay_triangulation_2> S( DT );

    PlanarSpanner( S, epsilon );

    return 0;
}
