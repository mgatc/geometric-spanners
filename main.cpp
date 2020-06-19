#include "CGALComponents.h"
#include "Graph.h"
#include "PlanarSpanner.h"


using namespace gsnunf;


int main() {

    const double epsilon = 2;

    std::list<Point> points;
    Graph graph;

//    const double width = 50;
//    const int n = 50;
//    CGAL::Random_points_in_square_2<Point,Creator> g( width/2 );
//    std::copy_n( g, n, std::back_inserter(points) );

    points.push_back( Point(0,0) );
    points.push_back( Point(0,3) );
    points.push_back( Point(5,0) );
    points.push_back( Point(5,3) );
    points.push_back( Point(7,5) );
    points.push_back( Point(7,-2) );
    points.push_back( Point(-4,4) );
    points.push_back( Point(-2,1) );
    points.push_back( Point(7,-9) );
    points.push_back( Point(-11,-3) );
    points.push_back( Point(10,0) );
    points.push_back( Point(1,-10) );
    points.push_back( Point(5,2) );
    points.push_back( Point(8,8) );

    PlanarSpanner planar_spanner_builder( points, epsilon, graph );

    //graph.print();

    return 0;
}
