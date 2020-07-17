#include <chrono>
#include <list>
#include <memory>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // CGAL
#include <CGAL/Delaunay_triangulation_2.h>                      // Triangulations
#include <CGAL/point_generators_2.h>                            // Random point generation, testing

#include "DelaunayGraph.h"
#include "FloydWarshall.h"
#include "PlanarSpanner.h"
#include "StretchFactor.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation_2;
typedef CGAL::Creator_uniform_2<double,Point> Creator;

using namespace gsnunf;
using namespace std;



int main() {

    const double epsilon = 2;
    int i = 2;

    //for( i=1; i<=5; i++ ) {
        // RANDOM POINT SET
        list<Point> points;
        const double width = 10*i;
        const int n = pow(5, i);
        CGAL::Random_points_in_square_2<Point,Creator> g( width/2 );
        std::copy_n( g, n, std::back_inserter(points) );


        auto start = chrono::steady_clock::now();

        Delaunay_triangulation_2 DT( points.begin(), points.end() );
        DelaunayGraph<Delaunay_triangulation_2> S(DT);
        PlanarSpanner( S, epsilon );

        auto stop = chrono::steady_clock::now();

//        cout
//            <<i<<"--------------------------"
//            <<" n:"<<n
//            <<" w:"<<width
//            <<" t:"<<StretchFactor(S)
//            <<" runtime:"<<chrono::duration_cast<chrono::milliseconds>(stop - start).count()
//            <<"\n"<<endl;
    //}



//    for( Point p : points )
//        cout<<p<<endl;

    // POINT SET FOR STRETCH FACTOR TEST
//    list<Point> points{
//            {0,0},
//            {0,1},
//            {1,1},
//            {1,0}
//        };

    // POINT SET FROM PAPER, PAGE 253
//    list<Point> points = {
//        { -1, 0.1 },
//        { -0.9, 3 },
//        { -2, 6 },
//        { -7, 3.1 },
//        { -6, -0.1 },
//        { -9, -0.2 },
//        { -7.7, -1 },
//        { -6.1, -1.5 },
//        { -10, -4 },
//        { -4, -3 },
//        { -1.5, -6 },
//        { 1, -9 },
//        { 4, -4 },
//        { 4.1, 0 },
//        { 3.9, 5.9 },
//        { 5, 3 },
//        { 5, -2 },
//        { 9, 1 }
//    };

    // TESTING POINT SET
//    list<Point> points = {
//        {0,0},
//        {0,3},
//        {5,0},
//        {5,3},
//        {7,5},
//        {7,-2},
//        {-4,4},
//        {-2,1},
//        {7,-9},
//        {-11,-3},
//        {10,0},
//        {1,-10},
//        {5,2},
//        {8,8}
//    };


    return 0;
}
