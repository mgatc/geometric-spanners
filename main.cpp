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


    const double width = 100;

    /*
        g1-g4 are the random point generators. Currently, they must be
        manually changed in the code to effect the point set properties.
        I tried to put them in a vector of the base class and loop through
        it to change the generator. However, their base class doesn't
        implement the ++ operator, required by copy_n. Therefore, we need
        to create a random point set factory for this purpose, which will
        be useful throughout the project.
    */

    //auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( width/2 );
   // auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>( width/2 );
    //auto g3 = CGAL::Random_points_on_circle_2<Point,Creator>( width/2 );
    //auto g4 = CGAL::Random_points_on_square_2<Point,Creator>( width/2 );

    size_t i = 10;
    GeometricSpannerPrinter printer( .25f );

//
//    //for( i=2; i<=10; ++i ) {
        // SET POINT SET
        list<Point> points = {
            {
                0,0
            },
            {
                0,2
            },
            {
                3,1
            },
            {
                -3,1
            },
            {
                2,-1
            },
            {
                -2,-1
            }
        };
        const int n = pow(2, i);
        //std::copy_n( g2, n, std::back_inserter(points) );

        auto start = chrono::steady_clock::now();

        Delaunay_triangulation_2 DT( points.begin(), points.end() );
        DelaunayGraph<Delaunay_triangulation_2> S(DT);

        PlanarSpanner(S);

        auto stop = chrono::steady_clock::now();

        cout
            <<i<<"--------------------------"
            <<" n:"<<n
            <<" w:"<<width
            <<" t:"<<StretchFactor(S)
            <<" b:"<<(PI+1)*(2*PI/(3*cos(PI/6)))
            <<" runtime:"<<chrono::duration_cast<chrono::microseconds>(stop - start).count()<<"us"
            <<"\n";
        //printer.draw( S, "temp");
    //}



    return 0;
}
