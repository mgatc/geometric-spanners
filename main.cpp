#include <chrono>
#include <list>
#include <utility>

#include <CGAL/point_generators_2.h>                            // Random point generation, testing

#include "Timer.h"
#include "FloydWarshall.h"
#include "PlanarSpanner.h"
#include "StretchFactor.h"

using namespace gsnunf;
typedef CGAL::Creator_uniform_2<double,Point> Creator;

int main() {
    using namespace std;

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

    size_t i = 50;

    for( i=1; i<=17; ++i ) {
        auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( width*sqrt(i)/2 );
        auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   width*sqrt(i)/2 );
        auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( width*sqrt(i)/2 );
        auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( width*sqrt(i)/2 );
        // SET POINT SET
        list<Point> points;// = {
//            {
//                0,0
//            },
//            {
//                0,2
//            },
//            {
//                3,1
//            },
//            {
//                -3,1
//            },
//            {
//                2,-1
//            },
//            {
//                -2,-1
//            }
//        };
        const int n = i*250000;
        std::copy_n( g1, n/3, back_inserter(points) );
        std::copy_n( g2, n/3, back_inserter(points) );
        std::copy_n( g3, n/6, back_inserter(points) );
        std::copy_n( g4, n/6, back_inserter(points) );
        //points.emplace_back( 0, 0 );

        cout<<points.size();
        cout<<",";
        list< pair< Point, Point > > result;

        PlanarSpanner( points.begin(), points.end(), back_inserter(result) );

        cout<<"\n";
    }

    return 0;
}
