#include <chrono>
#include <list>
#include <utility>

#include <CGAL/point_generators_2.h>                            // Random point generation, testing

#include "Timer.h"
#include "FloydWarshall.h"
#include "BGS2002.h"
#include "LW2004_2.h"
#include "StretchFactor.h"

using namespace gsnunf;
typedef CGAL::Creator_uniform_2<double,Point> Creator;

void experiment();
void scratch();
void stretchFactorAndDegreeExperiment();

int main() {
    experiment();
    //scratch();

    return 0;
}

void scratch() {
    using namespace std;

    const double width = 100;
    size_t n = 30, i=n;

    //for( i=1; i<=17; ++i ) {
        auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( width*sqrt(i)/2 );
        auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   width*sqrt(i)/2 );
        auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( width*sqrt(i)/2 );
        auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( width*sqrt(i)/2 );
        // SET POINT SET
//        list<Point> points = {
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

            // POINT SET FROM PAPER, PAGE 253
    list<Point> points;// = {
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

        n = 3000000;

        std::copy_n( g1, n/3, back_inserter(points) );
        std::copy_n( g2, n/3, back_inserter(points) );
        std::copy_n( g3, n/6, back_inserter(points) );
        std::copy_n( g4, n/6, back_inserter(points) );

        cout<<points.size();
        cout<<",";

        double alpha = PI/2;
        list< pair< Point, Point > > result;

        LW2004_2( points.begin(), points.end(), back_inserter(result), alpha );
        //BGS2002( points.begin(), points.end(), back_inserter(result) );

        cout<<"\n";

}

void experiment() {
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

    GeometricSpannerPrinter printer;


    size_t i = 50;

    for( size_t trial=1; trial<=5; ++trial ) {
        for( i=1; i<=10; ++i ) {
            auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( width*sqrt(i)/2 );
            auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   width*sqrt(i)/2 );
            auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( width*sqrt(i)/2 );
            auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( width*sqrt(i)/2 );
            // SET POINT SET
            list<Point> points;
            const int n = i*25;
            std::copy_n( g1, n/3, back_inserter(points) );
            std::copy_n( g2, n/3, back_inserter(points) );
            std::copy_n( g3, n/6, back_inserter(points) );
            std::copy_n( g4, n/6, back_inserter(points) );
            //points.emplace_back( 0, 0 );

            cout<< points.size();
            cout<< ",";
            list< pair< Point, Point > > result;

            LW2004_2( points.begin(), points.end(), back_inserter(result) );

            pair<pair<Vertex_handle,Vertex_handle>,double> t = StretchFactor( result.begin(), result.end() );
            cout<< t.second;
            cout<<",";
            if(t.second>50){
                printer.drawEdges( result.begin(), result.end() );
                printer.drawVertexPair( t.first, {{"color","red"}} );
                printer.print( "big_ole_t" );
                return;
            }

            result.clear();

            BGS2002( points.begin(), points.end(), back_inserter(result) );
            cout<< StretchFactor( result.begin(), result.end() ).second;
            cout<<",";

            cout<<"\n";
        }
    }
}

void stretchFactorAndDegreeExperiment() {
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

    for( size_t trial=1; trial<=5; ++trial ) {
        for( i=1; i<=10; ++i ) {
            auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( width*sqrt(i)/2 );
            auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   width*sqrt(i)/2 );
            auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( width*sqrt(i)/2 );
            auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( width*sqrt(i)/2 );
            // SET POINT SET
            list<Point> points;
            const int n = i*250;
            std::copy_n( g1, n/3, back_inserter(points) );
            std::copy_n( g2, n/3, back_inserter(points) );
            std::copy_n( g3, n/6, back_inserter(points) );
            std::copy_n( g4, n/6, back_inserter(points) );
            //points.emplace_back( 0, 0 );

            cout<<points.size();
            cout<<",";
            list< pair< Point, Point > > result;

            BGS2002( points.begin(), points.end(), back_inserter(result) );

            cout<<"\n";
        }
    }
}
