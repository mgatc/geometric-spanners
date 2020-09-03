#include <chrono>
#include <list>
#include <utility>

#include <CGAL/point_generators_2.h>                            // Random point generation, testing

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "Timer.h"
#include "FloydWarshall.h"
#include "GeometricSpannerPrinter.h"
#include "BGS2002.h"
#include "LW2004_2.h"
#include "LW2004_3.h"
#include "metrics.h"

using namespace gsnunf;
typedef CGAL::Creator_uniform_2<double,Point> Creator;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

void experiment();
void scratch();
void stretchFactorAndDegreeExperiment();
void generateRandomPoints( vector<Point> &points, const size_t n, const string outputFileName = "" );

template< class OutputIterator >
void readPointsFromFile( OutputIterator out, const string outputFileName ) {
    ifstream in (outputFileName);
    if (in.is_open()) {
        double x,y;
        while ( in >> x >> y ) {
            *out = Point(x,y);
            ++out;
        }
        in.close();
    }
}

template< class OutputIterator >
void generateRandomPoints( size_t n, double size, OutputIterator pointsOut ) {
    typedef CGAL::Creator_uniform_2<double,Point> Creator;
    //Random_points_in_disc_2<Point,Creator> g(10);
    //Random_points_in_square_2<Point,Creator> g(10);

    auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( size );
    auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( size );

    //random_convex_set_2(n,std::back_inserter(points), g);
    // Random_points_on_circle_2<Point,Creator> gen(50);
    // std::copy_n( gen, n, std::back_inserter(points));

    vector<Point> points;
    points.reserve(n);

    std::copy_n( g3, n/2, back_inserter(points) );
    std::copy_n( g4, n/2, back_inserter(points) );

    // copy points to output iterator
    for( Point p : points )
        *(pointsOut++) = p;

    // copy points to file
    ofstream out;
    out.open( to_string(n) + "_" + to_string(size) + "x" + to_string(size) + ".txt", ios::trunc );
    for( Point p : points )
        out << p << endl;
    out.close();
}

int main() {
    //experiment();
    scratch();

    return 0;
}

void scratch() {
    using namespace std;

    GeometricSpannerPrinter printer(0.1);

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

        n = 300;

//        std::copy_n( g1, n/3, back_inserter(points) );
//        std::copy_n( g2, n/3, back_inserter(points) );
//        std::copy_n( g3, n/6, back_inserter(points) );
        std::copy_n( g4, n/6, back_inserter(points) );

        points.emplace_back( 0,0 );

        //readPointsFromFile( back_inserter( points ), "in2.txt" );


        cout<< points.size();
        cout<< ",";
        list< pair< Point, Point > > result;
        pair<pair<Vertex_handle,Vertex_handle>,double> t;

        // Get t of Delaunay triangulation
            { // scope it so it doesn't stay in memory
                DelaunayGraph Del( points.begin(), points.end() );
                cout<<degree(Del._DT);
                cout<<",";
                cout << weight( Del._DT );
                cout <<",";
//               Del.add_all_edges();
//                t = StretchFactor(Del);
//                cout<< t.second;
//                cout<<",";
            }

        {
//                Timer tim;
            LW2004_3( points.begin(), points.end(), back_inserter(result) );
            //BGS2002( points.begin(), points.end(), back_inserter(result) );
        }
        cout << degree( result.begin(), result.end() );
        cout <<",";
        cout << weight( result.begin(), result.end() )/2;
        cout <<",";
//            t = StretchFactor( result.begin(), result.end() );
//            cout<< t.second;
//            cout<<",";
//        result.clear();

//            {
//                Timer tim;
//                BGS2002( points.begin(), points.end(), back_inserter(result) );
//            }
//            t = StretchFactor( result.begin(), result.end() );
//            cout<< t.second;
//            cout<<",";
//            result.clear();

        cout<<"\n";

        printer.drawEdges( result.begin(), result.end() );
        //printer.drawVertexPair( t.first, {{"color","red"}} );
        printer.print( "big_t" );
        //printer.print("bgs2002");


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

    for( size_t trial=1; trial<=10; ++trial ) {
        for( i=1; i<=12; ++i ) {
            double size = width*sqrt(i)/2;
            auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( size );
            auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   size );
            auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( size );
            auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( size );
            // SET POINT SET
            list<Point> points;
            const int n = i*10000;
//            std::copy_n( g1, n/3, back_inserter(points) );
//            std::copy_n( g2, n/3, back_inserter(points) );
//            std::copy_n( g3, n/6, back_inserter(points) );
//            std::copy_n( g4, n/6, back_inserter(points) );
            //points.emplace_back( 0, 0 );

            generateRandomPoints( n, size, back_inserter(points) );

            //readPointsFromFile( back_inserter( points ), "40000_100.000000x100.000000.txt" );

            cout<< points.size();
            cout<< ",";
            list< pair< Point, Point > > result;
            pair<pair<Vertex_handle,Vertex_handle>,double> t;

            // Delaunay triangulation
            { // scope it so it doesn't stay in memory
                CGAL::Delaunay_triangulation_2<K> DT( points.begin(), points.end() );
                cout << degree(DT);
                cout << ",";
                cout << weight(DT);
                cout << ",";
            }

            {
                Timer tim;
                LW2004_3( points.begin(), points.end(), back_inserter(result) );
            }
            cout << degree( result.begin(), result.end() );
            cout <<",";
            cout << weight( result.begin(), result.end() );
            cout <<",";
////            t = StretchFactor( result.begin(), result.end() );
////            cout<< t.second;
////            cout<<",";
            result.clear();

            {
                Timer tim;
                BGS2002( points.begin(), points.end(), back_inserter(result) );
            }
            cout << degree( result.begin(), result.end() );
            cout << ",";
            cout << weight( result.begin(), result.end() );
            cout << ",";
//            t = StretchFactor( result.begin(), result.end() );
//            cout<< t.second;
//            cout<<",";
            result.clear();

            cout<<"\n";
        }
    }
}

