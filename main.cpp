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

    /*
        Delaunay Triangulation Test

        The Delaunay triangulation is known to be a t-spanner
        where t <= 2pi/(3cos(pi/6)). Therefore, in output from
        this section, t<=b verifies correct operation.
    */


    cout
        <<"--------------------------\n"
        <<"--------------------------\n"
        <<"Delaunay Triangulation Test\n"
        <<"--------------------------\n"
        <<"The Delaunay triangulation is known to be a t-spanner\n"
        <<"where t <= 2pi/(3cos(pi/6)). Therefore, in output from\n"
        <<"this section, t<=b verifies correct operation.\n"
        <<"--------------------------\n";


    const double epsilon = 2;
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

    auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( width/2 );
    auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>( width/2 );
    auto g3 = CGAL::Random_points_on_circle_2<Point,Creator>( width/2 );
    auto g4 = CGAL::Random_points_on_square_2<Point,Creator>( width/2 );

    size_t i = 2;

    for( i=1; i<=10; ++i ) {
        // SET POINT SET
        list<Point> points;
        const int n = pow(2, i);
        std::copy_n( g4, n, std::back_inserter(points) );

        auto start = chrono::steady_clock::now();

        Delaunay_triangulation_2 DT( points.begin(), points.end() );
        DelaunayGraph<Delaunay_triangulation_2> S(DT);
        S.add_all_edges();

        //PlanarSpanner( S, epsilon );

        auto stop = chrono::steady_clock::now();

        cout
            <<i<<"--------------------------"
            <<" n:"<<n
            <<" w:"<<width
            <<" t:"<<StretchFactor(S)
            <<" b:"<<(2*PI/(3*cos(PI*6)))
            <<" runtime:"<<chrono::duration_cast<chrono::microseconds>(stop - start).count()<<"us"
            <<"\n";
    }

    /*
        Unit square test

        Four points placed in a 1x1 square with edges completing the square will
        yield t = sqrt2 ~ 1.41
    */

    list<Point> points {
        { 0,0 },
        { 0,1 },
        { 1,1 },
        { 1,0 }
    };

    Delaunay_triangulation_2 DT;
    vector<DelaunayGraph<Delaunay_triangulation_2>::Vertex_handle> V;

    for( auto p : points )
        V.emplace_back( DT.insert(p) );

    DelaunayGraph<Delaunay_triangulation_2> S(DT);

    for( size_t i=0; i<points.size(); ++i )
        S.add_edge( V.at(i), V.at( (i+1)%points.size() ) );


    cout
        <<"--------------------------\n"
        <<"--------------------------\n"
        <<"Unit Square Test\n"
        <<"--------------------------\n"
        <<"Four points placed in a 1x1 square with edges completing \n"
        <<" the square will yield t = sqrt2 ~ 1.41\n"
        <<"--------------------------\n"
        <<" t:"<<StretchFactor(S)
        <<" b:"<<sqrt(2)
        <<"\n";


    return 0;
}
