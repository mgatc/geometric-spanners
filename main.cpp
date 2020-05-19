#include "Greedy.h"

#include "GeometricSpannerExperiment.h"
#include "GeometricSpannerPrinter.h"
#include "CgalComponents.h"
#include <chrono>

#define PI 3.14159265

const int million = 1000000;
const int billion = 1000*million;
const int thousand = 1000;

void executeExp(const vector<Point> &points);
void pisExp();
void pidExp();
void convexSetExp();
void executeFromFile(string fName);

void simpleExperiment();
void runExperiment();
void scratch();

int main() {

    scratch();
    //simpleExperiment();
    //runExperiment();

    return EXIT_SUCCESS;
}


// The scratch function is used to keep temporary code out of main
void scratch() {

    double a = 10;   // size of point bound
    int n = 200;    // num points

    std::vector<Point> P; // problem
    std::list<Point> C;   // solution

    //Load from file
//    ifstream input("./experiments/problems/synthetic/square_200_6.txt");
//
//    long double x,y;
//    while(input >> x >> y )
//      P.emplace_back(Point(x,y));

    GeometricSpannerPrinter pointPrinter( P, C, 1, "points" );
    //pointPrinter.displayPDF();

    GeometricSpanner::GeometricSpannerSolutionChecker checker( GeometricSpanner::EPSILON );

    Random_points_in_square_2<Point,Creator> g1(a);
    std::copy_n( g1, n, back_inserter(P) ); // generate n random points in square of length size

    Greedy alg( P, C );

    checker.check( P, C );

    GeometricSpannerPrinter printer( P, C, 1, "pointsanddisks" );
    printer.displayPDF();
}

// conducts a single algorithm, single trial experiment
void simpleExperiment() {
    GeometricSpanner::simpleExperiment<
            Greedy,
            Random_points_in_disc_2< Point, Creator >
    >( 100, 5 );
}

// conducts a multi-algorithm, multi-trial experiment
void runExperiment() {

    GeometricSpanner::GeometricSpannerExperiment exp;

    std::vector<int> N {    // define num points for each trial
    //     10*thousand,
         50*thousand,
        100*thousand,
        500*thousand,
       1000*thousand
    };

    std::vector<int> A {    // define bounding box denominator
 //       thousand,           // note that N is divided by A to get bounding box size
        5*thousand,
        10*thousand,
        thousand,
        100
    };

    // Match N values with A values
    std::vector< std::pair< int, int > > squarePairs {
        std::make_pair( N[0], N[0]*2/A[0] ),
        std::make_pair( N[1], N[1]*2/A[1] ),
        std::make_pair( N[2], N[2]*2/A[2] )     // add more pairs if you want more point sets
    };
    // add square point sets to experiment
    GeometricSpanner::GeometricSpannerProblemSynthetic< CGAL::Random_points_in_square_2< Point, Creator > > squares(
        squarePairs,
        "square"
    );

    // circle pairs will be used for ellipses as well
    std::vector< std::pair< int, int > > circlePairs {
        std::make_pair( N[0], N[0]/A[0] ),
        std::make_pair( N[1], N[1]/A[1] ),
        std::make_pair( N[2], N[2]/A[2] )
    };
    GeometricSpanner::GeometricSpannerProblemSynthetic< CGAL::Random_points_in_disc_2< Point, Creator > > circles(
        circlePairs,
        "circle"
    );

    GeometricSpanner::GeometricSpannerProblemSynthetic< GeometricSpanner::Random_points_in_ellipse_2_hor< Point, Creator > > horEllipse(
        circlePairs,
        "horEllipse"
    );

    GeometricSpanner::GeometricSpannerProblemSynthetic< GeometricSpanner::Random_points_in_ellipse_2_vert< Point, Creator > > vertEllipse(
        circlePairs,
        "vertEllipse"
    );

     // circle pairs will be used for ellipses as well
    std::vector< std::pair< int, int > > onShapePairs {
        std::make_pair( N[0]/A[0], N[0]/A[0] ),
        std::make_pair( 2*N[1]/A[1], N[1]/A[1] ),
        std::make_pair( N[2]/A[2], N[2]/A[2] )
    };

    GeometricSpanner::GeometricSpannerProblemSynthetic< CGAL::Random_points_on_circle_2< Point, Creator > > onCircle(
        onShapePairs,
        "onCircle"
    );

    GeometricSpanner::GeometricSpannerProblemSynthetic< CGAL::Random_points_on_square_2< Point, Creator > > onSquare(
        onShapePairs,
        "onSquare"
    );

    // add problem from file
    //exp.addProblem( GeometricSpanner::GeometricSpannerProblem( "PIS-Small-1-0", "./experiments/problems/PIS-Small-1-0.txt" ) );
    //exp.addProblem( GeometricSpanner::GeometricSpannerProblem( "PID-Medium-0-5", "./experiments/problems/PID-Medium-0-5.txt" ) );

    // add synthetic problems
    exp.addProblem( squares );
    exp.addProblem( circles );
    //exp.addProblem( vertEllipse );
    exp.addProblem( horEllipse );
    exp.addProblem( onCircle );
    exp.addProblem( onSquare );

    // specify algorithms as template arguments
    exp.run<
        Greedy
    >(1);

}

