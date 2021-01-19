#include <chrono>
#include <fstream> //Reading and writing point sets.
#include <list>
#include <limits>
#include <optional>
#include <utility>

//Random point generation, testing.
#include <CGAL/point_generators_2.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "FloydWarshall.h"
#include "GeometricSpannerPrinter.h"
#include "GraphAlgoTV.h"
#include "BGS2005.h"
#include "LW2004.h"
#include "BSX2009.h"
#include "KPX2010.h"
#include "BCC2012.h"
#include "BHS2017.h"
#include "metrics.h"
//#include "utilities.h"

using namespace gsnunf;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Creator_uniform_2<double,Point> Creator;

bool experiment(size_t,size_t,size_t,size_t=1);
bool singleRun(size_t,double,string,optional<string> = nullopt,bool=false,bool=false);
void scratch();
void stretchScratch();
void stretchFactorAndDegreeExperiment();
void algoTVScratch();

template< class OutputIterator >
void readPointsFromFile( OutputIterator out, const string outputFileName, const size_t n=SIZE_T_MAX ) {
    ifstream in(outputFileName);
    if (in.is_open()) {
        double x,y;
        size_t i = 0;
        while ( i<n && in >> x >> y ) {
            *out = Point(x,y);
            ++out;
            ++i;
        }
        in.close();
    }
}

template< class OutputIterator >
string generateRandomPoints( size_t n, double size, OutputIterator pointsOut ) {
    typedef CGAL::Creator_uniform_2<double,Point> Creator;

    auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( size );
    auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   size );
    auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( size );
    auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( size );


    auto g1s = CGAL::Random_points_in_square_2<Point,Creator>( size/4 );
    auto g2s = CGAL::Random_points_in_disc_2<Point,Creator>(   size/4 );
    auto g3s = CGAL::Random_points_on_square_2<Point,Creator>( size/4 );
    auto g4s = CGAL::Random_points_on_circle_2<Point,Creator>( size/4 );

    vector<Point> points;
    points.reserve(n);

    std::copy_n( g1, n*2/9, back_inserter(points) );
    std::copy_n( g2, n/9, back_inserter(points) );
    std::copy_n( g3, n*2/18, back_inserter(points) );
    std::copy_n( g4, n/18, back_inserter(points) );

    std::copy_n( g1s, n/9, back_inserter(points) );
    std::copy_n( g2s, n*2/9, back_inserter(points) );
    std::copy_n( g3s, n/18, back_inserter(points) );
    std::copy_n( g4s, n*2/18, back_inserter(points) );

    //points.emplace_back(0,0);

    // copy points to output iterator
    for( Point p : points )
        *(pointsOut++) = p;

    // copy points to file
    ofstream out;
    string fName;
    fName = "data-" + to_string(n) + "_" + to_string(size) + "x" + to_string(size) + ".txt";
    out.open( fName, ios::trunc );
    for( Point p : points )
        out << p << endl;

    out.close();

    return fName;
}

int main() {
//    size_t n = 20;
//    for( size_t i=4; i<=n; ++i )
//        if( !experiment( 1, i, i*1000, i*100 ) )
//            break;

<<<<<<< HEAD
    //singleRun( 0, 0, "bhsTestResult", "data-100_5000.000000x5000.000000.txt" );
    experiment( 100000, 100000, 1000000, 10000 );
=======
    //singleRun( 0, 0, "bsxTestResult", "data-200_7071.067812x7071.067812.txt", true, true );
    //singleRun( 0, 0, "bsxTestResult", "data-200_7071.067812x7071.067812 (copy).txt", true, true );
    experiment( 100000, 100000, 2000000, 100000 );

>>>>>>> e44eb03bd22a93ae9c103039942c4a1909f89eb1
    //scratch();

    return 0;
}

void scratch() {
    using namespace std;

    const double width = 5;
    size_t i=30;
//
//    //for( i=1; i<=17; ++i ) {
        auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( width*sqrt(i)/2 );
        auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   width*sqrt(i)/2 );
        auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( width*sqrt(i)/2 );
        auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( width*sqrt(i)/2 );
        vector<Point> points;
        // SET POINT SET
//        points = {
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
//    points = {
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

//    points = {
//        {
//            2.73,4.13
//        },
//        {
//            2.87,0
//        },
//        {
//            3.5,-1.4
//        },
//        {
//            -6.3,-0.14
//        },
//        {
//            -4.27,-1.05
//        },
//        {
//            -2.8,-2.1
//        },
//        {
//            -7.2,-2.8
//        },
//        {
//            -5.39,-0.7
//        },
//        {
//            -4.2,-0.07
//        },
//        {
//            -1.4,4.2
//        },
//        {
//            -4.9,2.17
//        },
//        {
//            6.3,0.7
//        },
//        {
//            3.5,2.1
//        },
//        {
//            -0.63,2.1
//        },
//        {
//            -0.7,0.07
//        },
//        {
//            2.8,-2.8
//        },
//        {
//            0.7,-6.3
//        },
//        {
//            -1.05,-4.2
//        },
//    };
        size_t n = 40;

//        std::copy_n( g1, n/3, back_inserter(points) );
//        std::copy_n( g2, n/3, back_inserter(points) );
//        std::copy_n( g3, n/6, back_inserter(points) );
//        std::copy_n( g4, n, back_inserter(points) );
//        points.emplace_back( 0,0 );
//        points.emplace_back(1, 13 );


        string filename = "data-30_2738.612788x2738.612788.txt";
        readPointsFromFile( back_inserter( points ), filename );

        //generateRandomPoints( n, width/2, back_inserter(points) );
        cout<< points.size();
        cout<< "\n";
        list< pair< Point, Point > > result;

        // Delaunay triangulation
        lw2004::Delaunay Del( points.begin(), points.end() );
        // Add IDs
        size_t id=0;
        for( auto v=Del.finite_vertices_begin(); v!=Del.finite_vertices_end(); ++v )
            v->info() = id++;

//                cout<<degree(Del);
//                cout<<",";
//                cout << weight( Del );
//                cout <<",";
//               Del.add_all_edges();
//                t = StretchFactor(Del);
//                cout<< t.second;
//                cout<<",";

        {

//                Timer tim;
            //LW2004_3( points.begin(), points.end(), back_inserter(result), PI/2, true );
            //BSX2009( points.begin(), points.end(), back_inserter(result), 2*PI/3, true );
            //BGS2002( points.begin(), points.end(), back_inserter(result) );
            //KPX2010( points.begin(), points.end(), back_inserter(result), 18, true );
            BCC2012<6>( points.begin(), points.end(), back_inserter(result), true );

        }

       //Johnsons( result.begin(), result.end() );

        cout << degree( result.begin(), result.end() );
        cout <<",";

//        cout << weight( result.begin(), result.end() )/2;
//        cout <<",";

        // Add edges to result

//        result = {
//            {
//                points[0],points[9]
//            },
//            {
//                points[0],points[12]
//            },
//            {
//                points[0],points[16]
//            },
//            {
//                points[1],points[2]
//            },
//            {
//                points[1],points[12]
//            },
//            {
//                points[1],points[13]
//            },
//            {
//                points[1],points[8]
//            },
//            {
//                points[2],points[11]
//            },
//            {
//                points[2],points[14]
//            },
//            {
//                points[2],points[15]
//            },
//            {
//                points[3],points[7]
//            },
//            {
//                points[3],points[8]
//            },
//            {
//                points[3],points[10]
//            },
//            {
//                points[4],points[5]
//            },
//            {
//                points[4],points[6]
//            },
//            {
//                points[4],points[7]
//            },
//            {
//                points[4],points[8]
//            },
//            {
//                points[5],points[6]
//            },
//            {
//                points[5],points[8]
//            },
//            {
//                points[5],points[13]
//            },
//            {
//                points[5],points[17]
//            },
//            {
//                points[5],points[14]
//            },
//            {
//                points[6],points[7]
//            },
//            {
//                points[6],points[16]
//            },
//            {
//                points[8],points[10]
//            },
//            {
//                points[8],points[12]
//            },
//            {
//                points[8],points[13]
//            },
//            {
//                points[9],points[13]
//            },
//            {
//                points[9],points[10]
//            },
//            {
//                points[11],points[12]
//            },
//            {
//                points[11],points[15]
//            },
//            {
//                points[12],points[13]
//            },
//            {
//                points[13],points[14]
//            },
//            {
//                points[14],points[15]
//            },
//            {
//                points[15],points[16]
//            },
//            {
//                points[16],points[17]
//            },
//        };
//
//        size_t n_p = result.size();
//        list<pair<Point,Point>> result2(result);
//
//        for( auto it=result.begin();it!=result.end();it++ ) {
//            result2.emplace_back(it->second, it->first);
//        }

//            double t = StretchFactorDijkstraReduction( result.begin(), result.end() );
//            cout<< t;
//            cout<<",";
//            size_t deg = degree( result.begin(), result.end() );
//            cout << deg;
//            cout <<",";
//            cout << " Dumping edge set..."<<result.size()<<" edges.\n\n";


//            for( auto e : result ) cout <<"("<< e.first << ", " << e.second<<")" << "\n";
//        result.clear();

//            {
//                Timer tim;
//                BGS2002( points.begin(), points.end(), back_inserter(result) );
//            }
//            t = StretchFactor( result.begin(), result.end() );
//            cout<< t.second;
//            cout<<",";
//            result.clear();

//        cout<<"\n";

//        printer.drawEdges( Del._DT );
//        printer.drawEdges( result.begin(), result.end(), {{"red",""}} );
//        printer.drawVertices( Del._DT );
//        printer.print( "kpx2010scratch" );
        //printer.print("bgs2002");


//        cout<<"\n";

//        string resultFileName = filename;
//        // strip file extension
//        const std::string ext(".txt");
//        if ( resultFileName != ext &&
//             resultFileName.size() > ext.size() &&
//             resultFileName.substr(resultFileName.size() - ext.size()) == ext )
//        {
//           // if so then strip them off
//           resultFileName = resultFileName.substr(0, resultFileName.size() - ext.size());
//        }
//        resultFileName += "_result-";
//        resultFileName += "redo";

       // singleRun( 30, 30, resultFileName, filename, true, true );
        GraphPrinter printer(0.7);
        GraphPrinter::OptionsList options;

        options = {
            { "color", printer.inactiveEdgeColor },
            { "line width", to_string(printer.inactiveEdgeWidth) }
        };
        printer.drawEdges( Del, options );

        options = { // active edge options
            { "color", printer.activeEdgeColor },
            { "line width", to_string(printer.activeEdgeWidth) }
        };
        printer.drawEdges( result.begin(), result.end(), options );


        options = {
            { "vertex", make_optional( to_string(printer.vertexRadius) ) }, // vertex width
            { "color", make_optional( printer.backgroundColor ) }, // text color
            { "fill", make_optional( printer.activeVertexColor ) }, // vertex color
            { "line width", make_optional( to_string(0) ) } // vertex border (same color as text)
        };
        GraphPrinter::OptionsList borderOptions = {
            { "border", make_optional( to_string(printer.vertexRadius) ) }, // choose shape of vertex
            { "color", printer.activeEdgeColor }, // additional border color
            { "line width", to_string(printer.inactiveEdgeWidth) }, // additional border width
        };
        printer.drawVerticesWithInfo( Del, options, borderOptions );

        string outputFilename = "rand";
        outputFilename += to_string(points.size());
        printer.print( outputFilename );
        cout<<"\n";

}

bool experiment( size_t trials, size_t n_start, size_t n_end, size_t increment ) {
    const double width = 1000;

    size_t invalid = 0;
    size_t trialNum = 1;

    for( size_t trial=0; trial<trials; ++trial ) {
        for( size_t n=n_start; n<=n_end; n+=increment ) {
            cout << "Trial:" << trialNum << "\n";
            trialNum++;
            if( !singleRun( n, width*sqrt(n), "output", nullopt, false, false ) ) {
                ++invalid;
                return false;
            }
        }
    }
    //cout<<"\nTesting complete. "<< invalid << " of "<<(trials*(n_end-n_start))<<" invalid results.\n\n";
    return invalid == 0;
}

bool singleRun( size_t n, double width, string resultFilename, optional<string> filename, bool forcePrint, bool printLog ) {
    double size = width/2; // cgal's generators produce width 2x given value
    size_t k = 14;
    // SET POINT SET
    list<Point> points;
    optional<string> generatedFile = nullopt;

    if( filename )
        readPointsFromFile( back_inserter( points ), *filename );
    else
        generatedFile = make_optional( generateRandomPoints( n, size, back_inserter(points) ) );

//    cout<< points.size();
//    cout<< ",";
//    cout<< size;
//    cout<< ",";

    list< pair< Point, Point > > result;
    size_t deg;

////    // Delaunay triangulation
////    //CGAL::Delaunay_triangulation_2<K> DT( points.begin(), points.end() );
//////                cout << degree(DT);
//////                cout << ",";
//////                cout << weight(DT);
//////                cout << ",";
////
////
////
////
//    {
//        Timer tim;
//        BGS2005( points.begin(), points.end(), back_inserter(result) );
//    }
//
//    deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
//
<<<<<<< HEAD
//    double t;
////    t = StretchFactorDijkstraReduction( result.begin(), result.end() );
////    cout << t;
//    cout <<"\n";
//
//    result.clear();
//
//
//
//
//
//    cout<< points.size();
//    cout<< ",";
//    cout<< size;
//    cout<< ",";
//    {
//        Timer tim;
//        LW2004( points.begin(), points.end(), back_inserter(result) );
//    }
//    deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
////
////    t = StretchFactorDijkstraReduction( result.begin(), result.end() );
////    cout << t;
//    cout <<"\n";
////
//    result.clear();
////
////
////
////
////
////
////
//    cout<< points.size();
//    cout<< ",";
//    cout<< size;
//    cout<< ",";
//    {
//        Timer tim;
//        BSX2009( points.begin(), points.end(), back_inserter(result) );
//    }
//    deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
//
////    t = StretchFactorDijkstraReduction( result.begin(), result.end() );
////    cout << t;
//    cout <<"\n";
//
//    result.clear();
//
//
//
//
//    cout<< points.size();
//    cout<< ",";
//    cout<< size;
//    cout<< ",";
//    {
//        Timer tim;
//        KPX2010( points.begin(), points.end(), back_inserter(result), k, printLog );
//    }
//    deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
////
////
////    t = StretchFactorDijkstraReduction( result.begin(), result.end() );
////    cout << t;
//    cout <<"\n";
////
//    result.clear();


=======
    {
        Timer tim;
        BGS2005( points.begin(), points.end(), back_inserter(result) );
    }

    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";


    double t;
//    t = StretchFactorDijkstraReduction( result.begin(), result.end() );
//    cout << t;
    cout <<"\n";

    result.clear();






    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    {
        Timer tim;
        LW2004( points.begin(), points.end(), back_inserter(result) );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

//    t = StretchFactorDijkstraReduction( result.begin(), result.end() );
//    cout << t;
    cout <<"\n";

    result.clear();




    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    {
        Timer tim;
        BSX2009( points.begin(), points.end(), back_inserter(result) );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

//    t = StretchFactorDijkstraReduction( result.begin(), result.end() );
//    cout << t;
    cout <<"\n";

    result.clear();




    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    {
        Timer tim;
        KPX2010( points.begin(), points.end(), back_inserter(result), k, printLog );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";
//
//
//    t = StretchFactorDijkstraReduction( result.begin(), result.end() );
//    cout << t;
    cout <<"\n";
>>>>>>> e44eb03bd22a93ae9c103039942c4a1909f89eb1

    result.clear();




//    result.clear();

    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    {
        Timer tim;
        BCC2012<7>( points.begin(), points.end(), back_inserter(result), printLog );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";


//    double t = StretchFactorDijkstraReduction( result.begin(), result.end() );
//    cout << t;
    cout <<"\n";

    result.clear();



<<<<<<< HEAD
    if( deg > 8 || forcePrint ) {
=======

    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    {
        Timer tim;
        BCC2012<6>( points.begin(), points.end(), back_inserter(result), printLog );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";


//    t = StretchFactorDijkstraReduction( result.begin(), result.end() );
//    cout << t;
    cout <<"\n";

    result.clear();
>>>>>>> e44eb03bd22a93ae9c103039942c4a1909f89eb1



    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    {
        Timer tim;
        BHS2017( points.begin(), points.end(), back_inserter(result), printLog );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";


//    t = StretchFactorDijkstraReduction( result.begin(), result.end() );
//    cout << t;
    cout <<"\n";

    result.clear();






//    if( deg > 6 || forcePrint ) {
//
//        string resultFileName = ( filename ? *filename : *generatedFile );
//        // strip file extension
//        const std::string ext(".txt");
//        if ( resultFileName != ext &&
//             resultFileName.size() > ext.size() &&
//             resultFileName.substr(resultFileName.size() - ext.size()) == ext )
//        {
//           // if so then strip them off
//           resultFileName = resultFileName.substr(0, resultFileName.size() - ext.size());
//        }
//        resultFileName += "_result-";
//        resultFileName += ( filename ? "redo" : "orig" );
//
//        if( generatedFile )
//            singleRun( n, width, resultFileName, *generatedFile, true, true );
//
//        return false;
//    }
//
//    result.clear();

    cout<<"\n";

    return true;
}
