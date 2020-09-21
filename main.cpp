#include <chrono>
#include <list>
#include <optional>
#include <utility>

#include <CGAL/point_generators_2.h>                            // Random point generation, testing

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "Timer.h"
#include "FloydWarshall.h"
#include "GeometricSpannerPrinter.h"
#include "BGS2002.h"
//#include "LW2004_2.h"
#include "LW2004_3.h"
#include "metrics.h"

using namespace gsnunf;
typedef CGAL::Creator_uniform_2<double,Point> Creator;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

bool experiment(size_t,size_t,size_t,size_t=1);
bool singleRun(size_t,double,string,optional<string> = nullopt,bool=false,bool=false);
void scratch();
void stretchScratch();
void stretchFactorAndDegreeExperiment();

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
string generateRandomPoints( size_t n, double size, OutputIterator pointsOut ) {
    typedef CGAL::Creator_uniform_2<double,Point> Creator;

    auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( size );
    auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   size );
    auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( size );
    auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( size );

    vector<Point> points;
    points.reserve(n);

//    std::copy_n( g1, n/3, back_inserter(points) );
//    std::copy_n( g2, n/3, back_inserter(points) );
    std::copy_n( g3, n/2, back_inserter(points) );
    std::copy_n( g4, n/2, back_inserter(points) );

    points.emplace_back(0,0);

    // copy points to output iterator
    for( Point p : points )
        *(pointsOut++) = p;

    // copy points to file
    ofstream out;
    string fName;
    fName = to_string(n) + "_" + to_string(size) + "x" + to_string(size) + ".txt";
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
    experiment( 10, 1000, 50000, 1000 );
    //scratch();
    //stretchScratch();

    return 0;
}

void stretchScratch() {
    // 1. Create a graph. In this case it will be the DT of the point set
    list<Point> points;
    points = {
        { -1, 0.1 },
        { -0.9, 3 },
        { -2, 6 },
        { -7, 3.1 },
        { -6, -0.1 },
        { -9, -0.2 },
        { -7.7, -1 },
        { -6.1, -1.5 },
        { -10, -4 },
        { -4, -3 },
        { -1.5, -6 },
        { 1, -9 },
        { 4, -4 },
        { 4.1, 0 },
        { 3.9, 5.9 },
        { 5, 3 },
        { 5, -2 },
        { 9, 1 }
    };
    list< pair< Point, Point > > result;

    DelaunayGraph Del( points.begin(), points.end() );
    Del.add_all_edges();

    // 2. Create a subgraph of that graph
    LW2004_3( points.begin(), points.end(), back_inserter(result), PI/2, false );

    // measure stretch factor using Floyd Warshall (StretchFactor function)
    pair<pair<Vertex_handle,Vertex_handle>,double> t_fw = StretchFactor( result.begin(), result.end() );
    // measure stretch factor using experimental method
    double t_exp = StretchFactorDjikstra( result.begin(), result.end() );
    // print measurements
    cout<< points.size();
    cout<< ",";
    cout<< t_fw.second;
    cout<<",";
    cout<< t_exp;
    cout<<",";

    GraphPrinter printer(1);
    GraphPrinter::OptionsList options;

    options = {
        { "color", printer.inactiveEdgeColor },
        { "line width", to_string(printer.inactiveEdgeWidth) }
    };
    printer.drawEdges( Del._DT, options );

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
//    GraphPrinter::OptionsList borderOptions = {
//        { "border", make_optional( to_string(printer.vertexRadius) ) }, // choose shape of vertex
//        { "color", printer.activeEdgeColor }, // additional border color
//        { "line width", to_string(printer.inactiveEdgeWidth) }, // additional border width
//    };
    printer.drawVertices( Del._DT, options );

    printer.print( "stretchscratch" );
}

void scratch() {
    using namespace std;

    //GraphPrinter printer(0.05);
//
//    const double width = 100;
//    size_t n = 30, i=n;
//
//    //for( i=1; i<=17; ++i ) {
//        auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( width*sqrt(i)/2 );
//        auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   width*sqrt(i)/2 );
//        auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( width*sqrt(i)/2 );
//        auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( width*sqrt(i)/2 );
        list<Point> points;
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
//
//        n = 60;

//        std::copy_n( g1, n/3, back_inserter(points) );
//        std::copy_n( g2, n/3, back_inserter(points) );
//        std::copy_n( g3, n/6, back_inserter(points) );
        //std::copy_n( g4, n/6, back_inserter(points) );

        //points.emplace_back( 0,0 );
        string filename = "11_1658.312395x1658.312395.txt";

        readPointsFromFile( back_inserter( points ), filename );


        cout<< points.size();
        cout<< ",";
        list< pair< Point, Point > > result;
        pair<pair<Vertex_handle,Vertex_handle>,double> t;

        // Get t of Delaunay triangulation
                DelaunayGraph Del( points.begin(), points.end() );
//                cout<<degree(Del._DT);
//                cout<<",";
//                cout << weight( Del._DT );
//                cout <<",";
//               Del.add_all_edges();
//                t = StretchFactor(Del);
//                cout<< t.second;
//                cout<<",";

        {
//                Timer tim;
            LW2004_3( points.begin(), points.end(), back_inserter(result), PI/2, true );
            //BGS2002( points.begin(), points.end(), back_inserter(result) );
        }
//        cout << degree( result.begin(), result.end() );
//        cout <<",";
//        cout << weight( result.begin(), result.end() )/2;
//        cout <<",";
            t = StretchFactor( result.begin(), result.end() );
            cout<< t.second;
            cout<<",";
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
//        printer.print( "lw2004scratch" );
        //printer.print("bgs2002");


//        cout<<"\n";

        string resultFileName = filename;
        // strip file extension
        const std::string ext(".txt");
        if ( resultFileName != ext &&
             resultFileName.size() > ext.size() &&
             resultFileName.substr(resultFileName.size() - ext.size()) == ext )
        {
           // if so then strip them off
           resultFileName = resultFileName.substr(0, resultFileName.size() - ext.size());
        }
        resultFileName += "_result-";
        resultFileName += "redo";

        singleRun( 30, 30, resultFileName, filename, true, true );

}

bool experiment( size_t trials, size_t n_start, size_t n_end, size_t increment ) {
    const double width = 1000;

    size_t invalid = 0;

    for( size_t trial=0; trial<trials; ++trial ) {
        for( size_t n=n_start; n<=n_end; n+=increment ) {
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

    // SET POINT SET
    list<Point> points;
    optional<string> generatedFile = nullopt;

    if( filename )
        readPointsFromFile( back_inserter( points ), *filename );
    else
        generatedFile = make_optional( generateRandomPoints( n, size, back_inserter(points) ) );

    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";

    list< pair< Point, Point > > result;
    pair<pair<Vertex_handle,Vertex_handle>,double> t_floydwarshall;
    size_t deg = 0;

    // Delaunay triangulation
    CGAL::Delaunay_triangulation_2<K> DT( points.begin(), points.end() );
//                cout << degree(DT);
//                cout << ",";
//                cout << weight(DT);
//                cout << ",";

    {
        Timer tim;
        LW2004_3( points.begin(), points.end(), back_inserter(result), PI/2, printLog );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";
//            cout << weight( result.begin(), result.end() );
//            cout <<",";
//    {
//        Timer tim;
//        t_floydwarshall = StretchFactor( result.begin(), result.end() );
//    }
//    cout<< t_floydwarshall.second;
//    cout<<",";

    // measure stretch factor using experimental method
    double t_exp = 0;
    {
        Timer tim;
        t_exp = StretchFactorDjikstra( result.begin(), result.end() );
    }
    cout<< t_exp;
    cout<<",";
    {
        Timer tim;
        t_exp = StretchFactorDjikstraParallel( result.begin(), result.end() );
    }
    cout<< t_exp;
    cout<<",";
    //result.clear();

    if( t_floydwarshall.second > 10.01602 || deg > 23 || forcePrint ) {
        string resultFileName = ( filename ? *filename : *generatedFile );
        // strip file extension
        const std::string ext(".txt");
        if ( resultFileName != ext &&
             resultFileName.size() > ext.size() &&
             resultFileName.substr(resultFileName.size() - ext.size()) == ext )
        {
           // if so then strip them off
           resultFileName = resultFileName.substr(0, resultFileName.size() - ext.size());
        }
        resultFileName += "_result-";
        resultFileName += ( filename ? "redo" : "orig" );

//        cout << " Dumping edge set..."<<result.size()<<" edges.\n\n";
//        for( auto e : result )
//            cout <<"("<< e.first << ", " << e.second<<")" << "\n";


        if( generatedFile )
            singleRun( n, width, resultFileName, *generatedFile, true, true );

        return false;
    }
    result.clear();
//            {
//                Timer tim;
//                BGS2002( points.begin(), points.end(), back_inserter(result) );
//            }
//            cout << degree( result.begin(), result.end() );
//            cout << ",";
//            cout << weight( result.begin(), result.end() );
//            cout << ",";
//            t = StretchFactor( result.begin(), result.end() );
//            cout<< t.second;
//            cout<<",";
    cout<<"\n";

    return true;
}

