#include <chrono>
#include <fstream> // reading and writing point sets
#include <list>
#include <optional>
#include <utility>

#include <CGAL/point_generators_2.h>                            // Random point generation, testing

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "FloydWarshall.h"
#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "BGS2005.h"
#include "LW2004.h"
#include "BSX2009.h"
#include "KPX2010.h"
#include "metrics.h"
#include "utilities.h"

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

int main() {
//    size_t n = 20;
//    for( size_t i=4; i<=n; ++i )
//        if( !experiment( 1, i, i*1000, i*100 ) )
//            break;
    //singleRun( 0, 0, "bsxTestResult", "250_7905.694150x7905.694150.txt" );
    experiment( 5, 1000, 10000, 1000 );
    //scratch();
    //stretchScratch();
    //algoTVScratch();

    return 0;
}

void scratch() {
    using namespace std;

    const double width = 100;
    size_t n = 30, i=n;
//
//    //for( i=1; i<=17; ++i ) {
        auto g1 = CGAL::Random_points_in_square_2<Point,Creator>( width*sqrt(i)/2 );
        auto g2 = CGAL::Random_points_in_disc_2<Point,Creator>(   width*sqrt(i)/2 );
        auto g3 = CGAL::Random_points_on_square_2<Point,Creator>( width*sqrt(i)/2 );
        auto g4 = CGAL::Random_points_on_circle_2<Point,Creator>( width*sqrt(i)/2 );
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
//        std::copy_n( g3, n/2, back_inserter(points) );
//        std::copy_n( g4, n/2, back_inserter(points) );

//        points.emplace_back( 0,0 );
        string filename = "data-30_50.000000x50.000000.txt";

        readPointsFromFile<Point>( back_inserter( points ), filename );

        //generateRandomPoints( n, width/2, back_inserter(points) );


        n = points.size();

        cout<< n;
        cout<< ",";
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
            //BGS2002( points.begin(), points.end(), back_inserter(result) );
            KPX2010( points.begin(), points.end(), back_inserter(result), 18, true );
        }
        vector<vector<double>> apdShortestPaths;
        vector<vector<double>> apdiShortestPaths;
        vector<vector<double>> fwShortestPaths;
        AllPairsDijkstra( result.begin(), result.end(), back_inserter(apdShortestPaths) );
        AllPairsDijkstraIterative( result.begin(), result.end(), back_inserter(apdiShortestPaths) );
        FloydWarshall( result.begin(), result.end(), back_inserter(fwShortestPaths) );

        assert( apdShortestPaths.size() == fwShortestPaths.size() );
        for( size_t i=0; i<n; ++i ) {
            for( size_t j=0; j<n; ++j ) {
                cout<< i <<","<<j<<": apd="<<apdShortestPaths.at(i).at(j)<<", fw="<< fwShortestPaths.at(i).at(j) << "  ";
//                if( abs( apdShortestPaths.at(i).at(j) - fwShortestPaths.at(i).at(j) ) > EPSILON ) {
//                    cout<<"THEY DONT MATCH";
//                }
            }
            cout<<"\n\n";
        }

//        cout << degree( result.begin(), result.end() );
//        cout <<",";
//        cout << weight( result.begin(), result.end() )/2;
//        cout <<",";
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
        GraphPrinter printer(0.2);
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

        printer.print( "bsx2009" );
        cout<<"\n";
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
    size_t k = 14;
    // SET POINT SET
    list<Point> points;
    optional<string> generatedFile = nullopt;

    if( filename )
        readPointsFromFile<Point>( back_inserter( points ), *filename );
    else
        generatedFile = make_optional( generateRandomPoints<Point>( n, size, back_inserter(points) ) );

    n = points.size();
    cout<< n;
    cout<< ",";
    cout<< size;
    cout<< ",";

    list< pair< Point, Point > > result;



    // ALGORITHM TESTING AREA

    // Delaunay triangulation
    //CGAL::Delaunay_triangulation_2<K> DT( points.begin(), points.end() );
//                cout << degree(DT);
//                cout << ",";
//                cout << weight(DT);
//                cout << ",";


//
//
//    {
//        Timer tim;
//        BGS2005( points.begin(), points.end(), back_inserter(result) );
//    }
//    size_t deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
//
//    double t;
//    t = StretchFactor( result.begin(), result.end() );
//    cout << t;
//    cout <<",";
//
//    result.clear();
//
//
//
//
//    {
//        Timer tim;
//        LW2004( points.begin(), points.end(), back_inserter(result) );
//    }
//    deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
//
//    t = StretchFactor( result.begin(), result.end() );
//    cout << t;
//    cout <<",";
//
//    result.clear();
//
//
//
//
//
//
//
//    {
//        Timer tim;
//        BSX2009( points.begin(), points.end(), back_inserter(result) );
//    }
//    deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
//
//    t = StretchFactor( result.begin(), result.end() );
//    cout << t;
//    cout <<",";
//
//    result.clear();




    {
        Timer tim;
        KPX2010( points.begin(), points.end(), back_inserter(result), k, printLog );
    }
//    deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
//
//    t = StretchFactor( result.begin(), result.end() );
//    cout << t;
//    cout <<",";
//
//    result.clear();






    // ALL PAIRS TESTING AREA


//        vector<vector<double>> apdShortestPaths;
//        vector<vector<double>> apdiShortestPaths;
//        vector<vector<double>> apdpShortestPaths;
//        vector<vector<double>> fwShortestPaths;
//        vector<vector<double>> fwbfShortestPaths;

//
//        {
//            Timer tim;
//            AllPairsDijkstra( result.begin(), result.end(), back_inserter(apdShortestPaths) );
//        }
//
//        {
//            Timer tim;
//            AllPairsDijkstraIterative( result.begin(), result.end(), back_inserter(apdiShortestPaths) );
//        }

//        {
//            Timer tim;
//            AllPairsDijkstraIterativeParallel( result.begin(), result.end(), back_inserter(apdpShortestPaths) );
//        }
//
//        {
//            Timer tim;
//            FloydWarshall( result.begin(), result.end(), back_inserter(fwShortestPaths) );
//        }
//
//        {
//            Timer tim;
//            FloydWarshall_BF( result.begin(), result.end(), back_inserter(fwbfShortestPaths) );
//        }
//
//        for( size_t i=0; i<n; ++i ) {
//            for( size_t j=0; j<n; ++j ) {
//                if( abs( fwShortestPaths.at(i).at(j) - fwbfShortestPaths.at(i).at(j) ) > EPSILON ) {
//                    cout<<"FAIL";
//                    return false;
//                }
//            }
//        }


    // STRETCH FACTOR TESTING AREA

    double t_standard = INF;
    {
        Timer tim;
        t_standard = StretchFactor( result.begin(), result.end() );
    }
    cout<< t_standard <<",";

    double t_bffw = INF;
    {
        Timer tim;
        t_bffw = StretchFactorBFFW( result.begin(), result.end() );
    }
    cout<< t_bffw <<",";

    if( abs( t_standard - t_bffw ) > EPSILON ) {
        cout<<"FAIL";
        return false;
    }



        cout<<"PASS,";



//    if( t > 2.92 || deg > k || forcePrint ) {
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
//    result.clear();
//
//    {
//        Timer tim;
//        LW2004( points.begin(), points.end(), back_inserter(result), PI/2, printLog );
//    }
//    deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
//
//    t = 0;
//
//    {
//        //Timer tim;
//        t = StretchFactorDijkstraParallel( result.begin(), result.end() );
//    }
//    cout<< t;
//    cout<<",";
//
//    if( t > 7.79 || deg > 23 || forcePrint ) {
//        pair<pair<Vertex_handle,Vertex_handle>,double> t_fw;
//        t_fw = StretchFactor( result.begin(), result.end() );
//        cout<<t_fw.second<<",";
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
//            singleRun( n, width, resultFileName, *generatedFile, true, false );
//
//        return false;
//    }
    result.clear();

    cout<<"\n";

    return true;
}

