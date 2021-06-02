#include <chrono>
#include <fstream> //Reading and writing point sets.
#include <list>
#include <limits>
#include <optional>
#include <set>
#include <utility>

//Random point generation, testing.
#include <CGAL/point_generators_2.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

//#include "FloydWarshall.h"
//#include "GeometricSpannerPrinter.h"
////#include "GraphAlgoTV.h"
//#include "BGS2005.h"
//#include "LW2004.h"
//#include "BSX2009.h"
#include "KPX2010.h"
#include "BCC2012.h"
//#include "BHS2017.h"
#include "KPT2017.h"
#include "metrics.h"
//#include "delaunay.h"
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

    std::set<Point> points;

    std::copy_n( g2, n/9, inserter(points) );
    std::copy_n( g3, n/9, inserter(points) );
    std::copy_n( g4, n/18, inserter(points) );

    std::copy_n( g1s, n/9, inserter(points) );
    std::copy_n( g2s, n*2/9, inserter(points) );
    std::copy_n( g3s, n/9, inserter(points) );
    std::copy_n( g4s, n/18, inserter(points) );

    int remaining;
    while( (remaining = int(n)-int(points.size())) > 0) {
        std::copy_n( g1, remaining, inserter(points) );
    }


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

template< class OutputIterator >
string generatePointsNovel( OutputIterator pointsOut, size_t rows = 10, size_t cols = 10 ) {

    vector<Point> points;

    const double skew = 0.01;
    for( size_t i=0; i<rows; ++i ) {
        bool rowIsOdd = i%2;
        for( size_t j=rowIsOdd; j<cols; j+=1+(rowIsOdd) ) {
            bool colIsOdd = j%2;
            double y = i;
            y += (rowIsOdd || colIsOdd) ?
                0 : (skew * ( (i+j)%4 == 0 ? -1 : 1 ) );
//            if( rowIsEven && j%2 == 0 ) {
//                if( (i+j)%4 == 0 ) {
//                    y -= skew;
//                } else {
//                    y += skew;
//                }
//            }
            Point p(j,y);
            //cout<<p<<"\n";
            points.push_back(p);
        }
    }

    // copy points to output iterator
    for( Point p : points )
        *(pointsOut++) = p;

    // copy points to file
    ofstream out;
    string fName;
    fName = "data-NOVEL-" + to_string(points.size()) + "_" + to_string(rows) + "x" + to_string(cols) + ".txt";
    out.open( fName, ios::trunc );
    for( Point p : points )
        out << p << endl;

    out.close();

    return fName;
}

int main( int argc, char *argv[] ) {

    const size_t runs      =   100;
    const size_t n_begin   =  5000;
    const size_t n_end     = 10000;
    const size_t increment =  1000;

    vector<size_t> experimentParameters = {
        runs, n_begin, n_end, increment
    };

    for( size_t arg=1;
         arg < min(size_t(argc),experimentParameters.size()+1);
         ++arg )
    {
        try {
            experimentParameters[arg-1] = stoul(argv[arg]);
            cout<<"Parameter "<<(arg-1)<<" = "<<experimentParameters[arg-1]<<"\n";
        }
        catch(invalid_argument ia)
        {
            cout<<"Invalid experiment parameter '"<<arg<<"', using default value = "<<experimentParameters[arg-1]<<"\n";
        }
    }

    experiment( experimentParameters[0],experimentParameters[1],experimentParameters[2],experimentParameters[3] );
    //scratch();

    //singleRun( 0, 0, "kptTestResult", "data-75_4330.127019x4330.127019.txt", true, true );

    return 0;
}

void scratch() {
    using namespace std;

    const double width = 5;
    size_t i=30;
//
//    //for( i=1; i<=17; ++i ) {
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
        size_t n = 1000;

//        std::copy_n( g1, n/3, back_inserter(points) );
//        std::copy_n( g2, n/3, back_inserter(points) );
//        std::copy_n( g3, n/6, back_inserter(points) );
//        std::copy_n( g4, n, back_inserter(points) );
//        points.emplace_back( 0,0 );
//        points.emplace_back( 0,1 );
//        points.emplace_back( 1,0 );
//        points.emplace_back( 1,1 );
//        points.emplace_back( 2,1 );



        string filename = "data-89_4716.990566x4716.990566.txt";
        readPointsFromFile( back_inserter( points ), filename );

        //generatePointsNovel(back_inserter(points));

//        generateRandomPoints( n, width/2, back_inserter(points) );
        cout<< points.size();
        cout<< "\n";
//        list< pair< Point, Point > > result;
        list< pair< size_t,size_t > > result;

        // Delaunay triangulation
//        lw2004::Delaunay Del( points.begin(), points.end() );
//        // Add IDs
//        size_t id=0;
//        for( auto v=Del.finite_vertices_begin(); v!=Del.finite_vertices_end(); ++v )
//            v->info() = id++;

//                cout<<degree(Del);
//                cout<<",";
//                cout << weight( Del );
//                cout <<",";
//               Del.add_all_edges();
//                t = StretchFactor(Del);
//                cout<< t.second;
//                cout<<",";

//        {

//                Timer tim;
            //LW2004( points.begin(), points.end(), back_inserter(result) );
            //BSX2009( points.begin(), points.end(), back_inserter(result), 2*PI/3 );
            //BGS2005( points.begin(), points.end(), back_inserter(result) );
            //KPX2010( points.begin(), points.end(), back_inserter(result), 18 );
            //BCC2012<6>( points.begin(), points.end(), back_inserter(result) );
            //BCC2012<7>( points.begin(), points.end(), back_inserter(result) );
            //BHS2017(points.begin(), points.end(), back_inserter(result) );
            KPT2017(points.begin(), points.end(), back_inserter(result), true );
            //delaunay_testing( points.begin(), points.end(), back_inserter(result) );
        //}

//        for( auto edge : result ) {
//            cout<<edge.first<<"("<<points.at(edge.first)<<") -- "
//                <<edge.second<<"("<<points.at(edge.second)<<")\n";
//        }

//        list< pair< Point, Point > > WorstPath;
//        SFWorstPath( result.begin(), result.end(),
//                     make_optional(inserter(WorstPath,WorstPath.begin())) );


        cout << degree( result.begin(), result.end() );
        cout<<",";
        double t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout<< t;
        cout<<",";

//        cout << weight( result.begin(), result.end() )/2;
//        cout <<",";

//        GraphPrinter printer(9);
//        GraphPrinter::OptionsList options;
//
//        options = {
//            { "color", printer.inactiveEdgeColor },
//            { "line width", to_string(printer.inactiveEdgeWidth) }
//        };
//        printer.drawEdges( Del, options );
//
//        options = { // active edge options
//            { "color", printer.activeEdgeColor },
//            { "line width", to_string(printer.inactiveEdgeWidth) }
//        };
//        printer.drawEdges( result.begin(), result.end(), options );
//
//        options = { // worst path edge options
//            { "color", printer.worstPathEdgeColor },
//            { "line width", to_string(printer.activeEdgeWidth) }
//        };
//        printer.drawEdges( WorstPath.begin(), WorstPath.end(), options );
//
//
//        options = {
//            { "vertex", make_optional( to_string(printer.vertexRadius/3) ) }, // vertex width
//            { "color", make_optional( printer.backgroundColor ) }, // text color
//            { "fill", make_optional( printer.activeEdgeColor ) }, // vertex color
//            { "line width", make_optional( to_string(0) ) } // vertex border (same color as text)
//        };
//        GraphPrinter::OptionsList borderOptions = {
//            { "border", make_optional( to_string(printer.vertexRadius/2) ) }, // choose shape of vertex
//            { "color", printer.activeVertexColor }, // additional border color
//            { "line width", to_string(printer.inactiveEdgeWidth) }, // additional border width
//        };
//        printer.drawVertices( Del, options );
//
//        string outputFilename = "lw";
//        outputFilename += to_string(points.size());
//        printer.print( outputFilename );
//        cout<<"\n";
}

bool experiment( size_t trials, size_t n_start, size_t n_end, size_t increment ) {
    const double width = 1000;

    size_t invalid = 0;
    size_t trialNum = 1;

    for( size_t trial=0; trial<trials; ++trial ) {
        for( size_t n=n_start; n<=n_end; n+=increment ) {
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

//    list< pair< Point, Point > > result;
    list< pair< size_t, size_t > > result;
    size_t deg;

    // Delaunay triangulation
    //CGAL::Delaunay_triangulation_2<K> DT( points.begin(), points.end() );
//                cout << degree(DT);
//                cout << ",";
//                cout << weight(DT);
//                cout << ",";


//    cout<< points.size();
//    cout<< ",";
//    cout<< size;
//    cout<< ",";
//
//    {
//        Timer tim;
//        BGS2005( points.begin(), points.end(), back_inserter(result) );
//    }
//
//    deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
////
    double t;
//    t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );////    cout << t;
//    cout <<"\n";
////
//    result.clear();





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
//
//    t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
//    cout << t;
//    cout <<"\n";
//
//    result.clear();







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
////
//    {
//        Timer tim;
//        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
//    }
//    cout << t;
//    cout <<"\n";
//
//    result.clear();




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


    t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
    cout << t;
    cout <<"\n";

    result.clear();




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


    t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
    cout << t;
    cout <<"\n";

    result.clear();



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

    t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
    cout << t;
    cout <<"\n";

    result.clear();



//    cout<< points.size();
//    cout<< ",";
//    cout<< size;
//    cout<< ",";
//    {
//        Timer tim;
//        BHS2017( points.begin(), points.end(), back_inserter(result), printLog );
//    }
//
//    deg = degree( result.begin(), result.end() );
//    cout << deg;
//    cout <<",";
//
//    t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
//        cout << t;
//    cout <<"\n";
//
//    result.clear();





    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    {
        Timer tim;
        KPT2017( points.begin(), points.end(), back_inserter(result), printLog );
    }

    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    cout <<"\n";

    result.clear();






    if( deg > 4 || t > 20 || forcePrint ) {

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

        cout << "DEGREE ERROR!!! DEGREE:" << deg << "\n"<<endl;
        if( generatedFile )
            singleRun( n, width, resultFileName, *generatedFile, true, true );

        return false;
    }

    result.clear();

    cout<<"\n";

    return true;
}
