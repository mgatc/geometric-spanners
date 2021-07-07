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
#include "BGS2005.h"
//#include "BGS2005_2.h"
#include "LW2004.h"
#include "BSX2009.h"
#include "KPX2010.h"
#include "KX2012.h"
#include "BCC2012.h"
#include "BHS2017.h"
#include "KPT2017.h"
#include "BKPX2015.h"
#include "BGHP2010.h"
#include "metrics.h"
//#include "delaunay.h"
#include "utilities.h"


using namespace unf_planespanners;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Creator_uniform_2<double,Point> Creator;

bool experiment(size_t,size_t,size_t,size_t=1);
bool singleRun(size_t,double,string,optional<string> = nullopt,bool=false,bool=false);
void scratch(size_t n1);
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
    std::copy_n( g4, n/9, inserter(points) );

    std::copy_n( g1s, n/9, inserter(points) );
    std::copy_n( g2s, n*2/9, inserter(points) );
    std::copy_n( g3s, n/9, inserter(points) );
    std::copy_n( g4s, n/18, inserter(points) );

    int remaining;
    while( (remaining = int(n)-int(points.size())) > 0) {
        std::copy_n( g1, remaining, inserter(points) );
    }


    //points.emplace(0,0);

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
        catch(invalid_argument &ia)
        {
            cout<<"Invalid experiment parameter '"<<arg<<"', using default value = "<<experimentParameters[arg-1]<<"\n";
        }
    }

    size_t N = experimentParameters.empty() ? 50 : experimentParameters[0];
    experiment( experimentParameters[0],experimentParameters[1],experimentParameters[2],experimentParameters[3] );
    //scratch(N);

 //   singleRun( 0, 0, "BKPXTestResult", "positive_points.txt", true, true );


    return 0;
}

void scratch(size_t n) {
    using namespace std;

    vector<Point> points;

    /////////// GET POINTS SOMEHOW //////////////////////////

    string filename = "data-89_4716.990566x4716.990566.txt";
    //readPointsFromFile( back_inserter( points ), filename );

    double width = 5;
    generateRandomPoints( n, width/2, back_inserter(points) );

    /////////////////////////////////////////////////////////

    cout<< points.size();
    cout<< "\n";

    list< pair< size_t,size_t > > result;

    { // RUN THE ALGORITHM(S) /////////////////////////////////////
//                Timer tim;
            //LW2004( points.begin(), points.end(), back_inserter(result) );
            //BSX2009( points.begin(), points.end(), back_inserter(result), 2*PI/3 );
            //BGS2005( points.begin(), points.end(), back_inserter(result) );
            //KPX2010( points.begin(), points.end(), back_inserter(result), 18 );
            //BCC2012<6>( points.begin(), points.end(), back_inserter(result) );
            //BCC2012<7>( points.begin(), points.end(), back_inserter(result) );
            //BHS2017(points.begin(), points.end(), back_inserter(result) );
            //KPT2017(points.begin(), points.end(), back_inserter(result), true );
            //BKPX2015(points.begin(), points.end(), back_inserter(result), true );
          //  BGHP2010(points.begin(), points.end(), back_inserter(result), true );
            KX2012( points.begin(), points.end(), back_inserter(result), true);
            //delaunay_testing( points.begin(), points.end(), back_inserter(result) );
    }

//        for( auto edge : result ) {
//            cout<<edge.first<<"("<<points.at(edge.first)<<") -- "
//                <<edge.second<<"("<<points.at(edge.second)<<")\n";
//        }

//        list< pair< Point, Point > > WorstPath;
//        SFWorstPath( result.begin(), result.end(),
//                     make_optional(inserter(WorstPath,WorstPath.begin())) );



        cout << degree( result.begin(), result.end() )<<endl;
        // cout<<",";
        // double t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        // cout<< t;
        // cout<<",";


    // PRODUCE A LaTeX / TiKz DOCUMENT AND DISPLAY

    GraphPrinter printer;

    double documentSizeInCm = 20;

    printer.autoscale( points.begin(), points.end(), documentSizeInCm );
    GraphPrinter::OptionsList options;

    options = { // active edge options
        {"color", printer.activeEdgeColor},
        {"line width", to_string(printer.activeEdgeWidth)}
    };

    printer.drawEdges(result.begin(), result.end(), points, options);

    options = {
        {"vertex", make_optional(to_string(printer.vertexRadius))}, // vertex width
        {"color", make_optional(printer.backgroundColor)}, // text color
        {"fill", make_optional(printer.activeVertexColor)}, // vertex color
        {"line width", make_optional(to_string(0))} // vertex border (same color as text)
    };
    GraphPrinter::OptionsList borderOptions = {
        {"border", make_optional(to_string(printer.vertexRadius))}, // choose shape of vertex
        {"color", printer.activeEdgeColor}, // additional border color
        {"line width", to_string(printer.inactiveEdgeWidth)}, // additional border width
    };
    printer.drawVerticesWithInfo(points.begin(), points.end(), options, borderOptions);

    printer.print("test");
    cout << "\n";

}

bool experiment( size_t trials, size_t n_start, size_t n_end, size_t increment ) {
    const double width = 10;

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
    //cout<<"\nTesting Complete. "<< invalid << " of "<<(trials*(n_end-n_start))<<" invalid results.\n\n";
    return invalid == 0;
}

bool singleRun( size_t n, double width, string resultFilename, optional<string> filename, bool forcePrint, bool printLog )
{
    double size = width/2; // cgal's generators produce width 2x given value

    // SET POINT SET
    list<Point> points;
    optional<string> generatedFile = nullopt;

    if( filename )
        readPointsFromFile( back_inserter( points ), *filename );
    else
        generatedFile = make_optional( generateRandomPoints( n, size, back_inserter(points) ) );

    bool measureStretchFactor = true;

//    list< pair< Point, Point > > result;
    list< pair< size_t, size_t > > result;
    size_t deg;
//
    double t;

    // STANDARD L2 DELAUNAY TRIANGULATION BASED ALGORITHMS ///////////////

    //CGAL::Delaunay_triangulation_2<K> DT( points.begin(), points.end() );
//                cout << degree(DT);
//                cout << ",";
//                cout << weight(DT);
//                cout << ",";


    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout<< "BGS2005,";

    {
        Timer tim;
        BGS2005( points.begin(), points.end(), back_inserter(result) );
    }

    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";

    result.clear();




    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout << "LW2004,   ";
    {
        Timer tim;
        LW2004( points.begin(), points.end(), back_inserter(result) );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";
    result.clear();



    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout<< "BSX2009,";
    {
        Timer tim;
        BSX2009( points.begin(), points.end(), back_inserter(result) );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";
    result.clear();



    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout<< "KPX2010,";
    const size_t k = 14; // the max degree of the spanner such that k>=14
    {
        Timer tim;
        KPX2010( points.begin(), points.end(), back_inserter(result), k, printLog );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";
    result.clear();



    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout<< "KX2012,";
    {
        Timer tim;
        KX2012( points.begin(), points.end(), back_inserter(result), printLog );
    }

    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";
    result.clear();




    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout<< "BCC-7,";
    {
        Timer tim;
        BCC2012<7>( points.begin(), points.end(), back_inserter(result), printLog );
    }

    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";
    result.clear();




    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout<< "BCC-6,";
    {
        Timer tim;
        BCC2012<6>( points.begin(), points.end(), back_inserter(result), printLog );
    }

    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";
    result.clear();




    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout<< "BHS2017,";
    {
        Timer tim;
        BHS2017( points.begin(), points.end(), back_inserter(result), printLog );
    }

    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";
    result.clear();




    // TD-DELAUNAY / HALF-THETA-6 GRAPH BASED ALGORITHMS /////////////

    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout<< "BGHP2010,";
    {
        Timer tim;
        BGHP2010( points.begin(), points.end(), back_inserter(result), printLog );
    }

    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";
    result.clear();



    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout<< "KPT2017,";
    {
        Timer tim;
        KPT2017( points.begin(), points.end(), back_inserter(result), printLog );
    }

    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";
    result.clear();



    // Linf DELAUNAY GRAPH BASED ALGORITHMS /////////////////////////////

    cout<< points.size();
    cout<< ",";
    cout<< size;
    cout<< ",";
    cout << "BKPX2015, ";
    {
        Timer tim;
        BKPX2015( points.begin(), points.end(), back_inserter(result), printLog );
    }
    deg = degree( result.begin(), result.end() );
    cout << deg;
    cout <<",";

    if(measureStretchFactor){
        t = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
        cout << t;
    }
    cout <<"\n";
    result.clear();






//    if( deg > 11 || t > 7 || forcePrint ) {
//
//       string resultFileName = ( filename ? *filename : *generatedFile );
//       // strip file extension
//       const std::string ext(".txt");
//       if ( resultFileName != ext &&
//            resultFileName.size() > ext.size() &&
//            resultFileName.substr(resultFileName.size() - ext.size()) == ext )
//       {
//          // if so then strip them off
//          resultFileName = resultFileName.substr(0, resultFileName.size() - ext.size());
//       }
//       resultFileName += "_result-";
//       resultFileName += ( filename ? "redo" : "orig" );
//
//       cout << "DEGREE ERROR!!! DEGREE:" << deg << "\n"<<endl;
//       cout << *generatedFile <<endl;
//
//       if( generatedFile )
//           singleRun( n, width, resultFileName, *generatedFile, true, true );
//
//       return false;
//    }

    cout<<endl;

    return true;
}
