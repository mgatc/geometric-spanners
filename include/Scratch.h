//
// Created by matt on 7/21/21.
//

#ifndef GEOMETRIC_SPANNERS_SCRATCH_H
#define GEOMETRIC_SPANNERS_SCRATCH_H



#include <iostream>
#include <iterator>
#include <list>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "LatexPrinter.h"
//#include "LineGraphPrinter.h"
#include "metrics.h"
#include "names.h"
#include "Result.h"
#include "utilities.h"


#include "BCC2012.h"
#include "BGHP2010.h"
#include "BGS2005.h"
#include "BHS2017.h"
#include "BKPX2015.h"
#include "BSX2009.h"
#include "KPT2017.h"
#include "KPX2010.h"
#include "KX2012.h"
#include "LW2004.h"

namespace unf_spanners {

void scratch(size_t n) {
    using namespace std;

    vector<Point> points;

    /////////// GET POINTS SOMEHOW //////////////////////////

    string filename = "data-89_4716.990566x4716.990566.txt";
    //readPointsFromFile( back_inserter( points ), filename );

    double width = 10;
    generateRandomPoints(n, width / 2, back_inserter(points));

    /////////////////////////////////////////////////////////

    cout << points.size();
    cout << "\n";

    list<pair<size_t, size_t> > result;

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
        KX2012(points.begin(), points.end(), back_inserter(result), true);
        //delaunay_testing( points.begin(), points.end(), back_inserter(result) );
    }

//        for( auto edge : result ) {
//            cout<<edge.first<<"("<<points.at(edge.first)<<") -- "
//                <<edge.second<<"("<<points.at(edge.second)<<")\n";
//        }

//        list< pair< Point, Point > > WorstPath;
//        SFWorstPath( result.begin(), result.end(),
//                     make_optional(inserter(WorstPath,WorstPath.begin())) );



    cout << degree(result.begin(), result.end()) << endl;
    // cout<<",";
    // double stretchFactor = StretchFactorDijkstraReduction( points.begin(), points.end(), result.begin(), result.end() );
    // cout<< stretchFactor;
    // cout<<",";


    // PRODUCE A LaTeX / TiKz DOCUMENT AND DISPLAY

    GraphPrinter tikz("temp-graph");

    tikz.autoscale(points.begin(), points.end());
    LatexPrinter::OptionsList options;

    options = { // active edge options
            {"color",      tikz.activeEdgeColor},
            {"line width", to_string(tikz.activeEdgeWidth)}
    };

    tikz.drawEdges(result.begin(), result.end(), points, options);

    options = {
            {"vertex",     to_string(tikz.vertexRadius)}, // vertex width
            {"color",      tikz.backgroundColor}, // text color
            {"fill",       (tikz.activeVertexColor)}, // vertex color
            {"line width", (to_string(0))} // vertex border (same color as text)
    };
    LatexPrinter::OptionsList borderOptions = {
            {"border",     to_string(tikz.vertexRadius)}, // choose shape of vertex
            {"color",      tikz.activeEdgeColor}, // additional border color
            {"line width", to_string(tikz.inactiveEdgeWidth)} // additional border width
    };
    tikz.drawVerticesWithInfo(points.begin(), points.end(), options, borderOptions);


    LatexPrinter latex("temp-latex");

    latex.addToDocument(tikz);

    latex.display();


//    graph.print("test");
    cout << "\n";

}

} // unf_spanners

#endif //GEOMETRIC_SPANNERS_SCRATCH_H
