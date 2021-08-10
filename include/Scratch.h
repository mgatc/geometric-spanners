#ifndef PLANESPANNERS_SCRATCH_H
#define PLANESPANNERS_SCRATCH_H



#include <iostream>
#include <iterator>
#include <list>
#include <optional>
#include <string>
#include <utility>
#include <vector>



#include "algorithms/BCC2012.h"
#include "algorithms/BGHP2010.h"
#include "algorithms/BGS2005.h"
#include "algorithms/BHS2018.h"
#include "algorithms/BKPX2015.h"
#include "algorithms/BSX2009.h"
#include "algorithms/KPT2017.h"
#include "algorithms/KPX2010.h"
#include "algorithms/KX2012.h"
#include "algorithms/LW2004.h"

#include "Names.h"

#include "printers/LatexPrinter.h"
#include "printers/PgfplotPrinter.h"
#include "printers/TablePrinter.h"

#include "tools/Metrics.h"
#include "tools/PointGenerators.h"
#include "tools/Results.h"
#include "tools/Utilities.h"

namespace planespanners {


    void scratch(const vector<Point>& points) {

        using namespace std;
        cout << points.size();
        cout << "\n";

        list<pair<size_t, size_t> > result;

        { // RUN THE ALGORITHM(S) /////////////////////////////////////
//            Timer tim;
//            LW2004( points.begin(), points.end(), back_inserter(result) );
//            BSX2009( points.begin(), points.end(), back_inserter(result), 2*PI/3 );
//            BGS2005( points.begin(), points.end(), back_inserter(result) );
//            KPX2010( points.begin(), points.end(), back_inserter(result), 18 );
//            BCC2012<6>( points.begin(), points.end(), back_inserter(result) );
//            BCC2012<7>( points.begin(), points.end(), back_inserter(result) );
//            BHS2018(points.begin(), points.end(), back_inserter(result) );
//            KPT2017(points.begin(), points.end(), back_inserter(result), true );
            BKPX2015(points.begin(), points.end(), back_inserter(result));
//            BGHP2010(points.begin(), points.end(), back_inserter(result), true );
//            KX2012(points.begin(), points.end(), back_inserter(result), true);
//            delaunay_testing( points.begin(), points.end(), back_inserter(result) );
        }

//        for( auto edge : result ) {
//            cout<<edge.first<<"("<<points.at(edge.first)<<") -- "
//                <<edge.second<<"("<<points.at(edge.second)<<")\n";
//        }

//        list< pair< Point, Point > > WorstPath;
//        SFWorstPath( result.begin(), result.end(),
//                     make_optional(inserter(WorstPath,WorstPath.begin())) );

        // PRODUCE A LaTeX / TiKz DOCUMENT AND DISPLAY

        GraphPrinter tikz("scratch-graph");
        tikz.autoscale(points.begin(), points.end());
        tikz.drawEdges(result.begin(), result.end(), points, tikz.activeEdgeOptions);
        tikz.drawVerticesWithInfo(points.begin(), points.end(), tikz.activeVertexOptions);

        LatexPrinter latex("scratch-latex");
        latex.addToDocument(tikz);
        latex.display();


        cout << degree(result.begin(), result.end()) << endl;
        cout<<",";


        double stretchFactor;
        {
            Timer tim;
            stretchFactor = StretchFactorDijkstraReduction(points.begin(), points.end(), result.begin(),
                                                                  result.end());
        }
        cout<< stretchFactor;
        cout<<",";

        double stretchFactor2;
        {
            Timer tim;
            stretchFactor2 = StretchFactorUsingHeuristic(points.begin(), points.end(), result.begin(),
                                                                result.end());
        }
        cout<< stretchFactor2;
        cout<<",";

        double stretchFactor3;
        {
            Timer tim;
            stretchFactor3 = StretchFactorUsingHeuristic2(points.begin(), points.end(), result.begin(),
                                                                 result.end());
        }
        cout<< stretchFactor3;
        cout<<",";





//    graph.print("test");
        cout << "\n";
    }

    void scratch(size_t n) {
        double width = 10;

        vector<Point> points;
        points.reserve(n);
        generateRandomPoints(n, width / 2, back_inserter(points));

        scratch(points);
    }
    void scratch(string filename) {
        vector<Point> points;
        readPointsFromFile( back_inserter( points ), filename );
        scratch(points);
    }

} // planespanners

#endif //PLANESPANNERS_SCRATCH_H
