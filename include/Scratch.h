#ifndef SPANNERS_SCRATCH_H
#define SPANNERS_SCRATCH_H



#include <iostream>
#include <iterator>
#include <list>
#include <optional>
#include <string>
#include <utility>
#include <vector>



#include "algorithms/BoundedDegreePlaneSpanners.h"

#include "printers/GraphPrinter.h"
#include "printers/LatexPrinter.h"
#include "printers/PgfplotPrinter.h"
#include "printers/TablePrinter.h"

#include "tools/Metrics.h"
#include "tools/PointGenerators.h"
#include "tools/Results.h"
#include "tools/Utilities.h"

namespace spanners {

    const string SCRATCH_DIRECTORY = "../scratch/";


    void scratch(const vector<Point>& points) {

        using namespace std;
        GraphPrinter tikz("/tmp/", "deg3");
        tikz.autoscale(points.begin(), points.end());


//        for(auto p : points)
//            if(p.x()>4.35)
//                spanners::bcc2012::POINT_COLLECTOR.emplace(p.x()-4.75,p.y());

        writePointsToFile(points.begin(), points.end(),"galaxy");

        cout << points.size();
        cout << "\n";

        list<pair<size_t, size_t> > result;

        { // RUN THE ALGORITHM(S) /////////////////////////////////////
//            Timer tim;
//            LW2004( points.begin(), points.end(), back_inserter(result));
//            BSX2009( points.begin(), points.end(), back_inserter(result), 2*PI/3);
//            BGS2005( points.begin(), points.end(), back_inserter(result));
//            KPX2010( points.begin(), points.end(), back_inserter(result), 18);
//            BCC2012<6>( points.begin(), points.end(), back_inserter(result));
//            BCC2012<7>( points.begin(), points.end(), back_inserter(result));
//            BHS2018(points.begin(), points.end(), back_inserter(result));
//            KPT2017(points.begin(), points.end(), back_inserter(result));
//            BKPX2015(points.begin(), points.end(), back_inserter(result));
//            BGHP2010(points.begin(), points.end(), back_inserter(result));
//            KX2012(points.begin(), points.end(), back_inserter(result));
            DEG3(points.begin(),points.end(),back_inserter(result));
//            delaunay_testing( points.begin(), points.end(), back_inserter(result));
        }

        cout << degree(result.begin(), result.end()) << endl;
//        for( auto edge : result ) {
//            cout<<edge.first<<"("<<points.at(edge.first)<<") -- "
//                <<edge.second<<"("<<points.at(edge.second)<<")\n";
//        }






        cout<< "EXP: time=";
        number_t t_exp = INF;
        {
            Timer tim;
            t_exp = StretchFactorExpDijk(points.begin(),points.end(),result.begin(),result.end());
        }
        cout<<"  t="<<t_exp<<";\n";

        GraphPrinter::OptionsList edgeOptions = { // active edge options
                {"color",      tikz.activeEdgeColor},
                {"line width", to_string(tikz.inactiveEdgeWidth/2)}
        };

        tikz.drawEdges(result.begin(),result.end(),points,edgeOptions);



        tikz.drawVertices(points.begin(), points.end(), tikz.activeVertexOptions);




        list<Edge> WorstPath;
        SFWorstPath( points.begin(),points.end(),result.begin(), result.end(),
                     make_optional(inserter(WorstPath,WorstPath.begin())) );

        // PRODUCE A LaTeX / TiKz DOCUMENT AND DISPLAY

        tikz.addLatexComment("worst path edges");
//        tikz.drawEdges(WorstPath.begin(),WorstPath.end(), points, tikz.highlightEdgeOptions);


        tikz.display();


        cout << "\n";
    }

    void scratchN(size_t n) {
        double width = 10;

        vector<Point> points;

        PointGenerator_2 generator;
        generator.inGalaxy(n,5,points);

        scratch(points);
    }
    void scratchFile(string filename) {
        vector<Point> points;
        readPointsFromFile( back_inserter( points ), filename );
        scratch(points);
    }

} // spanners

#endif //SPANNERS_SCRATCH_H
