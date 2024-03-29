#ifndef SPANNERS_SCRATCH_H
#define SPANNERS_SCRATCH_H



#include <iostream>
#include <iterator>
#include <list>
#include <optional>
#include <string>
#include <utility>
#include <vector>



#include "libspanner/BoundedDegreePlaneSpanners.h"
#include "libspanner/delaunay/DelaunayL2.h"
#include "libspanner/delaunay/DelaunayLinf.h"
#include "libspanner/delaunay/DelaunayTD.h"

#include "libspanner/measure/degree.h"
#include "libspanner/measure/stretchfactor.h"
#include "libspanner/measure/timer.h"
#include "libspanner/measure/weight.h"

#include "libspanner/points/generators.h"
#include "libspanner/types.h"
#include "libspanner/utilities.h"



#include "cpptex/cpptex.h"



#include "Results.h"

namespace bdps_experiment {

    void scratch(const spanner::bdps::input_t& points) {

        using namespace spanner;



//        for(auto p : points)
//            if(p.x()>4.35)
//                bdps_experiment::bcc2012::POINT_COLLECTOR.emplace(p.x()-4.75,p.y());

        writePointsToFile(points.begin(), points.end(),"scratch");

        std::cout << points.size();
        std::cout << "\n";

        bdps::output_t result;

        { // RUN THE ALGORITHM(S) /////////////////////////////////////
            Timer tim;
//            LW2004( points, result);
//            BSX2009( points, result, 2*PI/3);
//            BGS2005( points, result);
//            KPX2010( points, result, 18);
//            BCC2012<6>( points, result);
            BCC2012<7>( points, result);
//            BHS2018(points, result);
//            KPT2017(points, result);
//            BKPX2015(points, result);
//            BGHP2010(points, result);
//            KX2012(points, result);
//            DEG3(points,result);
//            DelaunayL2Spanner(points,result);
//            DelaunayLinfSpanner(points,result);
//            DelaunayTDSpanner(points,result);
        }

        std::cout << degree(result.begin(), result.end()) << std::endl<<std::endl;

        std::sort(result.begin(),result.end(),[](const auto& a, const auto& b) {
            return a.first < b.first || (a.first==b.first && a.second<b.second);
        });

        for( auto edge : result ) {
            std::cout<<edge.first+1<<" "
                     <<edge.second+1<<"\n";
        }
        std::cout<<std::endl;






        std::cout<< "EXP: time=";
        number_t t_exp = INF;
        {
            Timer tim;
            t_exp = spanner::StretchFactorExpDijk(points.begin(),points.end(),result.begin(),result.end());
        }
        std::cout<<"  t="<<t_exp<<";\n";



        cpptex::GraphPrinter tikz("/tmp/scratch", points.begin(), points.end(), 10.0);

        tikz.drawEdges(result.begin(),result.end(),points,tikz.activeEdgeOptions);
        tikz.drawVertices(points.begin(), points.end(), tikz.activeVertexOptions);




        std::list<Edge> WorstPath;
        spanner::SFWorstPath( points.begin(),points.end(),result.begin(), result.end(),
                     make_optional(inserter(WorstPath,WorstPath.begin())) );

        // PRODUCE A LaTeX / TiKz DOCUMENT AND DISPLAY

        tikz.addComment("worst path edges");
        tikz.drawEdges(WorstPath.begin(),WorstPath.end(), points, tikz.highlightEdgeOptions);


        tikz.display();


        std::cout << "\n";
    }

    void scratchN(size_t n) {
        double width = 10;

        spanner::bdps::input_t points;

        spanner::PointGenerator_2 generator;
        generator.inGalaxy(n,5,points);

        scratch(points);
    }
    void scratchFile(std::string filename) {
        spanner::bdps::input_t points;
        spanner::readPointsFromFile( back_inserter( points ), filename );
        scratch(points);
    }

} // bdps_experiment

#endif //SPANNERS_SCRATCH_H
