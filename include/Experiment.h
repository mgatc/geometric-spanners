#ifndef PLANESPANNERS_EXPERIMENT_H
#define PLANESPANNERS_EXPERIMENT_H

#include <iostream>
#include <list>
#include <optional>
#include <utility>

#include "LatexPrinter.h"
#include "LineGraphPrinter.h"
#include "metrics.h"
#include "names.h"
#include "PointGenerators.h"
#include "Result.h"
#include "TablePrinter.h"
#include "utilities.h"


#include "BCC2012.h"
#include "BGHP2010.h"
#include "BGS2005.h"
#include "BHS2018.h"
#include "BKPX2015.h"
#include "BSX2009.h"
#include "KPT2017.h"
#include "KPX2010.h"
#include "KX2012.h"
#include "LW2004.h"

namespace planespanners {

    using std::to_string;

    const bool PRINT_GEOMETRY = false;
    const bool PRINT_PGFPLOTS = true;
    const bool PRINT_IV_TABLES = true;
    const bool PRINT_SUMMARY_TABLES = true;

    size_t WRONG_COUNT_2 = 0;
    number_t WRONG_AMOUNT_2 = 0.0;
    size_t WRONG_COUNT_3 = 0;
    number_t WRONG_AMOUNT_3 = 0.0;
    size_t EXP_COUNT = 0;

    LatexPrinter latex("exp-main");
    GraphPrinter graph("exp-vis");
    PgfplotsPrinter pgfplots("exp-plots");

    vector<BoundedDegreeSpannerResultSet> RESULTS(DistributionTypeLast);

    void PlaneSpannerExperiment( const vector<Point>&,const DistributionType,const Algorithm);
    bool singleRun( const size_t,const DistributionType,const double);



    bool experiment(size_t numRuns, size_t n_start, size_t n_end, size_t increment ) {

        // EXPERIMENT PARAMETER OVERRIDE ----------------------------------------//
        numRuns = numRuns, n_start = n_start, n_end = n_end, increment = increment;
        // ----------------------------------------------------------------------//

        const number_t width = 1000;

        for( int dist=DistributionTypeFirst; dist!=DistributionTypeLast; ++dist ) {
            auto distibutionType = static_cast<DistributionType>(dist);
            auto& result = RESULTS.at(dist);
            for (size_t trial = 0; trial <= numRuns; ++trial) {
                for (size_t n = n_start; n <= n_end; n += increment) {
                    singleRun(n, distibutionType, width);//, "output", nullopt, false, false ) )
                }
            }
            result.computeStatistics();

            if (PRINT_PGFPLOTS) {
                plotResults(distibutionType, result, &latex);
            }
            if(PRINT_IV_TABLES){
                tabulateIVs(distibutionType, result, &latex);
            }
            if (PRINT_SUMMARY_TABLES) {
                tabulateSummaryResults(distibutionType, result, &latex);
            }

        }

        if (PRINT_GEOMETRY || PRINT_PGFPLOTS || PRINT_SUMMARY_TABLES || PRINT_IV_TABLES)
            latex.display();

//        cout<<"Astar Wrong Count="<<WRONG_COUNT_2<<" Average Amount="<<(WRONG_COUNT_2 > 0 ? WRONG_AMOUNT_2/WRONG_COUNT_2 : 0)<<endl;
//        cout<<"      Wrong Rate="<<(100*((double)WRONG_COUNT_2/EXP_COUNT))<<"%"<<endl;
//        cout<<"Dijkstra Wrong Count="<<WRONG_COUNT_3<<" Average Amount="<<(WRONG_COUNT_3 > 0 ? WRONG_AMOUNT_3/WRONG_COUNT_3 : 0)<<endl;
//        cout<<"      Wrong Rate="<<(100*((double)WRONG_COUNT_3/EXP_COUNT))<<"%"<<endl;
//        cout<<endl<<"Total exp="<<EXP_COUNT<<endl;

        //cout<<"\nTesting Complete. "<< invalid << " of "<<(numRuns*(n_end-n_start))<<" invalid results.\n\n";
        return true;
    }




    void PlaneSpannerExperiment(const vector<Point> &points,
                                const DistributionType dist,
                                const Algorithm algorithm ) {
        using namespace std;

        list<pair<size_t, size_t> > spanner;

        auto pointsBegin = points.begin(),
                pointsEnd = points.end();

        Timer tim;

        switch (algorithm) {
            case Algorithm::Bgs2005:
                BGS2005(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Lw2004:
                LW2004(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bsx2009:
                BSX2009(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Kpx2010:
                KPX2010(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Kx2012:
                KX2012(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bcc2012_7:
                BCC2012<7>(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bcc2012_6:
                BCC2012<6>(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bhs2018:
                BHS2018(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bghp2010:
                BGHP2010(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Kpt2017:
                KPT2017(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bkpx2015:
                BKPX2015(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::AlgorithmLast:
                assert(false);
        }

        number_t runtime = tim.stop();
        BoundedDegreeSpannerResult result(algorithm, runtime, points.begin(), points.end(), spanner.begin(), spanner.end());
        cout << result;

        ++EXP_COUNT;

//        double stretchFactor;
//        {
//            Timer tim;
//            stretchFactor = StretchFactorDijkstraReduction(points.begin(), points.end(), spanner.begin(),
//                                                           spanner.end());
//        }
//        cout<< stretchFactor;
//        cout<<",";
//
//        double stretchFactor2;
//        {
//            Timer tim;
//            stretchFactor2 = StretchFactorUsingHeuristic(points.begin(), points.end(), spanner.begin(),
//                                                         spanner.end());
//        }
//        cout<< stretchFactor2;
//        if(abs(stretchFactor2-stretchFactor) > EPSILON) {
//            cout<<"(WRONG)";
//            WRONG_COUNT_2++;
//            WRONG_AMOUNT_2 += abs(stretchFactor2-stretchFactor);
//        }
//        cout<<",";
//
//        double stretchFactor3;
////        list<Edge> WorstPathHeuristic;
////        SFWorstPath2( points.begin(), points.end(), spanner.begin(), spanner.end(),
////                     make_optional(inserter(WorstPathHeuristic,WorstPathHeuristic.begin())) );
//        {
//            Timer tim;
//            stretchFactor3 = StretchFactorUsingHeuristic2(points.begin(), points.end(), spanner.begin(),
//                                                          spanner.end());
//        }
//        cout<< stretchFactor3;
//        if(abs(stretchFactor3-stretchFactor) > EPSILON) {
//            cout<<"(WRONG)";
//            WRONG_COUNT_3++;
//            WRONG_AMOUNT_3 += abs(stretchFactor3-stretchFactor);
//            //print and return
//            GraphPrinter tikz("scratch-graph");
//            tikz.autoscale(points.begin(), points.end());
//            tikz.drawEdges(spanner.begin(), spanner.end(), points, tikz.inactiveEdgeOptions);
//            tikz.drawVerticesWithInfo(points.begin(), points.end(), tikz.activeVertexOptions);
//
//            list<Edge> WorstPath;
//            SFWorstPath( points.begin(), points.end(), spanner.begin(), spanner.end(),
//                         make_optional(inserter(WorstPath,WorstPath.begin())) );
//            tikz.drawEdges(WorstPath.begin(),WorstPath.end(), points, tikz.activeEdgeOptions);
//
//
//            LatexPrinter latex("scratch-latex");
//            latex.addToDocument(tikz);
//            latex.display();
//            assert(!"WRONG STRETCH FACTOR");
//        }
//        cout<<endl<<endl;


        RESULTS.at(dist).registerResult(result);

        if(PRINT_GEOMETRY){
            list<Edge> WorstPath;
            SFWorstPath( points.begin(), points.end(), spanner.begin(), spanner.end(),
                         make_optional(inserter(WorstPath,WorstPath.begin())) );

            string outputname = "visual-"
                                + DISTRIBUTION_NAMES.at(dist)
                                + to_string(points.size())
                                + "-"
                                + ALGORITHM_NAMES.at(algorithm);
            GraphPrinter printer(outputname);

            double documentSizeInCm = 10;
            printer.autoscale( points.begin(), points.end(), documentSizeInCm);

            result.setCaption(printer);

            // remove worst path edges from spanner edges
            spanner.erase( remove_if( spanner.begin(), spanner.end(), [&]( const Edge& edge ){
                auto spannerEdge = normalize_pair(edge);
                return find_if(WorstPath.begin(), WorstPath.end(), [&] ( const Edge& wpEdge ) {
                    return spannerEdge == normalize_pair(wpEdge);
                }) != WorstPath.end();
            }), spanner.end());

            printer.addLatexComment("spanner edges (except for worst path)");
            printer.drawEdges( spanner.begin(), spanner.end(), points, printer.activeEdgeOptions);

            printer.addLatexComment("worst path edges");
            printer.drawEdges(WorstPath.begin(),WorstPath.end(), points, printer.highlightEdgeOptions);

            printer.addLatexComment("vertices");
            printer.drawVertices( points.begin(), points.end(), printer.activeVertexOptions);

            printer.addLatexComment("worst pair vertices");
            vector<Point> worstPair = { points[WorstPath.front().first], points[WorstPath.back().second] };
            printer.drawVertices( worstPair.begin(), worstPair.end(), printer.highlightVertexOptions);

            latex.addToDocumentAsFigure(printer);
        }
    }

    bool singleRun( const size_t n, DistributionType dist, const double width ){ //}, string resultFilename, optional<string> filename, bool forcePrint, bool printLog )

        // SET POINTS
        vector<Point> points;

        switch(dist) {
        case UniformInsideSquare:
            generatePointsInsideASquare(n,width,points);
            break;
        case UniformInsideDisc:
            generatePointsInsideADisc(n,width,points);
            break;
        case NormalInsideSquare:
            generatePointsInsideASquareNormal(n,1,points);
            break;
        case NormalClustersInsideSquare:
            generatePointsInsideASquareNormal(n/pow(n,(1/3)),
                                              pow(n,(1/3)),points);
            break;
        case ContiguousGrid:
            generateContiguousPointsOnAGrid(n, points);
            break;
        case UniformRandomGrid:
            generateRandomPointsOnAGrid(n, points);
            break;
        case UniformInsideAnnulus:
            generateRandomInsideAnnulus(n, width, width*0.8, points);
            break;
//        case Real:
//            generatePointsFromFile(n, points);
//            break;
        case DistributionTypeLast:
        default:
            assert(!"Invalid distribution type!");
        }

        if(PRINT_GEOMETRY) {
            string outputname = string("RPIS-") + to_string(n);
            GraphPrinter printer(outputname);
            double documentSizeInCm = 10;
            printer.autoscale(points.begin(), points.end(), documentSizeInCm);
            GraphPrinter::clearCaptionFile();

            // draw point set
            string caption = string("$N=")
                             + to_string(n)
                             + "$";
            printer.setCaption(caption);

            printer.drawVertices(points.begin(), points.end(), printer.activeVertexOptions);

            latex.addToDocumentAsFigure(printer);
        }
        for(int alg=Algorithm::AlgorithmFirst;
             alg!=Algorithm::AlgorithmLast; ++alg )
        {
            PlaneSpannerExperiment( points, dist, static_cast<Algorithm>(alg) );
        }

        cout<<endl;

        return true;
    }







} // planespanners

#endif //PLANESPANNERS_EXPERIMENT_H
