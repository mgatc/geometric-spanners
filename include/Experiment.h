#ifndef GEOMETRIC_SPANNERS_EXPERIMENT_H
#define GEOMETRIC_SPANNERS_EXPERIMENT_H

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
#include "BHS2017.h"
#include "BKPX2015.h"
#include "BSX2009.h"
#include "KPT2017.h"
#include "KPX2010.h"
#include "KX2012.h"
#include "LW2004.h"

namespace unf_spanners {

    using std::to_string;

    const bool PRINT_GEOMETRY = false;
    const bool PRINT_PGFPLOTS = true;
    const bool PRINT_TABLES = true;

    const bool USE_EXACT_STRETCH_FACTOR = true;

    vector<BoundedDegreeSpannerResultSet> RESULTS(DistributionTypeLast);

    LatexPrinter latex("output-main");
    GraphPrinter graph("output-vis");
    PgfplotsPrinter pgfplots("output-plots");

    GraphPrinter::OptionsList activeEdgeOptions = { // active edge options
            {"color",      graph.activeEdgeColor},
            {"line width", to_string(graph.activeEdgeWidth)}
    },
    highlightEdgeOptions = { // active edge options
            {"densely dashed", ""},
            {"color",          graph.worstPathEdgeColor},
            {"line width",     to_string(graph.activeEdgeWidth)}
    },
    highlightVertexOptions = {
            {"diamond",    ""},
            {"vertex",     (to_string(graph.vertexRadius * 1.61))}, // vertex width
            {"color",      (graph.worstPathEdgeColor)}, // text color
            {"fill",       (graph.worstPathEdgeColor)}, // vertex color
            {"line width", (to_string(0))} // vertex border (same color as text)
    },
    activeVertexOptions = {
            {"circle",     ""},
            {"vertex",     (to_string(graph.vertexRadius))}, // vertex width
            {"color",      (graph.backgroundColor)}, // text color
            {"fill",       (graph.activeVertexColor)}, // vertex color
            {"line width", (to_string(0))} // vertex border (same color as text)
    },
    borderOptions = {
            {"border",     (to_string(graph.vertexRadius))}, // choose shape of vertex
            {"color",      graph.activeEdgeColor}, // additional border color
            {"line width", to_string(graph.inactiveEdgeWidth)}, // additional border width
    };
    template<class Container>
    void PlaneSpannerExperiment(const Container &points,
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
            case Algorithm::Bhs2017:
                BHS2017(pointsBegin, pointsEnd, back_inserter(spanner));
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
                //assert(false);
        }

        number_t runtime = tim.stop();
        size_t n = points.size();
        size_t deg = degree( spanner.begin(), spanner.end() );
        number_t degAvg = degreeAvg( spanner.begin(), spanner.end() );
        number_t lightness = getLightness( points.begin(), points.end(), spanner.begin(), spanner.end() );
        number_t t = USE_EXACT_STRETCH_FACTOR ?
            StretchFactorDijkstraReduction( points.begin(), points.end(), spanner.begin(), spanner.end() )
            : StretchFactorUsingHeuristic( points.begin(), points.end(), spanner.begin(), spanner.end() );

        list<Edge> WorstPath;
        SFWorstPath( points.begin(), points.end(), spanner.begin(), spanner.end(),
                     make_optional(inserter(WorstPath,WorstPath.begin())) );

        BoundedDegreeSpannerResult result(algorithm, n, runtime, deg, degAvg, lightness, t );

        RESULTS.at(dist).registerResult(result);
        cout << result << endl;
        if(PRINT_GEOMETRY){
            string outputname = string("RPIS-")
                                + to_string(n)
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
            printer.drawEdges( spanner.begin(), spanner.end(), points, activeEdgeOptions);

            printer.addLatexComment("worst path edges");
            printer.drawEdges(WorstPath.begin(),WorstPath.end(), points, highlightEdgeOptions);

            printer.addLatexComment("vertices");
            printer.drawVertices( points.begin(), points.end(), activeVertexOptions);

            printer.addLatexComment("worst pair vertices");
            vector<Point> worstPair = { points[WorstPath.front().first], points[WorstPath.back().second] };
            printer.drawVertices( worstPair.begin(), worstPair.end(), highlightVertexOptions);

            latex.addToDocumentAsFigure(printer);
        }
    }

    bool singleRun( const size_t n, DistributionType dist, const double width ){ //}, string resultFilename, optional<string> filename, bool forcePrint, bool printLog )

        // SET POINT SET
        vector<Point> points;

        switch(dist) {
        case UniformInsideSquare:
            generatePointsInsideASquare(n,width,points);
            break;
        case UniformInsideDisc:
            generatePointsInsideADisc(n,width/2,points);
            break;
        case NormalInsideSquare:
            generatePointsInsideASquareNormal(n,1,points);
            break;
        case NormalClustersInsideSquare:
            generatePointsInsideASquareNormal(ceil(sqrt(n)),floor(sqrt(n)),points);
            break;
        case ContiguousGrid:
            generateContiguousPointsOnAGrid(n, width, points);
            break;
        case UniformRandomGrid:
            generateRandomPointsOnAGrid(n,width, points);
            break;
//        case UniformInsideAnnulus:
//            generateRandomInsideAnnulus(n, width, width*0.63, points);
//            break;
//        case Real:
//            generatePointsFromFile(n, points);
//            break;

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

            printer.drawVertices(points.begin(), points.end(), activeVertexOptions);

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

    bool experiment( size_t trials, size_t n_start, size_t n_end, size_t increment ) {
        const double width = 10;

        size_t invalid = 0;
        size_t trialNum = 1;

        for( int dist=DistributionTypeFirst; dist!=DistributionTypeLast; ++dist ) {
            auto& result = RESULTS.at(dist);
            for (size_t trial = 0; trial <= trials; ++trial) {
                for (size_t n = n_start; n <= n_end; n += increment) {
                    trialNum++;
                    singleRun(n, static_cast<DistributionType>(dist), width);//, "output", nullopt, false, false ) )
                }
            }
            result.computeStatistics();

            if (PRINT_PGFPLOTS) {
                // Create plot names
                vector<string> plotNames;
                transform(PGFPLOT_NAMES.begin(),
                          PGFPLOT_NAMES.end(),
                          back_inserter(plotNames),
                          [&dist](const string& str) {
                              string plotName = str + " (" + DISTRIBUTION_NAMES.at(dist) + ")";
                              //cout<<plotName;
                              return plotName;
                          });
                plotResults(result, &latex, plotNames);
            }
            string tableFilename = string("output-table")
                    + DISTRIBUTION_NAMES.at(dist);
            TablePrinter table(tableFilename);
            if (PRINT_TABLES) {
                table.ignoreIV(0); // ignore runtime
                table.addColumn(DEGREE_BOUND_SYMBOL, DEGREE_BOUND_PER_ALGORITHM);
                table.addColumn(STRETCH_FACTOR_BOUND_SYMBOL, STRETCH_FACTOR_BOUND_PER_ALGORITHM);
                table.tabulateResults(result);
                latex.addToDocumentAsFigure(table);
            }

        }

        if (PRINT_GEOMETRY || PRINT_PGFPLOTS || PRINT_TABLES)
            latex.display();

        //cout<<"\nTesting Complete. "<< invalid << " of "<<(trials*(n_end-n_start))<<" invalid results.\n\n";
        return invalid == 0;
    }





} // unf_spanners

#endif //GEOMETRIC_SPANNERS_EXPERIMENT_H
