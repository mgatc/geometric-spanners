#ifndef SPANNERS_EXPERIMENT_H
#define SPANNERS_EXPERIMENT_H

#include <iostream>
#include <list>
#include <optional>
#include <utility>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "algorithms/Degree3.h"

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


namespace spanners {

    const string INPUT_DATA_DIRECTORY = "../input/";
    const string OUTPUT_DATA_DIRECTORY = "../output/";

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

    LatexPrinter latex(OUTPUT_DATA_DIRECTORY, "exp-main");
    GraphPrinter graph(OUTPUT_DATA_DIRECTORY, "exp-vis");
    PgfplotPrinter pgfplots(OUTPUT_DATA_DIRECTORY, "exp-plots");

    map<string,BoundedDegreeSpannerResultSet> RESULTS;

    void SingleTrial (const vector<Point>& points, const string dist );
    void PlaneSpannerTest( const vector<Point>&,const DistributionType,const Algorithm);
    void SyntheticTrial(const size_t n, DistributionType dist, const double width);

    void ExperimentFromConfigurationFile(size_t numRuns, string configFilename) {
        boost::property_tree::ptree config;

        using std::string, std::vector;
        namespace pt = boost::property_tree;
        pt::read_xml(configFilename,config);

        auto& result = RESULTS.emplace(configFilename, BoundedDegreeSpannerResultSet()).first->second;

        map<index_t, string> PointsetNames;

        for( auto pointset : config.get_child("pointsets") ) {

            string filename = pointset.second.get_child("filename").data(),
                   fullname = INPUT_DATA_DIRECTORY + filename,
                   filenameNoExtension = filename;
            boost::erase_all(filenameNoExtension, ".xy");

            std::ifstream in(fullname);

            if (in.is_open()) {
                vector<Point> P;
                number_t x,y;
                while ( in >> x >> y ) {
                    P.emplace_back(x,y);
                }
                in.close();

                const index_t n = P.size();
                string pointsetName = pointset.second.get_child("nicename").data();
                PointsetNames.emplace(n,pointsetName);

                cout<< "!! Starting  "<< pointsetName << " trials !!\n"
                    << "Added "<< n <<" points from file\n"

                    << "NOTE: one extra trial is performed because trial 0 will be thrown out!"<<endl<<endl;

                for (size_t trial = 0; trial <= numRuns; ++trial) {
                    SingleTrial(P, configFilename);
                }

                cout<< "!! Ending  "<< pointsetName << " trials !!\n"
                    << "-------------------------------------------"<<endl;

            } else {
                cout<<"Error opening file!\n";
            }

        }

        result.computeStatistics();

        string experimentName = configFilename;
        boost::erase_all(experimentName, ".xml");
        boost::erase_all(experimentName, ".");
        boost::erase_all(experimentName, "/");

        if (PRINT_PGFPLOTS) {
            plotResults(experimentName, result, &latex);
        }
        if(PRINT_IV_TABLES){
            tabulateIVsFromConfigExperiment(experimentName, result, PointsetNames, &latex);
        }
        if (PRINT_SUMMARY_TABLES) {
            tabulateSummaryResults(experimentName, result, &latex);
        }

        if (PRINT_GEOMETRY || PRINT_PGFPLOTS || PRINT_SUMMARY_TABLES || PRINT_IV_TABLES)
            latex.display();

    }


    void SyntheticExperiment(size_t numRuns, size_t n_start, size_t n_end, size_t increment ) {

        // EXPERIMENT PARAMETER OVERRIDE ----------------------------------------//
        numRuns = numRuns, n_start = n_start, n_end = n_end, increment = increment;
        // ----------------------------------------------------------------------//

        const number_t width = 1000;

        for( int dist=DistributionTypeFirst; dist!=DistributionTypeLast; ++dist ) {
            auto distributionType = static_cast<DistributionType>(dist);
            string distName = DISTRIBUTION_NAMES.at(dist);

            cout<< "!! Starting  "<< distName << "distribution trials !!\n"
                << "NOTE: one extra trial is performed because trial 0 will be thrown out!"<<endl<<endl;

            auto& result = (RESULTS.emplace(distName, BoundedDegreeSpannerResultSet()).first)->second;
            for (size_t trial = 0; trial <= numRuns; ++trial) {
                cout<< "Starting trial "<< trial << "..."<<endl<<endl;
                for (size_t n = n_start; n <= n_end; n += increment) {
                    SyntheticTrial(n, distributionType, width);//, "output", nullopt, false, false ) )
                }
                cout<<"\n\n";
            }
            result.computeStatistics();

            if (PRINT_PGFPLOTS) {
                plotResults(distName, result, &latex);
            }
            if(PRINT_IV_TABLES){
                tabulateIVs(distName, result, &latex);
            }
            if (PRINT_SUMMARY_TABLES) {
                tabulateSummaryResults(distName, result, &latex);
            }
            cout<< "!! Ending  "<< distName << "distribution trials !!\n"
                << "-------------------------------------------"<<endl;


        }

        if (PRINT_GEOMETRY || PRINT_PGFPLOTS || PRINT_SUMMARY_TABLES || PRINT_IV_TABLES)
            latex.display();

//        cout<<"Astar Wrong Count="<<WRONG_COUNT_2<<" Average Amount="<<(WRONG_COUNT_2 > 0 ? WRONG_AMOUNT_2/WRONG_COUNT_2 : 0)<<endl;
//        cout<<"      Wrong Rate="<<(100*((double)WRONG_COUNT_2/EXP_COUNT))<<"%"<<endl;
//        cout<<"Dijkstra Wrong Count="<<WRONG_COUNT_3<<" Average Amount="<<(WRONG_COUNT_3 > 0 ? WRONG_AMOUNT_3/WRONG_COUNT_3 : 0)<<endl;
//        cout<<"      Wrong Rate="<<(100*((double)WRONG_COUNT_3/EXP_COUNT))<<"%"<<endl;
//        cout<<endl<<"Total exp="<<EXP_COUNT<<endl;

        //cout<<"\nTesting Complete. "<< invalid << " of "<<(numRuns*(n_end-n_start))<<" invalid results.\n\n";
    }




    void PlaneSpannerTest(const vector<Point> &points,
                          const string dist,
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
            case Algorithm::Bkpx2015:
                BKPX2015(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Kpt2017:
                KPT2017(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
//            case Algorithm::Degree3:
//                DEG3(pointsBegin, pointsEnd, back_inserter(spanner));
//                break;
            case Algorithm::AlgorithmLast:
                assert(false);
        }

        number_t runtime = tim.stop();
        BoundedDegreeSpannerResult result(algorithm, runtime, points.begin(), points.end(), spanner.begin(), spanner.end());
        cout << result;

        ++EXP_COUNT;


        RESULTS.at(dist).registerResult(result);

        if(PRINT_GEOMETRY){
            list<Edge> WorstPath;
            SFWorstPath( points.begin(), points.end(), spanner.begin(), spanner.end(),
                         make_optional(inserter(WorstPath,WorstPath.begin())) );

            string outputname = "visual-"
                                + dist
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

    void SingleTrial (const vector<Point>& points, const string dist ){
        const size_t n = points.size();

        cout<< "Starting trial..."<<endl<<endl;
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
          alg!=Algorithm::AlgorithmLast; ++alg ) {
            PlaneSpannerTest(points, dist, static_cast<Algorithm>(alg));
        }

        cout<< "Finished trial...\n"
            << "-----------------"<<endl;
    }


    void SyntheticTrial(const size_t n, DistributionType dist, const double width ){ //}, string resultFilename, optional<string> filename, bool forcePrint, bool printLog )


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

        SingleTrial(points, DISTRIBUTION_NAMES.at(dist));

    }







} // spanners

#endif //SPANNERS_EXPERIMENT_H
