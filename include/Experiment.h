#ifndef SPANNERS_EXPERIMENT_H
#define SPANNERS_EXPERIMENT_H

#include <iostream>
#include <list>
#include <optional>
#include <random>
#include <utility>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "algorithms/BoundedDegreePlaneSpanners.h"

#include "printers/LatexPrinter.h"
#include "printers/PgfplotPrinter.h"
#include "printers/TablePrinter.h"

#include "tools/Metrics.h"
#include "tools/PointGenerators.h"
#include "tools/Results.h"
#include "tools/Utilities.h"


namespace spanners {

    const bool USE_EXACT_STRETCH_FACTOR = false;

    const string OUTPUT_DATA_DIRECTORY = "./";
    const string OUTPUT_EXTENSION = ".csv";
    const string DELIMITER = ",";

    using std::to_string;

    const bool PRINT_GEOMETRY = false;
    const bool PRINT_PGFPLOTS = true;
    const bool PRINT_IV_TABLES = false;
    const bool PRINT_SUMMARY_TABLES = false;

    const number_t INPUT_WIDTH = 1000;

    size_t WRONG_COUNT_DIJKSTRA = 0;
    number_t WRONG_AMOUNT_DIJKSTRA = 0.0;
    size_t EXP_COUNT = 0;

    void BoundedDegreePlaneSpannerTest(const BoundedDegreePlaneSpannerAlgorithm algorithm,
                                       const DistributionType dist,
                                       const vector<Point> &points,
                                       ofstream& expOut,
                                       bool lite = false) {
        using namespace std;

        list<pair<size_t, size_t> > spanner;

        auto pointsBegin = points.begin(),
                pointsEnd = points.end();

        Timer tim;

        switch (algorithm) {
            case BoundedDegreePlaneSpannerAlgorithm::Bgs2005:
                BGS2005(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case BoundedDegreePlaneSpannerAlgorithm::Lw2004:
                LW2004(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case BoundedDegreePlaneSpannerAlgorithm::Bsx2009:
                BSX2009(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case BoundedDegreePlaneSpannerAlgorithm::Kpx2010:
                KPX2010(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case BoundedDegreePlaneSpannerAlgorithm::Kx2012:
                KX2012(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case BoundedDegreePlaneSpannerAlgorithm::Bcc2012_7:
                BCC2012<7>(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case BoundedDegreePlaneSpannerAlgorithm::Bcc2012_6:
                BCC2012<6>(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case BoundedDegreePlaneSpannerAlgorithm::Bhs2018:
                BHS2018(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case BoundedDegreePlaneSpannerAlgorithm::Bghp2010:
                BGHP2010(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case BoundedDegreePlaneSpannerAlgorithm::Bkpx2015:
                BKPX2015(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case BoundedDegreePlaneSpannerAlgorithm::Kpt2017:
                KPT2017(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
//            case BoundedDegreePlaneSpannerAlgorithm::Degree3:
//                DEG3(pointsBegin, pointsEnd, back_inserter(spanner));
//                break;
            case BoundedDegreePlaneSpannerAlgorithm::AlgorithmLast:
                assert(false);
        }

        number_t runtime = tim.stop();

        BoundedDegreeSpannerResult result(dist, algorithm, runtime,
                                          points.begin(), points.end(),
                                          spanner.begin(), spanner.end(),
                                          !USE_EXACT_STRETCH_FACTOR);
        cout << result;
        expOut << result;

        if(!result.verify()) {
            string filename = "breaks";
            filename += ALGORITHM_NAMES.at(algorithm);
            filename += ".xy";
            writePointsToFile(points.begin(),points.end(),filename);

            cout<< "\n"
                << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                << "!!!!!!!!!!!!!!!! ALGORITHM ERROR !!!!!!!!!!!!!!!!!!!!"
                << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";

            char userContinues = 'x';
            while( userContinues != 'y' && userContinues != 'n' ) {
                cout<< "\n Do you wish to continue the experiment? (y/n) ";
                cin>>userContinues;
            }
            if(userContinues=='n') assert(!"Experiment terminated by user.");
        }


        ++EXP_COUNT;
    }

    // An experiment featuring multiple graphs on the same point set
    void BoundedDegreePlaneSpannerAlgorithmLoop(const DistributionType dist,
                                                const vector<Point>& points,
                                                ofstream& expOut ){
        const size_t n = points.size();

        cout<< "Starting trial..."<<endl<<endl;

        for(int alg=BoundedDegreePlaneSpannerAlgorithm::AlgorithmFirst;
            alg != BoundedDegreePlaneSpannerAlgorithm::AlgorithmLast; ++alg ) {
            auto algorithm = static_cast<BoundedDegreePlaneSpannerAlgorithm>(alg);
            list<pair<size_t, size_t> > spanner;

            BoundedDegreePlaneSpannerTest(algorithm, dist, points, expOut);

            ++EXP_COUNT;
        }

        cout<< "Finished trial...\n"
            << "-----------------"<<endl;
    }

    // Generates a random point set from the given distribution and runs an experiment with multiple graphs on the point set
    void generateRandomPointSet(const DistributionType dist,
                                const size_t n,
                                const double width,
                                vector<Point>& points){
        RandomPointGenerator_2 pointGenerator;

        switch(dist) {
            case UniformInsideSquare:
                pointGenerator.generatePointsInsideASquare(n,width,points);
                break;
            case UniformInsideDisc:
                pointGenerator.generatePointsInsideADisc(n,width,points);
                break;
//        case UniformOnSquare:
//            pointGenerator.generatePointsOnASquare(n,width,points);
//            break;
            case UniformOnDisc:
                pointGenerator.generatePointsOnADisc(n,width,points);
                break;
            case NormalInsideSquare:
                pointGenerator.generatePointsInsideASquareNormal(n,1,points);
                break;
            case NormalClustersInsideSquare:
                pointGenerator.generatePointsInsideASquareNormal(n/10, 10, points);
                break;
            case ContiguousGrid:
                pointGenerator.generateContiguousPointsOnAGrid(n, points);
                break;
            case UniformRandomGrid:
                pointGenerator.generateRandomPointsOnAGrid(n, points);
                break;
            case UniformInsideAnnulus:
                pointGenerator.generateRandomInsideAnnulus(n, width, width*0.8, points);
                break;
            case Galaxy:
                pointGenerator.generatePointsInGalaxy(n, 5, points);
                break;
            case DistributionTypeLast:
            default:
                assert(!"Invalid distribution type!");
        }

    }

    // An experiment from n_start to n_end with a single distribution
    void SyntheticExperimentInputSizeLoop(DistributionType dist, size_t n_start, size_t n_end, size_t increment, ofstream& expOut) {
        for (size_t n = n_start; n <= n_end; n += increment) {
            // SET POINTS
            vector<Point> points;
            generateRandomPointSet(dist, n, INPUT_WIDTH, points);
            BoundedDegreePlaneSpannerAlgorithmLoop(dist,points,expOut);
        }
    }

    // An experiment of numRuns trials for a single distribution
    void SyntheticExperimentRepetitionLoop(DistributionType dist,
                                           size_t numRuns,
                                           size_t n_start,
                                           size_t n_end,
                                           size_t increment,
                                           ofstream& expOut ) {
        string distName = DISTRIBUTION_NAMES.at(dist);
        srand(0);
        for (size_t trial = 0; trial < numRuns; ++trial) {
            cout<< "Starting trial "<< trial << " of "<<distName<<"\n\n";
            SyntheticExperimentInputSizeLoop(dist, n_start, n_end, increment, expOut);
            cout<<"\n"<<endl;
        }
    }

    void SyntheticExperimentDistributionLoop(size_t numRuns,
                                             size_t n_start,
                                             size_t n_end,
                                             size_t increment,
                                             ofstream& expOut) {
        const number_t width = 10;

        for( int dist=DistributionTypeFirst; dist!=DistributionTypeLast; ++dist ) {
            auto distributionType = static_cast<DistributionType>(dist);
            string distName = DISTRIBUTION_NAMES.at(dist);

            cout<< "!! Starting  "<< distName << "distribution trials !!\n"<<endl<<endl;

            SyntheticExperimentRepetitionLoop(distributionType, numRuns, n_start, n_end, increment, expOut);

            cout<< "!! Ending  "<< distName << "distribution trials !!\n"
                << "-------------------------------------------"<<endl;
        }
    }

    // An experiment of numRuns trials for each distribution
    void SyntheticExperiment(size_t numRuns, size_t n_start, size_t n_end, size_t increment ) {

        // get unix timestamp to use as experiment file name
        auto time = std::time(nullptr);
        auto filename = to_string(time) + OUTPUT_EXTENSION;

        ofstream expOut;
        expOut.open(filename,ios_base::out);

        if(!expOut.is_open()) assert(!"ERROR OPENING OUTPUT FILE\n\n");

        expOut << "distribution" << DELIMITER
               << "n" << DELIMITER
               << "spannerAlg" << DELIMITER
               << "stretchFactorExact" << DELIMITER
               << "runtimeExact" << DELIMITER
               << "stretchFactorExperimental" << DELIMITER
               << "runtimeExperimental" << DELIMITER
               << "\n";

        SyntheticExperimentDistributionLoop(numRuns,n_start,n_end,increment,expOut);

        expOut.close();
    }



//    void ExperimentFromConfigurationFile(size_t numRuns, string configFilename) {
//        boost::property_tree::ptree config;
//
//        using std::string, std::vector;
//        namespace pt = boost::property_tree;
//        pt::read_xml(configFilename,config);
//
////        auto& result = RESULTS.emplace(configFilename, BoundedDegreeSpannerResultSet()).first->second;
//
//        map<index_t, string> PointsetNames;
//
//
//        for( auto pointset : config.get_child("pointsets") ) {
//
//            string filename = pointset.second.get_child("filename").data(),
//                   fullname = INPUT_DATA_DIRECTORY + filename,
//                   filenameNoExtension = filename;
//            boost::erase_all(filenameNoExtension, ".xy");
//
//            std::ifstream in(fullname);
//
//            if (in.is_open()) {
//                vector<Point> P;
//                number_t x,y;
//                while ( in >> x >> y ) {
//                    P.emplace_back(x,y);
//                }
//                in.close();
//
//                const index_t n = P.size();
//                string pointsetName = pointset.second.get_child("nicename").data();
//                PointsetNames.emplace(n,pointsetName);
//
//                cout<< "!! Starting  "<< pointsetName << " trials !!\n"
//                    << "Added "<< n <<" points from file\n"
//
//                    << "NOTE: one extra trial is performed because trial 0 will be thrown out!"<<endl<<endl;
//
//                for (size_t trial = 0; trial <= numRuns; ++trial) {
//                    SingleTrial(P, DistributionType::Real, trial!=1);
//                }
//
//                cout<< "!! Ending  "<< pointsetName << " trials !!\n"
//                    << "-------------------------------------------"<<endl;
//
//            } else {
//                cout<<"Error opening file!\n";
//            }
//
//        }
//
////        result.computeStatistics(true);
//
//        string experimentName = configFilename;
//        boost::erase_all(experimentName, ".xml");
//        boost::erase_all(experimentName, ".");
//        boost::erase_all(experimentName, "/");
//
//        if (PRINT_PGFPLOTS) {
//            plotResults(experimentName, result, &latex);
//        }
//        if(PRINT_IV_TABLES){
//            tabulateIVsFromConfigExperiment(experimentName, result, PointsetNames, &latex);
//        }
//        if (PRINT_SUMMARY_TABLES) {
//            tabulateSummaryResults(experimentName, result, &latex);
//        }
//
//        if (PRINT_GEOMETRY || PRINT_PGFPLOTS || PRINT_SUMMARY_TABLES || PRINT_IV_TABLES)
//            latex.display();
//
//    }


} // spanners

#endif //SPANNERS_EXPERIMENT_H
