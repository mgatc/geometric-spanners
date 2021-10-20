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

    const bool USE_EXPERIMENTAL_STRETCH_FACTOR = true;

    const string INPUT_DATA_DIRECTORY = "../input/";
    const string OUTPUT_DATA_DIRECTORY = "./";
    const string OUTPUT_EXTENSION = ".csv";
    const string DELIMITER = ",";

    using std::to_string;

    const number_t INPUT_WIDTH = 1000;

    size_t EXP_COUNT = 0;


    /**
     * Constructs on spanner using the given algorithm on the given point set,
     * then writes the result to the expOut file. Due to its
     */
    template <typename DistributionSubTypeEnum>
    void BoundedDegreePlaneSpannerTest(const BoundedDegreePlaneSpannerAlgorithm algorithm,
                                       const DistributionType distributionType,
                                       const DistributionSubTypeEnum distribution,
                                       const vector<Point> &points,
                                       ofstream& expOut,
                                       bool measureStretchFactor = true) {
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
//            case BoundedDegreePlaneSpannerAlgorithm::Bghp2010:
//                BGHP2010(pointsBegin, pointsEnd, back_inserter(spanner));
//                break;
            case BoundedDegreePlaneSpannerAlgorithm::Bkpx2015:
                BKPX2015(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
//            case BoundedDegreePlaneSpannerAlgorithm::Kpt2017:
//                KPT2017(pointsBegin, pointsEnd, back_inserter(spanner));
//                break;
//            case BoundedDegreePlaneSpannerAlgorithm::Degree3:
//                DEG3(pointsBegin, pointsEnd, back_inserter(spanner));
//                break;
            case BoundedDegreePlaneSpannerAlgorithm::AlgorithmLast:
                assert(false);
        }

        number_t runtime = tim.stop();

        BoundedDegreeSpannerResult result(distributionType, distribution,
                                          algorithm, runtime,
                                          points.begin(), points.end(),
                                          spanner.begin(), spanner.end(),
                                          !measureStretchFactor);
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
    template<typename DistributionSubTypeEnum>
    void BoundedDegreePlaneSpannerAlgorithmLoop(const DistributionType distributionType,
                                                const DistributionSubTypeEnum distribution,
                                                const vector<Point>& points,
                                                ofstream& expOut,
                                                const bool measureStretchFactor = true ){
        const size_t n = points.size();

        cout<< "Starting trial..."<<endl<<endl;

        for(int alg=BoundedDegreePlaneSpannerAlgorithm::AlgorithmFirst;
            alg != BoundedDegreePlaneSpannerAlgorithm::AlgorithmLast; ++alg ) {
            auto algorithm = static_cast<BoundedDegreePlaneSpannerAlgorithm>(alg);
            list<pair<size_t, size_t> > spanner;

            BoundedDegreePlaneSpannerTest(algorithm, distributionType, distribution, points, expOut, measureStretchFactor);

            ++EXP_COUNT;
        }

        cout<< "Finished trial...\n"
            << "-----------------"<<endl;
    }

    // Generates a random point set from the given distribution and runs an experiment with multiple graphs on the point set
    void generateRandomPointSet(const SyntheticDistribution dist,
                                const size_t n,
                                const double width,
                                vector<Point>& points){
        PointGenerator_2 pointGenerator;

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
//            case UniformOnCircle:
//                pointGenerator.generatePointsOnACircle(n, width, points);
//                break;
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
            case SyntheticDistributionLast:
            default:
                assert(!"Invalid distribution type!");
        }

    }

    // An experiment from n_start to n_end with a single distribution
    void SyntheticExperimentInputSizeLoop(SyntheticDistribution dist,
                                          size_t n_start, size_t n_end, size_t increment, ofstream& expOut, bool measureStretchFactor = true ) {
        measureStretchFactor = false;
        for (size_t n = n_start; n <= n_end; n += increment) {
            // SET POINTS
            vector<Point> points;
            generateRandomPointSet(dist, n, INPUT_WIDTH, points);
            BoundedDegreePlaneSpannerAlgorithmLoop(DistributionType::Synthetic,dist,points,expOut,measureStretchFactor);
        }
    }

    // An experiment of numRuns trials for a single distribution
    void SyntheticExperimentRepetitionLoop(SyntheticDistribution dist,
                                           size_t numRuns,
                                           size_t n_start,
                                           size_t n_end,
                                           size_t increment,
                                           ofstream& expOut ) {
        string distName = SYNTHETIC_DISTRIBUTION_NAMES.at(dist);
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

        for(int dist=SyntheticDistributionFirst; dist != SyntheticDistributionLast; ++dist ) {
            auto distributionType = static_cast<SyntheticDistribution>(dist);
            string distName = SYNTHETIC_DISTRIBUTION_NAMES.at(dist);

            cout<< "!! Starting  "<< distName << "distribution trials !!\n"<<endl<<endl;

            SyntheticExperimentRepetitionLoop(distributionType, numRuns, n_start, n_end, increment, expOut);

            cout<< "!! Ending  "<< distName << "distribution trials !!\n"
                << "-------------------------------------------"<<endl;
        }
    }

    // An experiment of numRuns trials for each distribution
    void SyntheticExperiment(size_t numRuns, size_t n_start, size_t n_end, size_t increment ) {

        // get unix timestamp to use as experiment file name
        string filename = "synthetic-";
        auto time = std::time(nullptr);
        filename += std::to_string(time)
                    + OUTPUT_EXTENSION;

        ofstream expOut;
        expOut.open(filename,ios_base::out);

        if(!expOut.is_open()) assert(!"ERROR OPENING OUTPUT FILE\n\n");

        expOut << "distribution" << DELIMITER
               << "n" << DELIMITER
               << "spannerAlg" << DELIMITER
               << "runtime" << DELIMITER
               << "degree" << DELIMITER
               << "degreeAvg" << DELIMITER
               << "stretchFactor" << DELIMITER
               << "lightness" << DELIMITER
               << "\n";

        SyntheticExperimentDistributionLoop(numRuns,n_start,n_end,increment,expOut);

        expOut.close();
    }



    void ExperimentFromConfigurationFile(string configFilename, size_t numRuns=1) {
        using std::string, std::vector;
        namespace pt = boost::property_tree;

        configFilename = INPUT_DATA_DIRECTORY + configFilename;
        pt::ptree config;
        pt::read_xml(configFilename,config);


        // get unix timestamp to use as experiment file name
        auto time = std::time(nullptr);
        string filename = "real-";
        filename += std::to_string(time)
                 + OUTPUT_EXTENSION;

        ofstream expOut;
        expOut.open(filename,ios_base::out);

        if(!expOut.is_open()) assert(!"ERROR OPENING OUTPUT FILE\n\n");

        unsigned pointSetID = 0;

        for( auto pointset : config.get_child("pointsets") ) {

            string filename = pointset.second.get_child("filename").data(),
                   fullname = INPUT_DATA_DIRECTORY + filename,
                   filenameNoExtension = filename;
            boost::erase_all(filenameNoExtension, ".xy");

            vector<Point> P;
            PointGenerator_2 generator;
            generator.loadFromFile(fullname,P);
            const index_t n = P.size();

            string pointsetName = pointset.second.get_child("nicename").data();
            REAL_POINTSET_NAMES.push_back(pointsetName);

            cout<< "!! Starting  "<< pointsetName << " trials !!\n"
                << "Added "<< n <<" points from file\n\n";

            for (size_t trial = 0; trial < numRuns; ++trial) {
                BoundedDegreePlaneSpannerAlgorithmLoop(DistributionType::Real, pointSetID, P, expOut, trial == 0);
            }

            cout<< "!! Ending  "<< pointsetName << " trials !!\n"
                << "-------------------------------------------"<<endl;

            pointSetID++;
        }

        expOut.close();
    }


} // spanners

#endif //SPANNERS_EXPERIMENT_H
