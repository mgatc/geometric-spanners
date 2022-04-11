#ifndef SPANNERS_EXPERIMENT_H
#define SPANNERS_EXPERIMENT_H

#include <iostream>
#include <list>
#include <optional>
#include <random>
#include <string>
#include <utility>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "libspanner/BoundedDegreePlaneSpanners.h"

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

    const bool USE_EXPERIMENTAL_STRETCH_FACTOR = true;

    const char * INPUT_DATA_DIRECTORY = "../input/";
    const char * OUTPUT_DATA_DIRECTORY = "./";
    const char *OUTPUT_EXTENSION = ".csv";
    const char *DELIMITER = ",";

    using std::to_string;

    const spanner::number_t INPUT_WIDTH = 1000;

    size_t EXP_COUNT = 0;


    /**
     * Constructs on spanner using the given algorithm on the given point set,
     * then writes the result to the expOut file. Due to its
     */
    template <typename DistributionSubTypeEnum>
    void BoundedDegreePlaneSpannerTest(const spanner::BoundedDegreePlaneSpannerAlgorithm algorithm,
                                       const spanner::DistributionType distributionType,
                                       const DistributionSubTypeEnum distribution,
                                       const spanner::bdps::input_t &points,
                                       std::ofstream& expOut,
                                       bool measureStretchFactor = true) {
        using namespace std;

        spanner::number_t t(1.25);

        spanner::bdps::output_t output;

        auto pointsBegin = points.begin(),
                pointsEnd = points.end();

        spanner::Timer tim;

        switch (algorithm) {
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bgs2005:
                spanner::BGS2005(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Lw2004:
                spanner::LW2004(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bsx2009:
                spanner::BSX2009(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Kpx2010:
                spanner::KPX2010(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Kx2012:
                spanner::KX2012(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bcc2012_7:
                spanner::BCC2012<7>(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bcc2012_6:
                spanner::BCC2012<6>(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bhs2018:
                spanner::BHS2018(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bkpx2015:
                spanner::BKPX2015(points,output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bghp2010:
                spanner::BGHP2010(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Kpt2017:
                spanner::KPT2017(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Degree3:
                spanner::DEG3(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::BdpsAlgorithmLast:
            default:
//                std::cout<< "Skip"<<std::endl;
                tim.stop();
                return;
        }

        spanner::number_t runtime = tim.stop();



        BoundedDegreePlaneSpannerResult result(distributionType, distribution,
                                          algorithm, runtime,
                                          points.begin(), points.end(),
                                          output.begin(), output.end(),
                                          !measureStretchFactor);
        std::cout << result;
        expOut << result;

        if(!result.verify(t)) {
            std::string filename("breaks");
            filename += spanner::bdps::ALGORITHM_NAMES.at(algorithm);
            filename += ".xy";
            spanner::writePointsToFile(points.begin(),points.end(),filename);


            cpptex::GraphPrinter tikz("/tmp/deg3", points.begin(), points.end());

            cpptex::GraphPrinter::OptionsList edgeOptions = { // active edge options
                    {"color",      tikz.activeEdgeColor},
                    {"line width", to_string(tikz.inactiveEdgeWidth/2)}
            };

            tikz.drawEdges(output.begin(),output.end(),points,edgeOptions);



            tikz.drawVerticesWithInfo(points.begin(), points.end(), tikz.activeVertexOptions);


            tikz.display();

            std::cout<< "\n"
                << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                << "!!!!!!!!!!!!!!!! ALGORITHM ERROR !!!!!!!!!!!!!!!!!!!!"
                << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";

            char userContinues = 'x';
            while( userContinues != 'y' && userContinues != 'n' ) {
                std::cout<< "\n Do you wish to continue the experiment? (y/n) ";
                cin>>userContinues;
            }
            if(userContinues=='n') assert(!"Experiment terminated by user.");
        }


        ++EXP_COUNT;
    }

    // An experiment featuring multiple graphs on the same point set
    template<typename DistributionSubTypeEnum>
    void BoundedDegreePlaneSpannerAlgorithmLoop(const spanner::DistributionType distributionType,
                                                const DistributionSubTypeEnum distribution,
                                                const std::vector<spanner::Point>& points,
                                                std::ofstream& expOut,
                                                const bool measureStretchFactor = true ){
        const size_t n = points.size();

        std::cout<< "Starting trial...\n"<<std::endl;

        for(int alg=spanner::BoundedDegreePlaneSpannerAlgorithm::BdpsAlgorithmFirst;
            alg != spanner::BoundedDegreePlaneSpannerAlgorithm::BdpsAlgorithmLast; ++alg ) {
            auto algorithm = static_cast<spanner::BoundedDegreePlaneSpannerAlgorithm>(alg);
            std::list<std::pair<size_t, size_t> > spanner;

            BoundedDegreePlaneSpannerTest(algorithm, distributionType, distribution, points, expOut, measureStretchFactor);

            ++EXP_COUNT;
        }

        std::cout<< "Finished trial...\n"
            << "-----------------"<<std::endl;
    }

    // Generates a random point set from the given distribution and runs an experiment with multiple graphs on the point set
    void generateRandomPointSet(const spanner::SyntheticDistribution dist,
                                const size_t n,
                                const double width,
                                std::vector<spanner::Point>& points){
        using namespace spanner;

        PointGenerator_2 pointGenerator;

        switch(dist) {
//            case UniformInsideSquare:
//                pointGenerator.insideSquare(n,width,points);
//                break;
//            case UniformInsideDisc:
//                pointGenerator.insideDisc(n,width,points);
//                break;
////            case UniformOnSquare:
////                pointGenerator.onSquare(n,width,points);
////                break;
////            case UniformOnCircle:
////                pointGenerator.onCircle(n, width, points);
////                break;
//            case NormalInsideSquare:
//                pointGenerator.insideSquareNormal(n,1,points);
//                break;
            case NormalClustersInsideSquare:
                pointGenerator.insideSquareNormal(n/5, 5, points);
                break;
//            case ContiguousGrid:
//                pointGenerator.contiguousOnGrid(n, points);
//                break;
//            case UniformRandomGrid:
//                pointGenerator.randomOnGrid(n, points);
//                break;
//            case UniformInsideAnnulus:
//                pointGenerator.insideAnnulus(n, width, width*0.8, points);
//                break;
//            case Galaxy:
//                pointGenerator.inGalaxy(n, 5, points);
//                break;
//            case ConvexHullInDisc:
//                pointGenerator.onConvexHullInDisc(n, width, points);
//                break;
            case SyntheticDistributionLast:
            default:
                assert(!"Invalid distribution type!");
        }

    }

    // An experiment from n_start to n_end with a single distribution
    void SyntheticExperimentInputSizeLoop(spanner::SyntheticDistribution dist,
                                          size_t n_start, size_t n_end, size_t increment, std::ofstream& expOut, bool measureStretchFactor = true ) {
        //measureStretchFactor = false;
        for (size_t n = n_start; n <= n_end; n += increment) {
            // SET POINTS
            std::vector<spanner::Point> points;
            generateRandomPointSet(dist, n, INPUT_WIDTH, points);
            BoundedDegreePlaneSpannerAlgorithmLoop(spanner::DistributionType::Synthetic,dist,points,expOut,measureStretchFactor);
        }
    }

    // An experiment of numRuns trials for a single distribution
    void SyntheticExperimentRepetitionLoop(spanner::SyntheticDistribution dist,
                                           size_t numRuns,
                                           size_t n_start,
                                           size_t n_end,
                                           size_t increment,
                                           std::ofstream& expOut ) {
        std::string distName = spanner::SYNTHETIC_DISTRIBUTION_NAMES.at(dist);
        srand(0);
        for (size_t trial = 0; trial < numRuns; ++trial) {
            std::cout<< "Starting trial "<< trial << " of "<<distName<<"\n\n";
            SyntheticExperimentInputSizeLoop(dist, n_start, n_end, increment, expOut);
            std::cout<<"\n"<<std::endl;
        }
    }

    void SyntheticExperimentDistributionLoop(size_t numRuns,
                                             size_t n_start,
                                             size_t n_end,
                                             size_t increment,
                                             std::ofstream& expOut) {
        const spanner::number_t width = 10;

        for(int dist=spanner::SyntheticDistributionFirst; dist != spanner::SyntheticDistributionLast; ++dist ) {
            auto distributionType = static_cast<spanner::SyntheticDistribution>(dist);
            std::string distName = spanner::SYNTHETIC_DISTRIBUTION_NAMES.at(dist);

            std::cout<< "!! Starting  "<< distName << "distribution trials !!\n"<<std::endl;

            SyntheticExperimentRepetitionLoop(distributionType, numRuns, n_start, n_end, increment, expOut);

            std::cout<< "!! Ending  "<< distName << "distribution trials !!\n"
                << "-------------------------------------------"<<std::endl;
        }
    }

    // An experiment of numRuns trials for each distribution
    void SyntheticExperiment(size_t numRuns, size_t n_start, size_t n_end, size_t increment ) {

        // get unix timestamp to use as experiment file name
        std::string filename = "synthetic-";
        auto time = std::time(nullptr);
        filename += std::to_string(time)
                    + OUTPUT_EXTENSION;

        std::ofstream expOut;
        expOut.open(filename,std::ios_base::out);

        if(!expOut.is_open()) assert(!"ERROR OPENING OUTPUT FILE\n\n");

        expOut << "distribution" << DELIMITER
               << "n" << DELIMITER
               << "spannerAlg" << DELIMITER
               << "runtime" << DELIMITER
               << "degree" << DELIMITER
               << "avgDegreePerPoint" << DELIMITER
               << "avgStretchFactor" << DELIMITER
               << "lightness" << DELIMITER
               << "\n";

        SyntheticExperimentDistributionLoop(numRuns,n_start,n_end,increment,expOut);

        expOut.close();
    }



    void ExperimentFromConfigurationFile(std::string configFilename, size_t numRuns=1) {
        using std::string, std::vector;
        namespace pt = boost::property_tree;

        configFilename = configFilename;
        pt::ptree config;
        pt::read_xml(configFilename,config);


        // get unix timestamp to use as experiment file name
        auto time = std::time(nullptr);
        string filename = "real-";
        filename += std::to_string(time)
                 + OUTPUT_EXTENSION;

        std::ofstream expOut;
        expOut.open(filename,std::ios_base::out);

        if(!expOut.is_open()) assert(!"ERROR OPENING OUTPUT FILE\n\n");

        unsigned pointSetID = 0;

        for( auto pointset : config.get_child("pointsets") ) {

            string filename = pointset.second.get_child("filename").data(),
                   fullname = INPUT_DATA_DIRECTORY + filename,
                   filenameNoExtension = filename;
            boost::erase_all(filenameNoExtension, ".xy");

            std::cout<<fullname<<std::endl;

            vector<spanner::Point> P;
            spanner::PointGenerator_2 generator;
            generator.loadFromFile(fullname,P);
            const spanner::index_t n = P.size();

            string pointsetName = pointset.second.get_child("nicename").data();
            spanner::REAL_POINTSET_NAMES.push_back(pointsetName);

            std::cout<< "!! Starting  "<< pointsetName << " trials !!\n"
                << "Added "<< n <<" points from file\n\n";

            for (size_t trial = 0; trial < numRuns; ++trial) {
                BoundedDegreePlaneSpannerAlgorithmLoop(spanner::DistributionType::Real, pointSetID, P, expOut, trial == 0);
            }

            std::cout<< "!! Ending  "<< pointsetName << " trials !!\n"
                << "-------------------------------------------"<<std::endl;

            pointSetID++;
        }

        expOut.close();
    }


} // bdps_experiment

#endif //SPANNERS_EXPERIMENT_H
