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

#include "libspanner/BoundedDegreePlaneSpanners.h"

#include "tools/printers/LatexPrinter.h"
#include "tools/printers/PgfplotPrinter.h"
#include "tools/printers/TablePrinter.h"

#include "tools/Metrics.h"
#include "tools/PointGenerators.h"
#include "tools/Results.h"
#include "tools/Utilities.h"


namespace bdps_experiment {

    const bool USE_EXPERIMENTAL_STRETCH_FACTOR = true;

    const char * INPUT_DATA_DIRECTORY = "../input/";
    const char * OUTPUT_DATA_DIRECTORY = "./";
    const char *OUTPUT_EXTENSION = ".csv";
    const char *DELIMITER = ",";

    using std::to_string;

    const number_t INPUT_WIDTH = 1000;

    size_t EXP_COUNT = 0;


    /**
     * Constructs on spanner using the given algorithm on the given point set,
     * then writes the result to the expOut file. Due to its
     */
    template <typename DistributionSubTypeEnum>
    void BoundedDegreePlaneSpannerTest(const spanner::BoundedDegreePlaneSpannerAlgorithm algorithm,
                                       const DistributionType distributionType,
                                       const DistributionSubTypeEnum distribution,
                                       const spanner::bdps::input_t &points,
                                       ofstream& expOut,
                                       bool measureStretchFactor = true) {
        using namespace std;

        spanner::bdps::output_t output;

        auto pointsBegin = points.begin(),
                pointsEnd = points.end();

        Timer tim;

        switch (algorithm) {
//            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bgs2005:
//                spanner::BGS2005(points, output);
//                break;
//            case spanner::BoundedDegreePlaneSpannerAlgorithm::Lw2004:
//                spanner::LW2004(points, output);
//                break;
//            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bsx2009:
//                spanner::BSX2009(points, output);
//                break;
//            case spanner::BoundedDegreePlaneSpannerAlgorithm::Kpx2010:
//                spanner::KPX2010(points, output);
//                break;
//            case spanner::BoundedDegreePlaneSpannerAlgorithm::Kx2012:
//                spanner::KX2012(points, output);
//                break;
//            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bcc2012_7:
//                spanner::BCC2012<7>(points, output);
//                break;
//            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bcc2012_6:
//                spanner::BCC2012<6>(points, output);
//                break;
//            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bhs2018:
//                spanner::BHS2018(points, output);
//                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bkpx2015:
                spanner::BKPX2015(points,output);
                break;
//            case spanner::BoundedDegreePlaneSpannerAlgorithm::Bghp2010:
//                spanner::BGHP2010(points, output);
//                break;
//            case spanner::BoundedDegreePlaneSpannerAlgorithm::Kpt2017:
//                spanner::KPT2017(points, output);
//                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::Degree3:
                spanner::DEG3(points, output);
                break;
            case spanner::BoundedDegreePlaneSpannerAlgorithm::AlgorithmLast:
                assert(false);
        }

        number_t runtime = tim.stop();



        BoundedDegreeSpannerResult result(distributionType, distribution,
                                          algorithm, runtime,
                                          points.begin(), points.end(),
                                          output.begin(), output.end(),
                                          !measureStretchFactor);
        cout << result;
        expOut << result;

        if(!result.verify()) {
            std::string filename("breaks");
            filename += spanner::bdps::ALGORITHM_NAMES.at(algorithm);
            filename += ".xy";
            writePointsToFile(points.begin(),points.end(),filename);


            GraphPrinter tikz("/tmp/", "deg3");
            tikz.autoscale(points.begin(), points.end());

            GraphPrinter::OptionsList edgeOptions = { // active edge options
                    {"color",      tikz.activeEdgeColor},
                    {"line width", to_string(tikz.inactiveEdgeWidth/2)}
            };

            tikz.drawEdges(output.begin(),output.end(),points,edgeOptions);



            tikz.drawVerticesWithInfo(points.begin(), points.end(), tikz.activeVertexOptions);


            tikz.display();

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

        for(int alg=spanner::BoundedDegreePlaneSpannerAlgorithm::AlgorithmFirst;
            alg != spanner::BoundedDegreePlaneSpannerAlgorithm::AlgorithmLast; ++alg ) {
            auto algorithm = static_cast<spanner::BoundedDegreePlaneSpannerAlgorithm>(alg);
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
    void SyntheticExperimentInputSizeLoop(SyntheticDistribution dist,
                                          size_t n_start, size_t n_end, size_t increment, ofstream& expOut, bool measureStretchFactor = true ) {
        //measureStretchFactor = false;
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
               << "avgStretchFactor" << DELIMITER
               << "lightness" << DELIMITER
               << "\n";

        SyntheticExperimentDistributionLoop(numRuns,n_start,n_end,increment,expOut);

        expOut.close();
    }



    void ExperimentFromConfigurationFile(string configFilename, size_t numRuns=1) {
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

        ofstream expOut;
        expOut.open(filename,ios_base::out);

        if(!expOut.is_open()) assert(!"ERROR OPENING OUTPUT FILE\n\n");

        unsigned pointSetID = 0;

        for( auto pointset : config.get_child("pointsets") ) {

            string filename = pointset.second.get_child("filename").data(),
                   fullname = INPUT_DATA_DIRECTORY + filename,
                   filenameNoExtension = filename;
            boost::erase_all(filenameNoExtension, ".xy");

            cout<<fullname<<endl;

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


} // bdps_experiment

#endif //SPANNERS_EXPERIMENT_H
