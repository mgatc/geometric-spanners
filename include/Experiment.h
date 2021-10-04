#ifndef SPANNERS_EXPERIMENT_H
#define SPANNERS_EXPERIMENT_H

#include <iostream>
#include <list>
#include <optional>
#include <utility>


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

#include "tools/Metrics.h"
#include "tools/PointGenerators.h"
#include "tools/Results.h"
#include "tools/Utilities.h"


namespace spanners {

    using std::to_string;

    void SingleTrial (const vector<Point>&, DistributionType, ofstream&, bool);
    void SyntheticTrial(const size_t, DistributionType, const double, ofstream& );



    void SyntheticExperiment(size_t numRuns, size_t n_start, size_t n_end, size_t increment ) {

        // get unix timestamp to use as experiment file name
        auto time = std::time(nullptr);
        auto filename = to_string(time) + ".csv";

        ofstream expOut;
        expOut.open(filename,ios_base::out);

        if(!expOut.is_open()) assert(!"ERROR OPENING OUTPUT FILE\n\n");

        const string DELIMITER = ",";

        expOut << "distribution" << DELIMITER
               << "n" << DELIMITER
               << "spannerAlg" << DELIMITER
               << "runtime" << DELIMITER
               << "degree" << DELIMITER
               << "degreeAvg" << DELIMITER
               << "stretchFactor" << DELIMITER
               << "lightness" << DELIMITER
               << "\n";

        const number_t width = 10;

        for( int dist=DistributionTypeFirst; dist!=DistributionTypeLast; ++dist ) {
            auto distributionType = static_cast<DistributionType>(dist);
            string distName = DISTRIBUTION_NAMES.at(dist);

            cout<< "!! Starting  "<< distName << "distribution trials !!\n"
                << "NOTE: one extra trial is performed because trial 0 will be thrown out!"<<endl<<endl;

            for (size_t trial = 0; trial <= numRuns; ++trial) {
                cout<< "Starting trial "<< trial << "..."<<endl<<endl;
                for (size_t n = n_start; n <= n_end; n += increment) {
                    SyntheticTrial(n, distributionType, width, expOut);
                }
                cout<<"\n\n";
            }

            cout<< "!! Ending  "<< distName << "distribution trials !!\n"
                << "-------------------------------------------"<<endl;

            if(!expOut.is_open()) assert(!"ERROR! OUTPUT FILE CLOSED!\n\n");
        }
        expOut.close();
    }


    void PlaneSpannerTest(const vector<Point> &points,
                          const DistributionType dist,
                          const Algorithm algorithm,
                          ofstream& expOut,
                          bool lite = false ) {
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
        BoundedDegreeSpannerResult result(dist, algorithm, runtime, points.begin(), points.end(), spanner.begin(), spanner.end(), lite);
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
    }

    void SingleTrial (const vector<Point>& points, const DistributionType dist, ofstream& expOut, bool lite = false ){
        cout<< "Starting trial..."<<endl<<endl;

        for(int alg=Algorithm::AlgorithmFirst;
          alg!=Algorithm::AlgorithmLast; ++alg ) {
            PlaneSpannerTest(points, dist, static_cast<Algorithm>(alg), expOut, lite);
        }

        cout<< "Finished trial...\n"
            << "-----------------"<<endl;
    }


    void SyntheticTrial(const size_t n, DistributionType dist, const double width, ofstream& expOut  ){

        // SET POINTS
        vector<Point> points;

        switch(dist) {
        case UniformInsideSquare:
            generatePointsInsideASquare(n,width,points);
            break;
        case UniformInsideDisc:
            generatePointsInsideADisc(n,width,points);
            break;
        case UniformOnDisc:
            generatePointsOnADisc(n,width,points);
            break;
        case NormalInsideSquare:
            generatePointsInsideASquareNormal(n,1,points);
            break;
        case NormalClustersInsideSquare:
            generatePointsInsideASquareNormal(static_cast<size_t>(n/pow(n,(1/3))),
                                              static_cast<size_t>(pow(n,(1/3))),points);
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
        case DistributionTypeLast:
        default:
            assert(!"Invalid distribution type!");
        }

        SingleTrial(points, dist,expOut );
    }

} // spanners

#endif //SPANNERS_EXPERIMENT_H
