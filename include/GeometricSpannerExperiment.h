#ifndef GEOMETRICSPANNEREXPERIMENT_H
#define GEOMETRICSPANNEREXPERIMENT_H

#include "CgalComponents.h"
#include "GeometricSpannerSolutionChecker.h"
#include <typeinfo>
#include <cxxabi.h>
#include <cstdlib>
#include <memory>
#include <string>

namespace GeometricSpanner {

const std::string DELIMITER = ",";
const double EPSILON = 0.001;

const std::string EXP_DIR = "experiments";
const std::string RES_DIR = "results";
const std::string FILE_EXT = ".csv" ;
const std::string PROB_DIR = "problems";
const std::string SYNTH_DIR = "synthetic";

bool makeDirectory( std::string dirName );
std::string demangle(char const* mangled);

/**************
 * GeometricSpannerProblem *
 **************/

    class GeometricSpannerProblem {
    public:
        std::string fileName;
        std::string niceName;
        std::optional<int> bound;

        GeometricSpannerProblem() {} // default constructor
        GeometricSpannerProblem( std::string niceName, std::string fileName, std::optional<int> bound = std::nullopt ); // construct from file
        std::vector<Point> getPoints();

    };



/***********************
 * GeometricSpannerProblemSynthetic *
 ***********************/

    template< class G >
    class GeometricSpannerProblemSynthetic : public GeometricSpannerProblem {
    public:
        std::string syntheticID;
        std::list<GeometricSpannerProblem> problems;

        GeometricSpannerProblemSynthetic( const std::vector< pair< int, int > >& problemDesc, const std::string& id )
            : syntheticID( id ) {
            std::cout << "Generating synthetic problem files...";
            for( auto desc : problemDesc ) {
                std::cout << ".";

                std::fstream file;
                std::string instanceID( std::to_string( desc.first ) + "_" + std::to_string( desc.second ) );
                std::string problemID( syntheticID + "_" + instanceID );

                // setup filesystem
                makeDirectory( "./" + EXP_DIR );
                makeDirectory( "./" + EXP_DIR + "/" + PROB_DIR );
                makeDirectory( "./" + EXP_DIR + "/" + PROB_DIR + "/" + SYNTH_DIR );

                // set path of output file
                std::string relPath(
                      "./" + EXP_DIR
                    +  "/" + PROB_DIR
                    +  "/" + SYNTH_DIR
                    +  "/" + problemID + ".txt"
                );
                //std::cout<<relPath<<std::endl;

                // create random point generator with second int as argument
                G g1( desc.second );

                file.open( relPath, std::fstream::trunc | std::fstream::out );
                //std::copy_n( g1, desc.first, back_inserter( points ) );

                if( !file.is_open() )
                    std::cout << "File error!" << std::endl;

                for( unsigned i=0; i<desc.first; i++ )
                    file << *g1++ << std::endl;

                file.close();

                // add new GeometricSpannerProblem list of problems
                problems.push_back( GeometricSpannerProblem( problemID, relPath ) );
            }
            std::cout << "done." << std::endl;
        }

    };



/***********************
 * GeometricSpannerExperimentResult *
 ***********************/

    class GeometricSpannerExperimentResult {
    public:
        std::string experiment;
        std::string problem; // problem
        std::string algorithm; // algorithm
        std::ofstream file; // file
        int bound; // bounding size
        int problemSize; // numPoints
        long int runtime; // runtime in microseconds
        int solutionSize; // solution size
        bool verified;


        GeometricSpannerExperimentResult(
            const std::string& experimentName,
            const std::string& problemName,
            const std::string& algorithmName,
            int bound,
            int problemSize,
            long int runtime,
            int solutionSize,
            bool verified
        );
        void writeToFile();
        bool passes();

    protected:

        void appendToFile( std::string );
    };

   inline ostream& operator<<( ostream& os, const GeometricSpannerExperimentResult& res ) {
        std::string checkResult = res.verified ? "pass" : "fail";

        os << res.experiment << DELIMITER
           << res.problem << DELIMITER
           << res.algorithm << DELIMITER
           << res.bound << DELIMITER
           << res.problemSize << DELIMITER
           << res.runtime << DELIMITER
           << res.solutionSize << DELIMITER
           << checkResult;

        return os;
    }



/***************************
 * GeometricSpannerExperiment *
 ***************************/

    class GeometricSpannerExperiment {
    public:
        std::string experimentID;
        std::list<GeometricSpannerProblem> problems;
        GeometricSpannerSolutionChecker checker;
        int experimentsTotal = 0;
        int experimentsFails = 0;

        GeometricSpannerExperiment( double epsilon = EPSILON );

        void addProblem( GeometricSpannerProblem prob ); // Add problem from GeometricSpannerProblem

        template< class G >                 // Add synthetic problem(s) for GeometricSpannerProblemSynthetic
        void addProblem( GeometricSpannerProblemSynthetic<G> synthProb ) {
            // add problems in prob.problems to problems
            for( auto p : synthProb.problems )
                addProblem( GeometricSpannerProblem( p.niceName, p.fileName ) );
        }

        // run experiment with given number of trials on the given algorithms
        template< class Alg, class ...Algs >
        void run( int numTrials=1 ) {
            std::cout << "\nBeginning Experiment...\n" << std::endl;
            for( unsigned trial = 0; trial < numTrials; trial++ ) {
                std::cout << "Starting trial " << trial << std::endl;
                for( auto problem : problems ) {
                    startTrial< Alg, Algs... >( problem );
                }
                std::cout << "Trial complete." << std::endl;
                std::cout<<std::endl;
            }
            std::cout << "\nExperiment complete." << std::endl;
        }

        // perform a single run on a single alg
        template< class Alg >      // Start a single trial with algorithm specified in template arguments
        void startTrial( GeometricSpannerProblem prob ) {
            //std::cout << ".";
            // get class name of Alg
            std::string algName = demangle( typeid( Alg ).name() );
            int bound = prob.bound ? *prob.bound : -1;

            // load problem from file
            std::vector<Point> P = prob.getPoints();
            std::list<Point> C; // list of disk centers

            auto start = chrono::high_resolution_clock::now();      // start time

            Alg( P, C );                                            // run alg

            auto stop = chrono::high_resolution_clock::now();       // stop time
            auto runtime = chrono::duration_cast<chrono::microseconds>(stop - start);

            // build result
            GeometricSpannerExperimentResult result(
                    experimentID,
                    prob.niceName,
                    algName, // algorithm name
                    bound, // bound, if present
                    P.size(),
                    runtime.count(),
                    C.size(),
                    checker.check( P, C )
            );

            std::cout << result << std::endl;

            result.writeToFile();
            if( !result.passes() )
                std::cout << "\nfailed: " << result << std::endl;
        }
        // perform trial on next alg in list
        template< class Alg, class NextAlg, class ...Algs >
        void startTrial( GeometricSpannerProblem prob ) {
            startTrial< Alg >( prob );              // Single trial for first algorithm in template arguments
            startTrial< NextAlg, Algs... >( prob ); // Call this function for the rest of the algorithms
        }

    protected:

        double epsilon;

    };

    template< class A, class G >
    void simpleExperiment( int n, double a ) {

        GeometricSpanner::GeometricSpannerExperiment exp;

        std::vector< std::pair< int, int > > specs {
            std::make_pair( n, a )
        };

        GeometricSpanner::GeometricSpannerProblemSynthetic< G > prob(
            specs,
            "single"
        );

        exp.addProblem( prob );

        exp.run<
            A
        >(1);

    }



/*************************
 * Random_points_in_ellipse *
 *************************/

    template< typename P, typename C >
    class Random_points_in_ellipse_2 : public CGAL::Random_points_in_disc_2< P,C > {
    public:

        Random_points_in_ellipse_2( double r, double x, double y )
            : g(r), x(x), y(y), r(r) {
            updatePoint();
        }

        P& operator*() {
            return p;                      // return the existing member point
        }

        const P& operator*() const {
            return p;                      // return the existing member point
        }

        Random_points_in_ellipse_2< P,C > operator++() {
            // update, then return updated object
            updatePoint();
            return *this;
        }

        Random_points_in_ellipse_2< P,C > operator++( int ) {
            // copy original to return, but update original
            Random_points_in_ellipse_2< P,C > temp( r, x, y );

            updatePoint();
            return temp;                           // return copy
        }

    protected:
        double x;
        double y;
        double r;
        CGAL::Random_points_in_disc_2< P,C > g;
        P p;

        void updatePoint() {
            Point temp = *(g++);
            p = P( x*temp.x(), y*temp.y() );
        }
    };

    template< typename P, typename C >
    class Random_points_in_ellipse_2_hor : public Random_points_in_ellipse_2< P,C > {
    public:
        Random_points_in_ellipse_2_hor( double r, double f=2 )
            : Random_points_in_ellipse_2< P,C >( r, f, 1 ) {
        }
    };

    template< typename P, typename C >
    class Random_points_in_ellipse_2_vert : public Random_points_in_ellipse_2< P,C > {
    public:
        Random_points_in_ellipse_2_vert( double r, double f=2 )
            : Random_points_in_ellipse_2< P,C >( r, 1, f ) {
        }
    };

}

#endif
