#include "GeometricSpannerExperiment.h"
#include <chrono>
#include <fstream>

#include <string>
#include <sys/stat.h>

#include <time.h>
#include <bits/stdc++.h>



namespace GeometricSpanner {

/***********************
 * Nonmember Functions *
 ***********************/

    bool makeDirectory( std::string dirName ) {
        // Adapted from https://www.geeksforgeeks.org/create-directoryfolder-cc-program/

        // Creating a directory, if necessary
        if( mkdir( dirName.c_str(), 0777 ) == -1 ) {
            if( errno != 17 ) {// file exists error, we don't care
                std::cerr << "Error :  " << strerror( errno ) << std::endl;
                return false;
            }
        }
        return true;
    }

    std::string demangle(char const* mangled) {
        auto ptr = std::unique_ptr<char, decltype(& std::free)>{
            abi::__cxa_demangle(mangled, nullptr, nullptr, nullptr),
            std::free
        };
        return {ptr.get()};
    }




/**************
 * GeometricSpannerProblem *
 **************/

    GeometricSpannerProblem::GeometricSpannerProblem( std::string niceName, std::string fileName, std::optional<int> bound )
        : niceName( niceName ), fileName( fileName ), bound( bound ) {

    }

    std::vector<Point> GeometricSpannerProblem::getPoints() {
        std::vector<Point> points;
        std::ifstream input( fileName );

        long double x,y;
        while(input >> x >> y )
            points.emplace_back(Point(x,y));

        input.close();
        return points;
    }




/***********************
 * GeometricSpannerExperimentResult *
 ***********************/

    GeometricSpannerExperimentResult::GeometricSpannerExperimentResult(
            const std::string& experimentName,
            const std::string& problemName,
            const std::string& algorithmName,
            int bound,
            int problemSize,
            long int runtime,
            int solutionSize,
            bool verified

        ):
            experiment( experimentName ),
            problem( problemName ),
            algorithm( algorithmName ),
            bound( bound ),
            problemSize( problemSize ),
            runtime( runtime ),
            solutionSize( solutionSize ),
            verified( verified )
    {

    }

    bool GeometricSpannerExperimentResult::passes() {
        return verified;
    }

    void GeometricSpannerExperimentResult::writeToFile() {
        // Make directories if necessary
        makeDirectory( "./" + EXP_DIR );
        makeDirectory( "./" + EXP_DIR + "/" + RES_DIR );
//        makeDirectory( "./" + EXP_DIR + "/" + RES_DIR + "/" + experiment + "/" + problem );

//        std::string dedicatedFilePath = "./" + EXP_DIR
//                            +  "/" + RES_DIR
//                            +  "/" + experiment
//                            +  "/" + problem
//                            +  "/" + algorithm + FILE_EXT;
//        appendToFile( dedicatedFilePath );

        std::string masterFilePath = "./" + EXP_DIR
                            +  "/" + RES_DIR
                            +  "/" + experiment + FILE_EXT;

        appendToFile( masterFilePath );
    }

    void GeometricSpannerExperimentResult::appendToFile( std::string path ) {
        // Open or create file
        file.open( path, std::fstream::app );

        if( !file.is_open() )
            std::cout << "File error!" << std::endl;

        std::string checkResult = verified ? "pass" : "fail";

        file << *this << std::endl;

        file.close();
    }



/***************************
 * GeometricSpannerExperiment *
 ***************************/

    GeometricSpannerExperiment::GeometricSpannerExperiment( double epsilon )
        : experimentID( std::to_string( static_cast<long int>( time( NULL ) ) ) ),
          checker( epsilon ) {

        std::cout << "Experiment ID: " << experimentID << std::endl;
        std::cout << "------------------------------------" << std::endl;
        std::cout << std::endl;
    }

    void GeometricSpannerExperiment::addProblem( GeometricSpannerProblem prob ) {
        // add prob to problems
        problems.push_back( GeometricSpannerProblem( prob.niceName, prob.fileName ) );
    }

}
