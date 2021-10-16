#include <algorithm> // min
#include <experimental/filesystem>
#include <iostream>
#include <vector>

#include "Analysis.h"
#include "Experiment.h"
#include "Scratch.h"

using namespace std;
using namespace spanners;

int main(int argc, char *argv[]) {

    // DEFAULT ARGUMENTS IN THE EVENT COMMAND LINE INPUT IS NOT GIVEN
    const size_t runs = 5;
    const size_t n_begin = 10000;
    const size_t n_end = 100000;
    const size_t increment = 10000;

    vector<size_t> experimentParameters = {
            runs, n_begin, n_end, increment
    };

    size_t N = experimentParameters[0];

    switch(argc) {
        case 2:
            try{
                spanners::scratch(stoi(argv[1]));
            } catch(invalid_argument &ia) {
                // check extension for .xy or .csv
                string filename = argv[1];
                auto dotIndex = filename.rfind('.');

                if(dotIndex == string::npos) {
                    cout<<"INVALID FILE TYPE\n\n";
                    return EXIT_FAILURE;
                }

                string extension = filename.substr(dotIndex+1);

                if(extension == "csv") {
                    BoundedDegreePlaneSpannerAnalysis(filename);
                } else if(extension == "xy") {
                    scratch(filename);
                }

            }
            break;
        case 3:
//            ignore = system( "rm ../output/real-*");
//            try{
//                spanners::ExperimentFromConfigurationFile(stoi(argv[1]), argv[2]);
//            } catch(invalid_argument &ia) {
//                cout << "Invalid parameter '" << argv[1] << "'... exiting\n";
//                return EXIT_FAILURE;
//            }
//            break;
        case 0: // run an experiment with default args
        case 5:
        default:
            for (size_t arg = 1;
                 arg < std::min(size_t(argc), experimentParameters.size() + 1);
                 ++arg) {
                try {
                    experimentParameters[arg - 1] = stoul(argv[arg]);
                    cout << "Parameter " << (arg - 1) << " = " << experimentParameters[arg - 1] << "\n";
                }
                catch (invalid_argument &ia) {
                    cout << "Invalid parameter '" << arg << "'... exiting\n";
                }
            }
            spanners::SyntheticExperiment(experimentParameters[0],
                                          experimentParameters[1],
                                          experimentParameters[2],
                                          experimentParameters[3]);
    }

    return EXIT_SUCCESS;
}
