#include <algorithm> // min
#include <iostream>
#include <vector>

#include "Analysis.h"
#include "Experiment.h"
//#include "Scratch.h"

using namespace std;
using namespace spanners;

int main(int argc, char *argv[]) {

    // DEFAULT ARGUMENTS IN THE EVENT COMMAND LINE INPUT IS NOT GIVEN
    const string defaultFilename = "../experiment.xml";
    const size_t runs = 3;
    const size_t n_begin = 10000;
    const size_t n_end = 100000;
    const size_t increment = 10000;

    vector<size_t> experimentParameters = {
            runs, n_begin, n_end, increment
    };

    size_t n = (argc > 2 ? stoi(argv[2]) : 5);
    string filename = argc > 1 ? string(argv[1]) : defaultFilename;
    auto dotIndex = filename.rfind('.');
    string extension = filename.substr(dotIndex+1);

    const int NO_ARGS_AMOUNT = 1;

    switch(argc) {
        case 2:
            if(extension == "csv") {
                BoundedDegreePlaneSpannerAnalysis(filename);
                break;
            }
//            else if(extension == "xy") {
//                scratch(filename);
//                break;
//            }
            [[fallthrough]];
//        case NO_ARGS_AMOUNT: // run a real-world experiment with default args
        case 3:
            if(extension == "xml") {
                ExperimentFromConfigurationFile(filename,n);
            }
            break;
        case NO_ARGS_AMOUNT: // run a synthetic experiment with default args
        case 5: // run a synthetic experiment with given args
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
            SyntheticExperiment(experimentParameters[0],
                                experimentParameters[1],
                                experimentParameters[2],
                                experimentParameters[3]);
            break;
        default:
            cout<<"Invalid arguments... try again."<<endl;
            return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
