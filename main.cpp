#include <algorithm> // min
#include <iostream>
#include <vector>

#include "Experiment.h"
#include "Scratch.h"

const bool MEASURE_STRETCH_FACTOR = true;

using namespace std;

int main(int argc, char *argv[]) {

    const size_t runs = 100;
    const size_t n_begin = 5000;
    const size_t n_end = 10000;
    const size_t increment = 1000;

    vector<size_t> experimentParameters = {
            runs, n_begin, n_end, increment
    };


    size_t N = experimentParameters[0];

    switch(argc) {
        case 2:
            ignore = system("rm ./output/scratch-*");
            try{
                spanners::scratch(stoi(argv[1]));
            } catch(invalid_argument &ia) {
                cout << "Invalid parameter '" << argv[1] << "'... exiting\n";
                return EXIT_FAILURE;
            }
            break;
        case 3:
            ignore = system( "rm ./output/real-*");
            try{
                spanners::ExperimentFromConfigurationFile(stoi(argv[1]), argv[2]);
            } catch(invalid_argument &ia) {
                cout << "Invalid parameter '" << argv[1] << "'... exiting\n";
                return EXIT_FAILURE;
            }
            break;
        case 5:
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
            ignore = system("rm ./output/exp-*");
            spanners::SyntheticExperiment(experimentParameters[0],
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
