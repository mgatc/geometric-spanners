#include <algorithm> // min
#include <iostream>
#include <vector>

#include "Experiment.h"

using namespace std;
using namespace spanners;

int main(int argc, char *argv[]) {

    // DEFAULT ARGUMENTS IN THE EVENT COMMAND LINE INPUT IS NOT GIVEN
    const size_t runs = 5;
    const size_t n_begin = 5000;
    const size_t n_end = 10000;
    const size_t increment = 1000;

    vector<size_t> experimentParameters = {
            runs, n_begin, n_end, increment
    };

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



    return EXIT_SUCCESS;
}
