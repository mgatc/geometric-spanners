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

    for (size_t arg = 1;
         arg < std::min(size_t(argc), experimentParameters.size() + 1);
         ++arg) {
        try {
            experimentParameters[arg - 1] = stoul(argv[arg]);
            cout << "Parameter " << (arg - 1) << " = " << experimentParameters[arg - 1] << "\n";
        }
        catch (invalid_argument &ia) {
            cout << "Invalid experiment parameter '" << arg << "', using default value = "
                 << experimentParameters[arg - 1] << "\n";
        }
    }

    size_t N = experimentParameters.empty() ? 50 : experimentParameters[0];

    if (argc == 2) {
        ignore = system("rm ./output/scratch-*");
        spanners::scratch(N);
    } else {
        ignore = system("rm ./output/exp-*");
        spanners::experiment(experimentParameters[0],
                   experimentParameters[1],
                   experimentParameters[2],
                   experimentParameters[3]);
    }
    return 0;
}
