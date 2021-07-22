#include <algorithm> // min
#include <fstream> //Reading and writing point sets.
#include <iostream>
#include <vector> // experiment paramenters

#include "Experiment.h"
#include "metrics.h"
#include "Scratch.h"
#include "utilities.h"

using namespace unf_spanners;

int main(int argc, char *argv[]) {

    const index_t runs = 100;
    const index_t n_begin = 5000;
    const index_t n_end = 10000;
    const index_t increment = 1000;

    vector<index_t> experimentParameters = {
            runs, n_begin, n_end, increment
    };

    for (size_t arg = 1;
         arg < min(size_t(argc), experimentParameters.size() + 1);
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

    if (argc == 2)
        scratch(N);
    else
        experiment(experimentParameters[0],
                   experimentParameters[1],
                   experimentParameters[2],
                   experimentParameters[3]);

    latex.display();
    pgfplots.plotResults(RESULTS.m_reduced);


    return 0;
}
