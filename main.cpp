#include <algorithm> // min
#include <fstream> //Reading and writing point sets.
#include <iostream>
#include <vector> // experiment paramenters

//Random point generation, testing.
#include <CGAL/point_generators_2.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

//#include "CGALComponents.h"
//#include "FloydWarshall.h"
//#include "GeometricSpannerPrinter.h"
////#include "GraphAlgoTV.h"
//#include "BGS2005.h"
#include "LW2004.h"
#include "BSX2009.h"
#include "KPX2010.h"
#include "BCC2012.h"
#include "BHS2017.h"
#include "KPT2017.h"
#include "BKPX2015.h"
//#include "degree_iii.h"
#include "BGHP2010.h"
#include "metrics.h"
//#include "delaunay.h"
#include "Experiment.h"
#include "Scratch.h"
#include "utilities.h"

const bool MEASURE_STRETCH_FACTOR = true;

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

    if (argc == 2)
        scratch(N);

    else
        experiment( experimentParameters[0],experimentParameters[1],experimentParameters[2],experimentParameters[3] );

    return 0;
}
