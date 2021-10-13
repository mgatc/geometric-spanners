#ifndef GEOMETRIC_SPANNERS_DELAUNAY_H
#define GEOMETRIC_SPANNERS_DELAUNAY_H

#include "DelaunayL2.h"
#include "DelaunayLinf.h"
#include "DelaunayTD.h"

namespace spanners {

    enum DelaunayAlgorithm {
        DelaunayAlgorithmFirst=0,
        L2Sorted = DelaunayAlgorithmFirst,
        LinfSorted,
        TDSorted,
        L2Unsorted,
        LinfUnsorted,
        TDUnsorted,
        DelaunayAlgorithmLast
    };
    const vector<string> DELAUNAY_ALGORITHM_NAMES = {
            "$L_2$ Delaunay (Sorted)",
            "$L_inf$ Delaunay (Sorted)",
            "TD Delaunay (Sorted)",
            "$L_2$ Delaunay (Unsorted)",
            "$L_inf$ Delaunay (Unsorted)",
            "TD Delaunay (Unsorted)"
    };

} // spanners
#endif //GEOMETRIC_SPANNERS_DELAUNAY_H
