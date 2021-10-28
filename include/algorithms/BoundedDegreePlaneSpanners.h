#ifndef GEOMETRIC_SPANNERS_BOUNDEDDEGREEPLANESPANNERS_H
#define GEOMETRIC_SPANNERS_BOUNDEDDEGREEPLANESPANNERS_H

#include <string>
#include <vector>

#include "algorithms/Degree3.h"

#include "algorithms/BCC2012.h"
#include "algorithms/BGHP2010.h"
#include "algorithms/BGS2005.h"
#include "algorithms/BHS2018.h"
#include "algorithms/BKPX2015.h"
#include "algorithms/BSX2009.h"
#include "algorithms/KPT2017.h"
#include "algorithms/KPX2010.h"
#include "algorithms/KX2012.h"
#include "algorithms/LW2004.h"

namespace spanners {

    using namespace std;

    enum BoundedDegreePlaneSpannerAlgorithm {
        AlgorithmFirst = 0,
//        Bgs2005 = AlgorithmFirst,
//        Lw2004,
//        Bsx2009,
//        Kpx2010,
//        Kx2012,
//        Bhs2018,
//        Bcc2012_7,
//        Bcc2012_6,
//        Bghp2010,
        Bkpx2015 = AlgorithmFirst,
//        Kpt2017,
        Degree3,
        AlgorithmLast
    };
    const string ALGORITHM_SYMBOL = "Algorithm";
    const vector<string> ALGORITHM_NAMES = {
//            "BGS2005",
//            "LW2004",
//            "BSX2009",
//            "KPX2010",
//            "KX2012",
//            "BHS2018",
//            "BCC2012-7",
//            "BCC2012-6",
//            "BGHP2010",
            "BKPX2015",
//            "KPT2017",
            "Degree3"
    };
    const string DEGREE_BOUND_SYMBOL = "$\\Delta_{\\mathrm{ub}}$";
    const vector<string> DEGREE_BOUND_PER_ALGORITHM = {
//            "27",
//            "23", "17", "14", "11", "8", "7", "6", "6",
//            "4",
            "4","3"
    };
    const string STRETCH_FACTOR_BOUND_SYMBOL = "$t_{\\mathrm{ub}}$";
    const vector<string> STRETCH_FACTOR_BOUND_PER_ALGORITHM = {
//            "8.27",
//            "6.44", "23.6", "2.92", "2.86", "4.41", "11.7", "81.7", "6",
            "157" ,
//            "20",
            "INF"
    };

} // spanners

#endif //GEOMETRIC_SPANNERS_BOUNDEDDEGREEPLANESPANNERS_H
