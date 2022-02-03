#ifndef GEOMETRIC_SPANNERS_BOUNDEDDEGREEPLANESPANNERS_H
#define GEOMETRIC_SPANNERS_BOUNDEDDEGREEPLANESPANNERS_H

#include <string>
#include <vector>

#include "Degree3.h"

#include "BCC2012.h"
#include "BGHP2010.h"
#include "BGS2005.h"
#include "BHS2018.h"
#include "BKPX2015.h"
#include "BSX2009.h"
#include "KPT2017.h"
#include "KPX2010.h"
#include "KX2012.h"
#include "LW2004.h"

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
    const string DEGREE_BOUND_SYMBOL = "$\\k$";
    const vector<string> DEGREE_BOUND_PER_ALGORITHM = {
//            "27",
//            "23", "17", "14", "11", "8", "7", "6", "6",
            "4",
//            "4",
            "3"
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
