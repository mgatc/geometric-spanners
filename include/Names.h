#ifndef SPANNERS_NAMES_H
#define SPANNERS_NAMES_H

#include <string>
#include <vector>

namespace spanners {

    using namespace std;

    const string N_SYMBOL = "$n$";

    enum DistributionType {
        DistributionTypeFirst=0,
        UniformInsideSquare = DistributionTypeFirst,
        UniformInsideDisc,
        UniformOnDisc,
        NormalInsideSquare,
        NormalClustersInsideSquare,
        ContiguousGrid,
        UniformRandomGrid,
        UniformInsideAnnulus,
        DistributionTypeLast,
        Real // special case, at the end to avoid using this value in synthetic experiments
    };
    vector<string> DISTRIBUTION_NAMES = {
        "Uniform Inside Square",
        "Uniform Inside Disc",
        "Uniform On Disc",
        "Normal Inside Square",
        "Normal Inside Square with Clusters",
        "Contiguous Grid",
        "Uniform Random Grid",
        "Uniform Inside Annulus",
        "Real"
    };
//    enum ExperimentType {
//        PlaneSpanner,
//        Degree3Spanner,
//        Delaunay
//    };

//    enum DelaunayAlgorithm {
//        DelaunayAlgorithmFirst=0,
//        L2Sorted = DelaunayAlgorithmFirst,
//        LinfSorted,
//        TDSorted,
//        L2Unsorted,
//        LinfUnsorted,
//        TDUnsorted,
//        DelaunayAlgorithmLast
//    };
//    const vector<string> DELAUNAY_ALGORITHM_NAMES = {
//            "\\texttt{$L_2$ Delaunay} (Sorted)",
//            "\\texttt{$L_inf$ Delaunay} (Sorted)",
//            "\\texttt{TD Delaunay} (Sorted)",
//            "\\texttt{$L_2$ Delaunay} (Unsorted)",
//            "\\texttt{$L_inf$ Delaunay} (Unsorted)",
//            "\\texttt{TD Delaunay} (Unsorted)",
//
//    };

    enum Algorithm {
        AlgorithmFirst=0,
        Bgs2005 = AlgorithmFirst,
        Lw2004,
        Bsx2009,
        Kpx2010,
        Kx2012,
        Bhs2018,
        Bcc2012_7,
        Bcc2012_6,
        Bghp2010,
        Bkpx2015,
        Kpt2017,
//        Degree3,
        AlgorithmLast
    };
    const string ALGORITHM_SYMBOL = "Algorithm";
    const vector<string> ALGORITHM_NAMES = {
        "BGS2005",
        "LW2004",
        "BSX2009",
        "KPX2010",
        "KX2012",
        "BHS2018",
        "BCC2012-7",
        "BCC2012-6",
        "BGHP2010",
        "BKPX2015",
        "KPT2017"
    //    "Degree3"
    };
    const string DEGREE_BOUND_SYMBOL = "$\\Delta_{\\mathrm{ub}}$";
    const vector<string> DEGREE_BOUND_PER_ALGORITHM = {
        "27",
        "23","17","14","11","8","7","6","6","4","4"//,"3"
    };
    const string STRETCH_FACTOR_BOUND_SYMBOL = "$t_{\\mathrm{ub}}$";
    const vector<string> STRETCH_FACTOR_BOUND_PER_ALGORITHM = {
        "8.27",
        "6.44","23.6","2.92","2.86","4.41","11.7","81.7","6","157","20"//,"INF"
    };

    const vector<string> IV_NAMES = {
        "runtime",
        "degreeMax",
        "degreeAvg",
        "stretchFactor",
        "lightness"
    };

    const vector<string> IV_SYMBOLS = {
        "runtime",
        "$\\Delta_{\\mathrm{obs}}$",
        "$\\Delta_{\\mathrm{avg}}$",
        "$t_{\\mathrm{obs}}$",
        "$\\lambda$"
    };
    const vector<string> IV_UNITS = {
        "seconds", // runtime in seconds
        "", // no unit for degree
        "", // no unit for degree
        "", // no unit for lightness
        "" // no unit for stretch factor
    };
    const vector<string> IV_NICE_NAMES = {
        "Average execution time",
        "Average maximum degree",
        "Average degree per vertex",
        "Average stretch factor",
        "Average lightness"
    };
    const vector<unsigned> IV_PRECISION = {
        2,2,2,2,2
    };
    const vector<string> PGFPLOT_NAMES(IV_NICE_NAMES);
//    = {
//        "Mean " + IV_NICE_NAMES[0],
//        "Mean maximum " + IV_NICE_NAMES[1],
//        "Mean " + IV_NICE_NAMES[2],
//        "Mean " + IV_NICE_NAMES[3],
//        "Mean " + IV_NICE_NAMES[4]
//    };

} // spanners

#endif //SPANNERS_NAMES_H
