//
// Created by matt on 7/16/21.
//

#ifndef PLANESPANNERS_NAMES_H
#define PLANESPANNERS_NAMES_H

#include <string>
#include <vector>

namespace planespanners {

    using namespace std;

    const string N_SYMBOL = "n";

    enum DistributionType {
        DistributionTypeFirst=0,
        UniformInsideSquare = DistributionTypeFirst,
        UniformInsideDisc,
        NormalInsideSquare,
        NormalClustersInsideSquare,
        ContiguousGrid,
        UniformRandomGrid,
        UniformInsideAnnulus,
  //      Real,
        DistributionTypeLast
    };
    vector<string> DISTRIBUTION_NAMES = {
        "Uniform Inside Square",
        "Uniform Inside Disc",
        "Normal Inside Square",
        "Normal Inside Square with Clusters",
        "Contiguous Grid",
        "Uniform Random Grid",
        "Uniform Inside Annulus",
        "Real-world"
    };

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
        Kpt2017,
        Bkpx2015,
        AlgorithmLast
    };
    const string ALGORITHM_SYMBOL = "Algorithm";
    const vector<string> ALGORITHM_NAMES = {
        "\\texttt{BGS2005}",
        "\\texttt{LW2004}",
        "\\texttt{BSX2009}",
        "\\texttt{KPX2010}",
        "\\texttt{KX2012}",
        "\\texttt{BHS2018}",
        "\\texttt{BCC2012-7}",
        "\\texttt{BCC2012-6}",
        "\\texttt{BGHP2010}",
        "\\texttt{BKPX2015}",
        "\\texttt{KPT2017}"
    };
    const string DEGREE_BOUND_SYMBOL = "$\\Delta_{\\mathrm{u.b.}}$";
    const vector<string> DEGREE_BOUND_PER_ALGORITHM = {
        "27","23","17","14","11","8","7","6","6","4","4"
    };
    const string STRETCH_FACTOR_BOUND_SYMBOL = "$t_{\\mathrm{u.b.}}$";
    const vector<string> STRETCH_FACTOR_BOUND_PER_ALGORITHM = {
        "8.27","6.44","23.6","2.92","2.86","4.41","11.7","81.7","6","157","20"
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
        "$\\Delta_{\\mathrm{obs.}}$",
        "$\\Delta_{\\mathrm{avg.}}$",
        "$t_{\\mathrm{obs.}}$",
        "$\\lambda$"
    };
    const vector<string> IV_UNITS = {
        "seconds", // runtime in seconds
        "", // no unit for degree
        "", // not unit for degree
        "", // not unit for lightness
        "" // not unit for stretch factor
    };
    const vector<string> IV_NICE_NAMES = {
        "Average execution time",
        "Average maximum degree",
        "Average degree per vertex",
        "Average stretch factor",
        "Average lightness"
    };
    const vector<unsigned> IV_PRECISION = {
        4,2,2,2,2
    };
    const vector<string> PGFPLOT_NAMES = {
        "Mean " + IV_NICE_NAMES[0],
        "Mean maximum " + IV_NICE_NAMES[1],
        "Mean " + IV_NICE_NAMES[2],
        "Mean " + IV_NICE_NAMES[3],
        "Mean " + IV_NICE_NAMES[4]
    };

} // planespanners

#endif //PLANESPANNERS_NAMES_H