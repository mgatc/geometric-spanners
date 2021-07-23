//
// Created by matt on 7/16/21.
//

#ifndef GEOMETRIC_SPANNERS_NAMES_H
#define GEOMETRIC_SPANNERS_NAMES_H

#include <string>
#include <vector>

namespace unf_spanners {

    using namespace std;

    enum Algorithm {
        First=0,
        Bgs2005 = First,
        Lw2004,
        Bsx2009,
        Kpx2010,
        Kx2012,
        Bhs2017,
        Bcc2012_7,
        Bcc2012_6,
        Bghp2010,
        Kpt2017,
        Bkpx2015,
        Last
    };
    const vector<string> ALGORITHM_NAMES = {
        "BGS2005",
        "LW2004",
        "BSX2009",
        "KPX2010",
        "KX2012",
        "BHS2017",
        "BCC2012-7",
        "BCC2012-6",
        "BGHP2010",
        "BKPX2015",
        "KPT2017"
    };
    const string DEGREE_BOUND_SYMBOL = "$\\Delta_{max}$";
    const vector<string> DEGREE_BOUND_PER_ALGORITHM = {
        "27","23","17","14","11","8","7","6","6","4","4"
    };
    const string STRETCH_FACTOR_BOUND_SYMBOL = "$t_{max}$";
    const vector<string> STRETCH_FACTOR_BOUND_PER_ALGORITHM = {
        "8.27","6.44","23.6","2.92","2.86","4.41","11.7","81.7","6","157","20"
    };

    const vector<string> IV_NAMES = {
        "runtime",
        "degreeMax",
        "degreeAvg",
        "lightness",
        "stretchFactor"
    };

    const vector<string> IV_SYMBOLS = {
        "$\\rho$",
        "$\\Delta$",
        "$\\Delta_{avg}$",
        "$\\lambda$",
        "$t$"
    };
    const vector<string> IV_UNITS = {
        "s", // runtime in seconds
        "", // no unit for degree
        "", // not unit for degree
        "", // not unit for lightness
        "" // not unit for stretch factor
    };
    const vector<string> IV_NICE_NAMES = {
        "runtime",
        "degree",
        "average degree",
        "lightness",
        "stretch factor"
    };
    const vector<unsigned> IV_PRECISION = {
        2,2,2,2,2
    };
    const vector<string> PGFPLOT_NAMES = {
        "Mean " + IV_NICE_NAMES[0] + " of resultant spanners",
        "Mean maximum " + IV_NICE_NAMES[1] + " of resultant spanners",
        "Mean " + IV_NICE_NAMES[2] + " of resultant spanners",
        "Mean " + IV_NICE_NAMES[3] + " of resultant spanners",
        "Mean " + IV_NICE_NAMES[4] + " of resultant spanners"
    };

} // unf_spanners

#endif //GEOMETRIC_SPANNERS_NAMES_H
