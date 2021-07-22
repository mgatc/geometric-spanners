//
// Created by matt on 7/16/21.
//

#ifndef GEOMETRIC_SPANNERS_NAMES_H
#define GEOMETRIC_SPANNERS_NAMES_H

#include <string>
#include <vector>

namespace unf_spanners {

    std::vector<std::string> ALGORITHM_NAMES = {
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

    std::vector<std::string> IV_NAMES = {
            "runtime",
            "degreeMax",
            "degreeAvg",
            "lightness",
            "stretchFactor"
    };

    std::vector<std::string> IV_NICE_NAMES = {
            "runtime",
            "degree $(\\Delta)$",
            "average degree $(\\Delta_{avg})$",
            "lightness $(\\lambda)$",
            "stretch factor $(t)$"
    };
    std::vector<std::string> PGFPLOT_NAMES = {
            "Mean " + IV_NICE_NAMES[0] + " of resultant spanners",
            "Mean maximum " + IV_NICE_NAMES[1] + " of resultant spanners",
            "Mean " + IV_NICE_NAMES[2] + " of resultant spanners",
            "Mean " + IV_NICE_NAMES[3] + " of resultant spanners",
            "Mean " + IV_NICE_NAMES[4] + " of resultant spanners"
    };

} // unf_spanners

#endif //GEOMETRIC_SPANNERS_NAMES_H
