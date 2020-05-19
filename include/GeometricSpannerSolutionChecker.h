#ifndef GEOMETRICSPANNERSOLUTIONCHECKER_H
#define GEOMETRICSPANNERSOLUTIONCHECKER_H

#include "CgalComponents.h"

namespace GeometricSpanner {

class GeometricSpannerSolutionChecker
{
    public:
        GeometricSpannerSolutionChecker( double epsilon );
        bool check( const vector<Point> &P, const list<Point> &C );
    protected:
        double epsilon;
};

}

#endif

