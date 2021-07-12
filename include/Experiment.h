#ifndef GEOMETRIC_SPANNERS_EXPERIMENT_H
#define GEOMETRIC_SPANNERS_EXPERIMENT_H

#include <iostream>
#include <list>
#include <optional>
#include <utility>

#include "LineGraphPrinter.h"
#include "metrics.h"
#include "utilities.h"

namespace unf_planespanners {

    template<class Container>
    void PlaneSpannerExperiment(const Container &points,
                                const Algorithm algorithm ) {
        using namespace std;

        list<pair<size_t, size_t> > spanner;

        Timer tim;

        auto pointsBegin = points.begin(),
                pointsEnd = points.end();

        switch (algorithm) {
            case Algorithm::Bgs2005:
                BGS2005(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Lw2004:
                LW2004(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bsx2009:
                BSX2009(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Kpx2010:
                KPX2010(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Kx2012:
                KX2012(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bcc2012_7:
                BCC2012<7>(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bcc2012_6:
                BCC2012<6>(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bhs2017:
                BHS2017(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bghp2010:
                BGHP2010(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Kpt2017:
                KPT2017(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bkpx2015:
                BKPX2015(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Last:
                assert(false);
        }

        size_t runtime = tim.stop();
        size_t n = points.size();
        size_t deg = degree(spanner.begin(), spanner.end());
        double t = StretchFactorUsingHeuristic(points.begin(), points.end(), spanner.begin(), spanner.end());

        Result result(algorithm, n, runtime, deg, t);
        lineGraph.registerResult(result);
        cout << result << endl;
    }

} // unf_planespanners

#endif //GEOMETRIC_SPANNERS_EXPERIMENT_H
