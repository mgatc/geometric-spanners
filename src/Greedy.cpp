#include "Greedy.h"
#include "CgalComponents.h"
#include "GeometricSpannerPrinter.h"

struct pair_hash {
    inline size_t operator()(const pair<int,int> &v) const {
        return v.first*31+v.second;
    }
};

inline bool isEven(int n) {
    return (n % 2 == 0);
}

inline bool isPresent(const unordered_set<pair<int,int>, pair_hash> &table, const int vertical, const int horizontal) {
    return table.find(make_pair(vertical,horizontal)) != table.end();
}

Greedy::Greedy(vector<Point> &P, list<Point> &optCenters) {

    unordered_set< pair<int,int>, pair_hash> hashTableForLatticeDiskCenters;
    const double sqrt2 = std::sqrt(2);
    const double additiveFactor = sqrt2/2;
  //  const double onePointFiveOverSqrt2 =  1.5/sqrt2;
//    const double oneOverTwoSqrt2 = 1/(2*sqrt2);
    const double sqrt2TimesOnePointFive = std::sqrt(2) * 1.5;
    const double sqrt2TimesZeroPointFive = std::sqrt(2) * 0.5;

    double verticalTimesSqrtTwo, horizontalTimesSqrt2;
    int vertical, horizontal;

    for( Point p : P) {
        vertical = floor(p.x()/sqrt2), horizontal = floor(p.y()/sqrt2);
        verticalTimesSqrtTwo  = vertical*sqrt2;
        horizontalTimesSqrt2  = horizontal*sqrt2;

        if( isPresent(hashTableForLatticeDiskCenters,vertical,horizontal))
           continue;

        if( (p.x() >= verticalTimesSqrtTwo + sqrt2TimesOnePointFive ) && isPresent(hashTableForLatticeDiskCenters,vertical+1,horizontal)
        && !(squared_distance(p,Point(sqrt2*(vertical+1)+additiveFactor,horizontalTimesSqrt2+additiveFactor))>1))
            continue;

        if( (p.x() <= verticalTimesSqrtTwo - sqrt2TimesZeroPointFive) && isPresent(hashTableForLatticeDiskCenters,vertical-1,horizontal)
           && !(squared_distance(p,Point(sqrt2*(vertical-1)+additiveFactor,horizontalTimesSqrt2+additiveFactor))>1))
            continue;

        if( (p.y() <= horizontalTimesSqrt2 - sqrt2TimesZeroPointFive) && isPresent(hashTableForLatticeDiskCenters,vertical,horizontal-1)
           && !(squared_distance(p,Point(verticalTimesSqrtTwo+additiveFactor,sqrt2*(horizontal-1)+additiveFactor))>1))
            continue;

        if( (p.y() >= horizontalTimesSqrt2 + sqrt2TimesOnePointFive) && isPresent(hashTableForLatticeDiskCenters,vertical,horizontal+1)
           && !(squared_distance(p,Point(verticalTimesSqrtTwo+additiveFactor,sqrt2*(horizontal+1)+additiveFactor))>1))
            continue;

        hashTableForLatticeDiskCenters.insert(pair<int,int>(vertical,horizontal));
        optCenters.push_back(Point(verticalTimesSqrtTwo+additiveFactor,horizontalTimesSqrt2+additiveFactor));
    }
}

