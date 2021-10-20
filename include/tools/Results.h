#ifndef SPANNERS_RESULT_H
#define SPANNERS_RESULT_H

#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "algorithms/BoundedDegreePlaneSpanners.h"

#include "tools/Metrics.h"
#include "tools/PointGenerators.h"
#include "tools/Utilities.h"

namespace spanners {
    using namespace std;

    const string N_SYMBOL = "$n$";


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
            2, 2, 2, 2, 2
    };

    double getDouble(const vector<string>& row, const size_t i) {
        return std::stod(row[i]);
    }
    int getInt(const vector<string>& row, const size_t i) {
        return std::stoi(row[i]);
    }
    string getSpannerAlgorithm(const vector<string>& row) {
        return row[2];
    }
    string getDistribution(const vector<string>& row) {
        return row[0];
    }
    size_t getLevel(const vector<string>& row) {
        return getInt(row,1);
    }
    double getRuntime(const vector<string>& row) {
        return getDouble(row,3);
    }
    size_t getDegree(const vector<string>& row) {
        return getInt(row,4);
    }
    double getAvgDegreePerPoint(const vector<string>& row) {
        return getDouble(row,5);
    }
    double getStretchFactor(const vector<string>& row) {
        return getDouble(row,6);
    }
    double getLightness(const vector<string>& row) {
        return getDouble(row,7);
    }

    struct BoundedDegreeSpannerAnalysisResult {
        string distribution = "";
        string algorithm = "";
        size_t n = 0;
        double runtime = 0.0;
        size_t degree = 0;
        double degreeAvg = 0.0;
        double avgDegreePerPoint = 0.0;
        double stretchFactor = 0.0;
        double lightness = 0.0;

        BoundedDegreeSpannerAnalysisResult() = default;
        BoundedDegreeSpannerAnalysisResult(const BoundedDegreeSpannerAnalysisResult& other) = default;

        void loadRow(const vector<string>& row) {
            distribution = getDistribution(row);
            algorithm = getSpannerAlgorithm(row);
            n = getLevel(row);
            runtime = getRuntime(row);
            degree = getDegree(row);
            degreeAvg = getDegree(row);
            avgDegreePerPoint = getAvgDegreePerPoint(row);
            stretchFactor = getStretchFactor(row);
            lightness = getLightness(row);
        }

        template<typename T>
        T getIV(const string& ivName) const {
            if(ivName == "runtime")
                return runtime;
            else if(ivName == "degree")
                return degree;
            else if(ivName == "degreeAvg")
                return degreeAvg;
            else if(ivName == "avgDegreePerPoint")
                return avgDegreePerPoint;
            else if(ivName == "stretchFactor")
                return stretchFactor;
            else if(ivName == "lightness")
                return lightness;

            assert(!"Invalid ivName");
            return static_cast<T>(0);
        }

        //friend ostream& operator<<(ostream &os, const BoundedDegreeSpannerResult &result);
        friend ostream &operator<<(ostream &os, const BoundedDegreeSpannerAnalysisResult &result);
        friend BoundedDegreeSpannerAnalysisResult operator+(const BoundedDegreeSpannerAnalysisResult&,const BoundedDegreeSpannerAnalysisResult&);
        friend BoundedDegreeSpannerAnalysisResult operator/(const BoundedDegreeSpannerAnalysisResult&,size_t);
    };


    BoundedDegreeSpannerAnalysisResult
    operator+(const BoundedDegreeSpannerAnalysisResult& lhs,
              const BoundedDegreeSpannerAnalysisResult& rhs) {
        return {
                lhs.distribution == rhs.distribution ? lhs.distribution : "",
                lhs.algorithm == rhs.algorithm ? lhs.algorithm : "",
                lhs.n + rhs.n,
                lhs.runtime + rhs.runtime,
                std::max(lhs.degree, rhs.degree),
                lhs.degreeAvg + rhs.degreeAvg,
                (lhs.avgDegreePerPoint * lhs.n + rhs.avgDegreePerPoint * rhs.n)/(lhs.n+rhs.n),
                lhs.stretchFactor + rhs.stretchFactor,
                lhs.lightness + rhs.lightness
        };
    }
    BoundedDegreeSpannerAnalysisResult
    operator/(const BoundedDegreeSpannerAnalysisResult& lhs,
              const size_t divisor ) {
        auto realDivisor = static_cast<double>(divisor);
        return {
                lhs.distribution,
                lhs.algorithm,
                static_cast<size_t>(lhs.n / realDivisor),
                lhs.runtime / realDivisor,
                lhs.degree,
                lhs.degreeAvg / realDivisor,
                lhs.avgDegreePerPoint,
                lhs.stretchFactor / realDivisor,
                lhs.lightness / realDivisor
        };
    }
    ostream&
    operator<<(ostream &os,
               const BoundedDegreeSpannerAnalysisResult &result) {
        os << result.distribution << ","
           << result.n << ","
           << result.algorithm << ","
           << result.runtime << ","
           << result.degree << ","
           << result.degreeAvg << ","
           << result.avgDegreePerPoint << ","
           << result.stretchFactor << ","
           << result.lightness;
        return os;
    }

    template <typename DistributionSubTypeEnum>
    struct BoundedDegreeSpannerResult {
        DistributionType distributionType;
        DistributionSubTypeEnum distribution;
        BoundedDegreePlaneSpannerAlgorithm algorithm;
        index_t n;
        number_t runtime;
        mixed_t degree;
        number_t degreeAvg;
        number_t stretchFactor;
        number_t lightness;
        //index_t numberOfIVs = 5;
        map<string, mixed_t> IV;

        BoundedDegreeSpannerResult() = default;

        template<class VertexIterator, class EdgeIterator>
        BoundedDegreeSpannerResult(const DistributionType distributionType,
                                   const DistributionSubTypeEnum distribution,
                                   const BoundedDegreePlaneSpannerAlgorithm algorithm,
                                   const number_t runtime,
                                   VertexIterator pointsBegin,
                                   VertexIterator pointsEnd,
                                   EdgeIterator edgesBegin,
                                   EdgeIterator edgesEnd,
                                   bool lite = false)
                : BoundedDegreeSpannerResult(distributionType,
                                             distribution,
                                             algorithm,
                                             std::distance(pointsBegin, pointsEnd),
                                             runtime,
                                             (lite ? 0 : spanners::degree(edgesBegin, edgesEnd)),
                                             (lite ? 0.0 : spanners::degreeAvg(edgesBegin, edgesEnd)),
                                             (lite ? 0.0 : StretchFactorExpDijk(pointsBegin, pointsEnd, edgesBegin, edgesEnd)),
                                             (lite ? 0.0 : getLightness(pointsBegin, pointsEnd, edgesBegin, edgesEnd)) ) {}

        BoundedDegreeSpannerResult(const DistributionType distributionType,
                                   const DistributionSubTypeEnum distribution,
                                   const BoundedDegreePlaneSpannerAlgorithm algorithm,
                                   const index_t n,
                                   number_t runtime,
                                   mixed_t degree,
                                   const number_t degreeAvg,
                                   const number_t stretchFactor,
                                   const number_t lightness,
                                   bool lite = false)
                : distributionType(distributionType),
                  distribution(distribution),
                  algorithm(algorithm),
                  n(n),
                  runtime(runtime),
                  degree(std::move(degree)),
                  degreeAvg(degreeAvg),
                  stretchFactor(stretchFactor),
                  lightness(lightness) {
            unsigned i = 0;
            IV.emplace(IV_NAMES[i++], runtime);
            IV.emplace(IV_NAMES[i++], degree);
            IV.emplace(IV_NAMES[i++], degreeAvg);
            IV.emplace(IV_NAMES[i++], stretchFactor);
            IV.emplace(IV_NAMES[i++], lightness);
        }

        bool verify() {
            const auto degreeBound = static_cast<size_t>(stoi(DEGREE_BOUND_PER_ALGORITHM.at(algorithm)));
            const auto sfBound = static_cast<double>(stod(STRETCH_FACTOR_BOUND_PER_ALGORITHM.at(algorithm)));

            const bool degreePasses = get<index_t>(degree) <= degreeBound;
            const bool stretchFactorPasses = stretchFactor < sfBound || abs(stretchFactor - sfBound) < EPSILON;
            const bool isBCC6 = algorithm == Bcc2012_6;

            return (degreePasses && stretchFactorPasses) || isBCC6;
        }

        template<class Printer>
        void setCaption(Printer &printer) {
            string caption = string("\\textsc{")
                             + ALGORITHM_NAMES.at(algorithm)
                             + "}: "
                             + "$\\Delta = "
                             + spanners::to_string(degree);

            caption += ",\\ stretchFactor = "
                       + to_string(stretchFactor);

            caption += "$";
            printer.setCaption(caption);
        }


        //friend ostream& operator<<(ostream &os, const BoundedDegreeSpannerResult &result);
        friend ostream &operator<<(ostream &os, const BoundedDegreeSpannerResult &result) {
            const vector<string>& distributionNames = result.distributionType == DistributionType::Synthetic ?
                    SYNTHETIC_DISTRIBUTION_NAMES : REAL_POINTSET_NAMES;

            os << distributionNames.at(result.distribution) << ","
               << result.n << ","
               << ALGORITHM_NAMES.at(result.algorithm) << ","
               << result.runtime << ","
               << result.degree << ","
               << result.degreeAvg << ","
               << result.stretchFactor << ","
               << result.lightness << ","
               << "\n";

            return os;
        }

    };



} // spanners

#endif //SPANNERS_RESULT_H
