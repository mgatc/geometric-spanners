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

#include "libspanner/GreedySpanners.h"

#include "libspanner/measure/degree.h"
#include "libspanner/measure/stretchfactor.h"
#include "libspanner/measure/timer.h"
#include "libspanner/measure/weight.h"

#include "libspanner/points/generators.h"
#include "libspanner/types.h"
#include "libspanner/utilities.h"

#include "cpptex/cpptex.h"

#include "tools/Results.h"

namespace bdps_experiment {

    const std::string N_SYMBOL = "$n$";


    const std::vector<std::string> IV_NAMES = {
            "runtime",
            "degreeMax",
            "avgDegreePerPoint",
            "avgStretchFactor",
            "lightness"
    };

    const std::vector<std::string> IV_SYMBOLS = {
            "runtime",
            "$\\Delta_{\\mathrm{obs}}$",
            "$\\Delta_{\\mathrm{avg}}$",
            "$t_{\\mathrm{obs}}$",
            "$\\lambda$"
    };
    const std::vector<std::string> IV_UNITS = {
            "seconds", // runtime in seconds
            "", // no unit for degree
            "", // no unit for degree
            "", // no unit for lightness
            "" // no unit for stretch factor
    };
    const std::vector<std::string> IV_NICE_NAMES = {
            "Average execution time",
            "Average maximum degree",
            "Average degree per vertex",
            "Average stretch factor",
            "Average lightness"
    };
    const std::vector<unsigned> IV_PRECISION = {
            2, 2, 2, 2, 2
    };

    double getDouble(const std::vector<std::string>& row, const size_t i) {
        return std::stod(row[i]);
    }
    int getInt(const std::vector<std::string>& row, const size_t i) {
        return std::stoi(row[i]);
    }
    std::string getSpannerAlgorithm(const std::vector<std::string>& row) {
        return row[2];
    }
    std::string getDistribution(const std::vector<std::string>& row) {
        return row[0];
    }
    size_t getLevel(const std::vector<std::string>& row) {
        return getInt(row,1);
    }
    double getRuntime(const std::vector<std::string>& row) {
        return getDouble(row,3);
    }
    size_t getDegree(const std::vector<std::string>& row) {
        return getInt(row,4);
    }
    double getAvgDegreePerPoint(const std::vector<std::string>& row) {
        return getDouble(row,5);
    }
    double getStretchFactor(const std::vector<std::string>& row) {
        return getDouble(row,6);
    }
    double getLightness(const std::vector<std::string>& row) {
        return getDouble(row,7);
    }

    struct GreedySpannerAnalysisResult {
        std::string distribution = "";
        std::string algorithm = "";
        size_t n = 0;
        double runtime = 0.0;
        size_t degree = 0;
        double degreeAvg = 0.0;
        double avgDegreePerPoint = 0.0;
        double avgStretchFactor = 0.0;
        double maxStretchFactor = 0.0;
        double lightness = 0.0;
        bool lite = false;

        GreedySpannerAnalysisResult() = default;
        GreedySpannerAnalysisResult(const GreedySpannerAnalysisResult& other) = default;

        void loadRow(const std::vector<std::string>& row) {
            distribution = getDistribution(row);
            algorithm = getSpannerAlgorithm(row);
            n = getLevel(row);
            runtime = getRuntime(row);
            degree = getDegree(row);
            degreeAvg = degree;
            avgDegreePerPoint = getAvgDegreePerPoint(row);
            avgStretchFactor = getStretchFactor(row);
            maxStretchFactor = avgStretchFactor;
            lightness = getLightness(row);
            lite = degree == 0
                    && abs(degreeAvg) < spanner::EPSILON
                    && abs(avgDegreePerPoint) < spanner::EPSILON
                    && abs(avgStretchFactor) < spanner::EPSILON
                    && abs(maxStretchFactor) < spanner::EPSILON
                    && abs(lightness) < spanner::EPSILON;
        }

        template<typename T>
        T getIV(const std::string& ivName) const {
            if(ivName == "runtime")
                return runtime;
            else if(ivName == "degree")
                return degree;
            else if(ivName == "degreeAvg")
                return degreeAvg;
            else if(ivName == "avgDegreePerPoint")
                return avgDegreePerPoint;
            else if(ivName == "avgStretchFactor")
                return avgStretchFactor;
            else if(ivName == "maxStretchFactor")
                return maxStretchFactor;
            else if(ivName == "lightness")
                return lightness;

            assert(!"Invalid ivName");
            return static_cast<T>(0);
        }

        //friend ostream& operator<<(ostream &os, const BoundedDegreeSpannerResult &result);
        friend std::ostream &operator<<(std::ostream &os, const GreedySpannerAnalysisResult &result);
        friend GreedySpannerAnalysisResult operator+(const GreedySpannerAnalysisResult&,const GreedySpannerAnalysisResult&);
        friend GreedySpannerAnalysisResult operator/(const GreedySpannerAnalysisResult&,size_t);
    };


    GreedySpannerAnalysisResult
    operator+(const GreedySpannerAnalysisResult& lhs,
              const GreedySpannerAnalysisResult& rhs) {
        return {
                lhs.distribution == rhs.distribution ? lhs.distribution : "",
                lhs.algorithm == rhs.algorithm ? lhs.algorithm : "",
                lhs.n + rhs.n,
                lhs.runtime + rhs.runtime,
                std::max(lhs.degree, rhs.degree),
                lhs.degreeAvg + rhs.degreeAvg,
                (lhs.avgDegreePerPoint * lhs.n + rhs.avgDegreePerPoint * rhs.n)
                    / (lhs.n+int(!(lhs.lite||rhs.lite))*rhs.n), // only include in denom if not lite
                lhs.avgStretchFactor + rhs.avgStretchFactor,
                std::max(lhs.maxStretchFactor, rhs.maxStretchFactor),
                lhs.lightness + rhs.lightness
        };
    }
    GreedySpannerAnalysisResult
    operator/(const GreedySpannerAnalysisResult& lhs,
              const size_t divisor ) {
        auto realDivisor = static_cast<double>(divisor);
        auto effectiveDivisor = lhs.lite ? 1.0 : realDivisor;

        return {
                lhs.distribution,
                lhs.algorithm,
                static_cast<size_t>(lhs.n / realDivisor),
                lhs.runtime / realDivisor,
                lhs.degree,
                lhs.degreeAvg / effectiveDivisor,
                lhs.avgDegreePerPoint,
                lhs.avgStretchFactor / effectiveDivisor,
                lhs.maxStretchFactor,
                lhs.lightness / effectiveDivisor
        };
    }
    std::ostream&
    operator<<(std::ostream &os,
               const GreedySpannerAnalysisResult &result) {
        os << result.distribution << ","
           << result.n << ","
           << result.algorithm << ","
           << result.runtime << ","
           << result.degree << ","
           << result.degreeAvg << ","
           << result.avgDegreePerPoint << ","
           << result.avgStretchFactor << ","
           << result.maxStretchFactor << ","
           << result.lightness;
        return os;
    }

    template <typename DistributionSubTypeEnum>
    struct GreedySpannerResult {
        spanner::DistributionType distributionType;
        DistributionSubTypeEnum distribution;
        spanner::GreedySpannerAlgorithm algorithm;
        spanner::index_t n;
        spanner::number_t runtime;
        spanner::mixed_t degree;
        spanner::number_t degreeAvg;
        spanner::number_t stretchFactor;
        spanner::number_t lightness;
        //index_t numberOfIVs = 5;
        std::map<std::string, spanner::mixed_t> IV;
        bool lite;

        GreedySpannerResult() = default;

        template<class VertexIterator, class EdgeIterator>
        GreedySpannerResult(const spanner::DistributionType distributionType,
                                   const DistributionSubTypeEnum distribution,
                                   const spanner::GreedySpannerAlgorithm algorithm,
                                   const spanner::number_t runtime,
                                   VertexIterator pointsBegin,
                                   VertexIterator pointsEnd,
                                   EdgeIterator edgesBegin,
                                   EdgeIterator edgesEnd,
                                   bool lite = false)
                : GreedySpannerResult(distributionType,
                                             distribution,
                                             algorithm,
                                             std::distance(pointsBegin, pointsEnd),
                                             runtime,
                                             (lite ? 0 : spanner::degree(edgesBegin, edgesEnd)),
                                             (lite ? 0.0 : spanner::avgDegreePerPoint(edgesBegin, edgesEnd)),
                                             (lite ? 0.0 :
                                                spanner::StretchFactorDijkstraReduction(pointsBegin, pointsEnd, edgesBegin, edgesEnd)),
//                                                StretchFactorExpDijk(pointsBegin, pointsEnd, edgesBegin, edgesEnd)),
                                             (lite ? 0.0 : spanner::getLightness(pointsBegin, pointsEnd, edgesBegin, edgesEnd)),
                                             lite) {}

        GreedySpannerResult(const spanner::DistributionType distributionType,
                                   const DistributionSubTypeEnum distribution,
                                   const spanner::GreedySpannerAlgorithm algorithm,
                                   const spanner::index_t n,
                                   spanner::number_t runtime,
                                   spanner::mixed_t degree,
                                   const spanner::number_t degreeAvg,
                                   const spanner::number_t stretchFactor,
                                   const spanner::number_t lightness,
                                   bool lite = false)
                : distributionType(distributionType),
                  distribution(distribution),
                  algorithm(algorithm),
                  n(n),
                  runtime(runtime),
                  degree(std::move(degree)),
                  degreeAvg(degreeAvg),
                  stretchFactor(stretchFactor),
                  lightness(lightness),
                  lite(lite) {
            unsigned i = 0;
            IV.emplace(IV_NAMES[i++], runtime);
            IV.emplace(IV_NAMES[i++], degree);
            IV.emplace(IV_NAMES[i++], degreeAvg);
            IV.emplace(IV_NAMES[i++], stretchFactor);
            IV.emplace(IV_NAMES[i++], lightness);
        }

        bool verify(const spanner::number_t sfBound = spanner::INF) {
            return stretchFactor < sfBound || abs(stretchFactor - sfBound) < spanner::EPSILON;
        }

        template<class Printer>
        void setCaption(Printer &printer) {
            std::string caption = std::string("\\textsc{")
                             + spanner::greedy::ALGORITHM_NAMES.at(algorithm)
                             + "}: "
                             + "$\\Delta = "
                             + spanner::to_string(degree);

            caption += ",\\ avgStretchFactor = "
                       + std::to_string(stretchFactor);

            caption += "$";
            printer.setCaption(caption);
        }


        //friend ostream& operator<<(ostream &os, const BoundedDegreeSpannerResult &result);
        friend std::ostream &operator<<(std::ostream &os, const GreedySpannerResult &result) {
            const std::vector<std::string>& distributionNames = result.distributionType == spanner::DistributionType::Synthetic ?
                                                                spanner::SYNTHETIC_DISTRIBUTION_NAMES : spanner::REAL_POINTSET_NAMES;
            using spanner::operator<<;

            os << distributionNames.at(result.distribution) << ","
               << result.n << ","
               << spanner::bdps::ALGORITHM_NAMES.at(result.algorithm) << ","
               << result.runtime << ","
               << result.degree << ","
               << result.degreeAvg << ","
               << result.stretchFactor << ","
               << result.lightness << ","
               << "\n";

            return os;
        }

    };

} // bdps_experiment

#endif //SPANNERS_RESULT_H
