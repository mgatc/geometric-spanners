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

#include "libspanner/BoundedDegreePlaneSpanners.h"
#include "libspanner/types.h"

#include "Metrics.h"
#include "PointGenerators.h"
#include "Utilities.h"

namespace bdps_experiment {

    const std::string N_SYMBOL = "$n$";


    const std::vector<std::string> IV_NAMES = {
            "runtime",
            "degreeMax",
            "degreeAvg",
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

    double getDouble(const std::vector<string>& row, const size_t i) {
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

    struct BoundedDegreeSpannerAnalysisResult {
        string distribution = "";
        string algorithm = "";
        size_t n = 0;
        double runtime = 0.0;
        size_t degree = 0;
        double degreeAvg = 0.0;
        double avgDegreePerPoint = 0.0;
        double avgStretchFactor = 0.0;
        double maxStretchFactor = 0.0;
        double lightness = 0.0;
        bool lite = false;

        BoundedDegreeSpannerAnalysisResult() = default;
        BoundedDegreeSpannerAnalysisResult(const BoundedDegreeSpannerAnalysisResult& other) = default;

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
                    && abs(degreeAvg) < EPSILON
                    && abs(avgDegreePerPoint) < EPSILON
                    && abs(avgStretchFactor) < EPSILON
                    && abs(maxStretchFactor) < EPSILON
                    && abs(lightness) < EPSILON;
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
        friend std::ostream &operator<<(std::ostream &os, const BoundedDegreeSpannerAnalysisResult &result);
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
                (lhs.avgDegreePerPoint * lhs.n + rhs.avgDegreePerPoint * rhs.n)
                    / (lhs.n+int(!(lhs.lite||rhs.lite))*rhs.n), // only include in denom if not lite
                lhs.avgStretchFactor + rhs.avgStretchFactor,
                std::max(lhs.maxStretchFactor, rhs.maxStretchFactor),
                lhs.lightness + rhs.lightness
        };
    }
    BoundedDegreeSpannerAnalysisResult
    operator/(const BoundedDegreeSpannerAnalysisResult& lhs,
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
               const BoundedDegreeSpannerAnalysisResult &result) {
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
    struct BoundedDegreeSpannerResult {
        DistributionType distributionType;
        DistributionSubTypeEnum distribution;
        spanner::BoundedDegreePlaneSpannerAlgorithm algorithm;
        index_t n;
        number_t runtime;
        mixed_t degree;
        number_t degreeAvg;
        number_t stretchFactor;
        number_t lightness;
        //index_t numberOfIVs = 5;
        std::map<std::string, mixed_t> IV;
        bool lite;

        BoundedDegreeSpannerResult() = default;

        template<class VertexIterator, class EdgeIterator>
        BoundedDegreeSpannerResult(const DistributionType distributionType,
                                   const DistributionSubTypeEnum distribution,
                                   const spanner::BoundedDegreePlaneSpannerAlgorithm algorithm,
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
                                             (lite ? 0 : bdps_experiment::degree(edgesBegin, edgesEnd)),
                                             (lite ? 0.0 : bdps_experiment::degreeAvg(edgesBegin, edgesEnd)),
                                             (lite ? 0.0 :
                                                StretchFactorDijkstraReduction(pointsBegin, pointsEnd, edgesBegin, edgesEnd)),
//                                                StretchFactorExpDijk(pointsBegin, pointsEnd, edgesBegin, edgesEnd)),
                                             (lite ? 0.0 : getLightness(pointsBegin, pointsEnd, edgesBegin, edgesEnd)),
                                             lite) {}

        BoundedDegreeSpannerResult(const DistributionType distributionType,
                                   const DistributionSubTypeEnum distribution,
                                   const spanner::BoundedDegreePlaneSpannerAlgorithm algorithm,
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
                  lightness(lightness),
                  lite(lite) {
            unsigned i = 0;
            IV.emplace(IV_NAMES[i++], runtime);
            IV.emplace(IV_NAMES[i++], degree);
            IV.emplace(IV_NAMES[i++], degreeAvg);
            IV.emplace(IV_NAMES[i++], stretchFactor);
            IV.emplace(IV_NAMES[i++], lightness);
        }

        bool verify() {
            const auto degreeBound = static_cast<size_t>(stoi(spanner::bdps::DEGREE_BOUND_PER_ALGORITHM.at(algorithm)));
            const auto sfBound = static_cast<double>(stod(spanner::bdps::STRETCH_FACTOR_BOUND_PER_ALGORITHM.at(algorithm)));

            const bool degreePasses = std::get<spanner::index_t>(degree) <= degreeBound;
            const bool stretchFactorPasses = stretchFactor < sfBound || abs(stretchFactor - sfBound) < EPSILON;
            const bool isBCC6 = false;//algorithm == Bcc2012_6;

            return lite || (degreePasses && stretchFactorPasses) || isBCC6;
        }

        template<class Printer>
        void setCaption(Printer &printer) {
            std::string caption = std::string("\\textsc{")
                             + spanner::bdps::ALGORITHM_NAMES.at(algorithm)
                             + "}: "
                             + "$\\Delta = "
                             + bdps_experiment::to_string(degree);

            caption += ",\\ avgStretchFactor = "
                       + std::to_string(stretchFactor);

            caption += "$";
            printer.setCaption(caption);
        }


        //friend ostream& operator<<(ostream &os, const BoundedDegreeSpannerResult &result);
        friend std::ostream &operator<<(std::ostream &os, const BoundedDegreeSpannerResult &result) {
            const std::vector<std::string>& distributionNames = result.distributionType == DistributionType::Synthetic ?
                    SYNTHETIC_DISTRIBUTION_NAMES : REAL_POINTSET_NAMES;

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
