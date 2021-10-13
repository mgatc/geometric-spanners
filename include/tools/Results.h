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
            2,2,2,2,2
    };

    struct BoundedDegreeSpannerResult {
        DistributionType distribution;
        BoundedDegreePlaneSpannerAlgorithm algorithm;
        index_t n;
        number_t runtime;
        mixed_t degree;
        number_t degreeAvg;
        number_t stretchFactor;
        number_t lightness;
        //index_t numberOfIVs = 5;
        map<string,mixed_t> IV;

        template <class VertexIterator,class EdgeIterator>
        BoundedDegreeSpannerResult(const DistributionType distribution,
                                   const BoundedDegreePlaneSpannerAlgorithm algorithm,
                                   const number_t runtime,
                                   VertexIterator pointsBegin,
                                   VertexIterator pointsEnd,
                                   EdgeIterator edgesBegin,
                                   EdgeIterator edgesEnd,
                                   bool lite = true)
            : BoundedDegreeSpannerResult(distribution,
                                         algorithm,
                                         std::distance(pointsBegin,pointsEnd),
                                         runtime,
                                         spanners::degree(edgesBegin, edgesEnd ),
                                         spanners::degreeAvg(edgesBegin, edgesEnd ),
                                         (lite ? StretchFactorExpDijk( pointsBegin, pointsEnd, edgesBegin, edgesEnd ) :
                                            StretchFactorDijkstraReduction( pointsBegin, pointsEnd, edgesBegin, edgesEnd )),
                                         getLightness( pointsBegin, pointsEnd, edgesBegin, edgesEnd ) ) {}

        BoundedDegreeSpannerResult(const DistributionType distribution,
                                   const BoundedDegreePlaneSpannerAlgorithm algorithm,
                                   const index_t n,
                                   number_t runtime,
                                   mixed_t degree,
                                   const number_t degreeAvg,
                                   const number_t stretchFactor,
                                   const number_t lightness,
                                   bool lite = true)
                : distribution(distribution),
                  algorithm(algorithm),
                  n(n),
                  runtime(runtime),
                  degree(std::move(degree)),
                  degreeAvg(degreeAvg),
                  stretchFactor(stretchFactor),
                  lightness(lightness) {
            unsigned i=0;
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

            return ( degreePasses && stretchFactorPasses) || isBCC6;
        }

        template<class Printer>
        void setCaption( Printer& printer ) {
            string caption = string("\\textsc{")
                             + ALGORITHM_NAMES.at(algorithm)
                             + "}: "
                             + "$\\Delta = "
                             + spanners::to_string(degree);

            caption += ",\\ stretchFactor = "
                       + to_string(stretchFactor);

            caption +="$";
            printer.setCaption(caption);
        }
        //friend ostream& operator<<(ostream &os, const BoundedDegreeSpannerResult &result);
        friend ostream& operator<<(ostream &os, const BoundedDegreeSpannerResult &result) {
            os << DISTRIBUTION_NAMES.at(result.distribution) << ","
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

//    class BoundedDegreeSpannerResultSet {
//    public:
//        template<typename T>
//        using AlgorithmMap = map<BoundedDegreePlaneSpannerAlgorithm,T>;
//        template<typename T>
//        using LevelMap = map<size_t, T>;
//        template<typename T>
//        using ResultMap = AlgorithmMap<LevelMap<T>>;
//
//        typedef vector<BoundedDegreeSpannerResult> IndividualResults;
//        typedef ResultMap<IndividualResults> IndividualResultMap;
//        typedef ResultMap<BoundedDegreeSpannerResult> ReducedSamplesResultMap;
//        typedef AlgorithmMap<BoundedDegreeSpannerResult> ReducedLevelsResultMap;
//
//        IndividualResultMap m_results;
//        ReducedSamplesResultMap m_reducedSamples;
//        ReducedLevelsResultMap m_reducedLevels;
//
//        void registerResult(const BoundedDegreeSpannerResult& result ) {
//            m_results[result.algorithm][result.n].push_back(result);
//        }
////        const ReducedSamplesResultMap& getReducedResults() const {
////            return m_reducedSamples;
////        }
//        void computeStatistics(bool lite = false) {
//            for( const auto& alg : m_results ) {
//                for( auto level : alg.second ) {
//                    auto numSamples = level.second.size()-1;
//                    // Calculate averages
//                    auto canonicalEnd = level.second.end();
//                    auto canonicalBegin = next(level.second.begin());// skip first run
//                    //next(canonicalBegin);
//
//                    auto runtime = std::accumulate(canonicalBegin,
//                                                   canonicalEnd,
//                                                   0.0,
//                                                  [&]( const auto& a, const auto& b ) {
//                                                      return a + b.runtime;
//                                                  }) / numSamples;
//                    // cast the result of accumulate for integral types to floating-point
//                    number_t degree = static_cast<number_t>(std::accumulate(canonicalBegin,
//                                                                        canonicalEnd,
//                                                                        0,
//                          []( const auto& a, const auto& b ) -> index_t {
//                              return a + get<index_t>(b.degree);
//                          })) / numSamples;
//                    // already a floating point, so no need to cast
//                    number_t degreeAvg = std::accumulate(canonicalBegin,
//                                                       canonicalEnd,
//                                                       0.0,
//                                                       [&]( const auto& a, const auto& b ) {
//                                                           return a + b.degreeAvg;
//                                                       }) / numSamples;
//                    number_t lightness = std::accumulate(canonicalBegin,
//                                                       canonicalEnd,
//                                                       0.0,
//                                                       [&]( const auto& a, const auto& b ) {
//                                                           return a + b.lightness;
//                                                       }) / numSamples;
//
//                    number_t stretchFactor = lite ? canonicalBegin->stretchFactor
//                        : (std::accumulate(canonicalBegin,
//                                          canonicalEnd,
//                                          0.0,
//                                          []( const auto& a, const auto& b ) {
//                                              return a + b.stretchFactor;
//                                          }) / numSamples);
//
//                    BoundedDegreeSpannerResult result(distribution,
//                                                      alg.first,
//                                                      level.first,
//                                                      runtime,
//                                                      degree,
//                                                      degreeAvg, stretchFactor, lightness,lite);
//                    m_reducedSamples[alg.first].emplace(level.first, result);
//                }
//            }
//
//            for(auto alg : m_reducedSamples) {
//                unsigned numSamples = alg.second.size();
//                auto canonicalBegin = alg.second.begin(),
//                    canonicalEnd = alg.second.end();
//                auto runtime = std::accumulate( canonicalBegin,
//                                                canonicalEnd,
//                                                0.0,
//                                                [&]( const auto& a, const auto& b ) {
//                                                    return a + b.second.runtime;
//                                                }) / numSamples;
//                auto degree = static_cast<number_t>(std::accumulate(canonicalBegin,
//                                                                    canonicalEnd,
//                                                                    0.0,
//                                                                    []( const auto& a, const auto& b ) {
//                                                                        return a + get<number_t>(b.second.degree);
//                                                                    })) / numSamples;
//                // already a floating point, so no need to cast
//                number_t degreeAvg = std::accumulate(canonicalBegin,
//                                                     canonicalEnd,
//                                                     0.0,
//                                                     [&]( const auto& a, const auto& b ) {
//                                                         return a + b.second.degreeAvg;
//                                                     }) / numSamples;
//                number_t lightness = std::accumulate(canonicalBegin,
//                                                     canonicalEnd,
//                                                     0.0,
//                                                     [&]( const auto& a, const auto& b ) {
//                                                         return a + b.second.lightness;
//                                                     }) / numSamples;
//
//                number_t stretchFactor = lite ? canonicalBegin->second.stretchFactor
//                     : (std::accumulate(canonicalBegin,
//                                      canonicalEnd,
//                                      0.0,
//                                      []( const auto& a, const auto& b ) {
//                                          return a + b.second.stretchFactor;
//                                      }) / numSamples);
//
//                BoundedDegreeSpannerResult result(distribution,alg.first,0,runtime,degree,degreeAvg,stretchFactor,lightness,lite);
//                m_reducedLevels.emplace(alg.first, result);
//            }
//
//        }
//    };
}


//    template<class ResultSet>
//    struct AverageMetric {
//        typedef typename ResultSet::BoundedDegreeSpannerResult BoundedDegreeSpannerResult;
//        typedef map< BoundedDegreePlaneSpannerAlgorithm, map<size_t,double>> ResultMap;
//
//        ResultMap sum;
//
//        ResultMap& operator()(const ResultSet& resultSet) {
//            size_t count = 0;
//            auto firstResult = resultSet.front();
//            auto countReference = make_pair(firstResult.algorithm, firstResult.n);
//
//            for( auto result : resultSet ) {
//                auto el1 = resultSet.begin();
//                bool inserted = false;
//                tie(el1, inserted) = sum.emplace(result.algorithm,map<size_t, double>());
//
//                auto el2 = el1->second.begin();
//                tie(el2, inserted) = el1->second.emplace(result.n, 0.0);
//
//                // Skip first run
//                if(!inserted) {
//                    el2->second += result.metric;
//                }
//                count += int(make_pair(result.algorithm,result.n) == countReference);
//            }
//            for( auto& result : sum ) {
//                for( auto& level : result.second ) {
//                    level.second /= count;
//                }
//            }
//            return sum;
//        }
//    };
//
//    template<class ResultSet>
//    struct AverageRuntime {
//        typedef typename ResultSet::BoundedDegreeSpannerResult BoundedDegreeSpannerResult;
//        typedef map< BoundedDegreePlaneSpannerAlgorithm, map<size_t,size_t>> InputMap;
//        typedef map< BoundedDegreePlaneSpannerAlgorithm, map<size_t,double>> ResultMap;
//
//        InputMap sum;
//        ResultMap avg;
//
//        ResultMap& operator()(const ResultSet& resultSet) {
//            size_t count = 0;
//            auto firstResult = resultSet.front();
//            auto countReference = make_pair(firstResult.algorithm, firstResult.n);
//
//            for( auto result : resultSet ) {
//                auto el1 = resultSet.begin();
//                bool inserted = false;
//                tie(el1, inserted) = sum.emplace(result.algorithm,map<size_t, size_t>());
//
//                auto el2 = el1->second.begin();
//                tie(el2, inserted) = el1->second.emplace(result.n, 0.0);
//
//                // Skip first run
//                if(!inserted) {
//                    el2->second += result.metric;
//                }
//                count += int(make_pair(result.algorithm,result.n) == countReference);
//            }
//            for( auto result : sum ) {
//                auto el1 = avg.begin();
//                el1 = avg.emplace(result.first,map<size_t,double>()).first;
//                for( auto level : result.second ) {
//                    el1->second.emplace( double(level.second) / count );
//                }
//            }
//            return avg;
//        }
//    };
//    template<class BoundedDegreeSpannerResult>
//    class ResultSet;

#endif //SPANNERS_RESULT_H
