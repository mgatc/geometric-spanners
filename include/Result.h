//
// Created by matt on 7/20/21.
//

#ifndef GEOMETRIC_SPANNERS_RESULT_H
#define GEOMETRIC_SPANNERS_RESULT_H

#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "utilities.h"

namespace unf_spanners {
    using namespace std;

//    template<class ResultSet>
//    struct AverageMetric {
//        typedef typename ResultSet::BoundedDegreeSpannerResult BoundedDegreeSpannerResult;
//        typedef map< Algorithm, map<size_t,double>> ResultMap;
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
//        typedef map< Algorithm, map<size_t,size_t>> InputMap;
//        typedef map< Algorithm, map<size_t,double>> ResultMap;
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
    struct BoundedDegreeSpannerResult {

        Algorithm algorithm;
        index_t n;
        mixed_t runtime;
        mixed_t degree;
        number_t degreeAvg;
        number_t lightness;
        std::optional <number_t> t;
        index_t numberOfIVs = 5;

        BoundedDegreeSpannerResult(const Algorithm algorithm,
                                   const index_t n,
                                   mixed_t runtime,
                                   mixed_t degree,
                                   const number_t degreeAvg,
                                   const number_t lightness,
                                   std::optional <number_t>  t = nullopt)
                : algorithm(algorithm),
                  n(n),
                  runtime(std::move(runtime)),
                  degree(std::move(degree)),
                  degreeAvg(degreeAvg),
                  lightness(lightness),
                  t(std::move(t)),
                  numberOfIVs(numberOfIVs - int(!t)) {
        }

        template<class Printer>
        void setCaption( Printer& printer ) {
            string caption = string("\\textsc{")
                             + Names.at(algorithm)
                             + "}: "
                             + "$\\Delta = "
                             + unf_spanners::to_string(degree);

            if( t ) {
                caption += ",\\ t = "
                           + to_string(*t);
            }

            caption +="$";
            printer.setCaption(caption);
        }
        //friend ostream& operator<<(ostream &os, const BoundedDegreeSpannerResult &result);
        friend ostream& operator<<(ostream &os, const BoundedDegreeSpannerResult &result) {
            os << result.n << ","
               << Names[result.algorithm] << ","
               << result.runtime << ","
               << result.degree << ","
               << result.degreeAvg << ","
               << result.lightness << ",";

            if (result.t)
                os << *(result.t);

            os << ",\n";

            return os;
        }
    };

    class BoundedDegreeSpannerResultSet {
    public:
        template<typename T>
        using ResultMap = map< Algorithm, map<size_t, T>>;

        typedef vector<BoundedDegreeSpannerResult> IndividualResults;
        typedef ResultMap<IndividualResults> IndividualResultMap;
        typedef ResultMap<BoundedDegreeSpannerResult> ReducedResultMap;

        void registerResult(const BoundedDegreeSpannerResult& result ) {
            m_results[result.algorithm][result.n].push_back(result);
        }
        void computeStatistics() {
            for( const auto& alg : m_results ) {
                for( auto level : alg.second ) {
                    auto numSamples = static_cast<double>(level.second.size());
                    // Calculate averages
                    auto canonicalEnd = level.second.end();
                    auto canonicalBegin = level.second.begin();
                    if( numSamples > 1 )
                        next(canonicalBegin);
                    // cast the result of accumulate for integral types to floating-point
                    auto runtime = static_cast<number_t>(std::accumulate(canonicalBegin,
                                                          canonicalEnd,
                                                          0,
                                                          [&]( const auto& a, const auto& b ) {
                                                              return a + get<index_t>(b.runtime);
                                                          })) / numSamples;
                    auto degree = static_cast<number_t>(std::accumulate(canonicalBegin,
                                                         canonicalEnd,
                                                         0,
                                                         []( const auto& a, const auto& b ) {
                                                             return a + get<index_t>(b.degree);
                                                         })) / numSamples;
                    // already a floating point, so no need to cast
                    number_t degreeAvg = std::accumulate(canonicalBegin,
                                                       canonicalEnd,
                                                       0.0,
                                                       [&]( const auto& a, const auto& b ) {
                                                           return a + b.degreeAvg;
                                                       }) / numSamples;
                    number_t lightness = std::accumulate(canonicalBegin,
                                                       canonicalEnd,
                                                       0.0,
                                                       [&]( const auto& a, const auto& b ) {
                                                           return a + b.lightness;
                                                       }) / numSamples;

                    optional<number_t> t = nullopt;
                    if(!level.second.empty() && level.second.front().t) {
                        t = make_optional(std::accumulate(canonicalBegin,
                                                          canonicalEnd,
                                                          0.0,
                                                          []( const auto& a, const auto& b ) {
                                                              return a + *b.t;
                                                          }) / numSamples );
                    }
                    BoundedDegreeSpannerResult result(alg.first, level.first, runtime, degree, degreeAvg, lightness, t);
                    m_reduced[alg.first].emplace(level.first,result);
                }
            }
        }
            private:
        IndividualResultMap m_results;
        ReducedResultMap m_reduced;
    };
}

#endif //GEOMETRIC_SPANNERS_RESULT_H
