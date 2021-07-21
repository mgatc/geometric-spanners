//
// Created by matt on 7/20/21.
//

#ifndef GEOMETRIC_SPANNERS_RESULT_H
#define GEOMETRIC_SPANNERS_RESULT_H

namespace unf_spanners {

    template<class ResultSet>
    struct AverageMetric {
        typedef typename ResultSet::ResultType ResultType;
        typedef map< Algorithm, map<size_t,double>> ResultMap;

        ResultMap sum;

        ResultMap& operator()(const ResultSet& resultSet) {
            size_t count = 0;
            auto firstResult = resultSet.front();
            auto countReference = make_pair(firstResult.algorithm, firstResult.n);

            for( auto result : resultSet ) {
                auto el1 = resultSet.begin();
                bool inserted = false;
                tie(el1, inserted) = sum.emplace(result.algorithm,map<size_t, double>());

                auto el2 = el1->second.begin();
                tie(el2, inserted) = el1->second.emplace(result.n, 0.0);

                // Skip first run
                if(!inserted) {
                    el2->second += result.metric;
                }
                count += int(make_pair(result.algorithm,result.n) == countReference);
            }
            for( auto& result : sum ) {
                for( auto& level : result.second ) {
                    level.second /= count;
                }
            }
            return sum;
        }
    };

    template<class ResultSet>
    struct AverageRuntime {
        typedef typename ResultSet::ResultType ResultType;
        typedef map< Algorithm, map<size_t,size_t>> InputMap;
        typedef map< Algorithm, map<size_t,double>> ResultMap;

        InputMap sum;
        ResultMap avg;

        ResultMap& operator()(const ResultSet& resultSet) {
            size_t count = 0;
            auto firstResult = resultSet.front();
            auto countReference = make_pair(firstResult.algorithm, firstResult.n);

            for( auto result : resultSet ) {
                auto el1 = resultSet.begin();
                bool inserted = false;
                tie(el1, inserted) = sum.emplace(result.algorithm,map<size_t, size_t>());

                auto el2 = el1->second.begin();
                tie(el2, inserted) = el1->second.emplace(result.n, 0.0);

                // Skip first run
                if(!inserted) {
                    el2->second += result.metric;
                }
                count += int(make_pair(result.algorithm,result.n) == countReference);
            }
            for( auto result : sum ) {
                auto el1 = avg.begin();
                el1 = avg.emplace(result.first,map<size_t,double>()).first;
                for( auto level : result.second ) {
                    el1->second.emplace( double(level.second) / count );
                }
            }
            return avg;
        }
    };
    template<class R>
    class ResultSet;

    struct BoundedDegreeSpannerResult {

        Algorithm algorithm;
        size_t n;
        size_t runtime;
        size_t degree;
        number_t degreeAvg;
        number_t lightness;
        std::optional <number_t> t;

        BoundedDegreeSpannerResult() = default;

        BoundedDegreeSpannerResult(const Algorithm algorithm,
                                   const size_t n,
                                   const size_t runtime,
                                   const size_t degree,
                                   const number_t degreeAvg,
                                   const number_t lightness,
                                   const std::optional <number_t> t = nullopt)
                : algorithm(algorithm),
                  n(n),
                  runtime(runtime),
                  degree(degree),
                  degreeAvg(degreeAvg),
                  lightness(lightness),
                  t(t) {
        }

        void operator()(ResultSet<BoundedDegreeSpannerResult>& results) {
            AverageRuntime avg;

        }
        friend ostream& operator<<(ostream &os, const BoundedDegreeSpannerResult &result);
    };

    ostream &operator<<(ostream &os, const BoundedDegreeSpannerResult &result) {
        os << result.n << ","
           << result.algorithm << ","
           << result.runtime << ","
           << result.degree << ","
           << result.degreeAvg << ","
           << result.lightness << ",";

        if (result.t)
            os << *(result.t);

        os << ",\n";

        return os;
    }

    template<class R>
    class ResultSet : public vector<R> {
    public:
        typedef R ResultType;
        //typedef map< Algorithm, map<size_t,vector<ResultType>>> ResultMap;

        void computeStatistics() {
            for( auto alg : m_results ) {
                for( auto level : alg.second ) {
                    double averageRuntime = std::accumulate(level.second.begin(), level.second.end(), 0.0,[&]( const auto& a, const auto& b ) {
                        return a + double(b.runtime);
                    } ) / level.second.size();
                    body[0] += "("
                               + to_string(level.first)
                               + ","
                               + to_string(averageRuntime)
                               + ")\n";

                    double averageDegree = std::accumulate(level.second.begin(), level.second.end(), 0,[]( const auto& a, const auto& b ) {
                        return a + b.degree;
                    } ) / double(level.second.size());
                    body[1] += "("
                               + to_string(level.first)
                               + ","
                               + to_string(averageDegree)
                               + ")\n";

                    if(!level.second.empty()
                       && level.second.front().t) {
                        double averageT = std::accumulate(level.second.begin(), level.second.end(), 0.0,[]( const auto& a, const auto& b ) {
                            return a + *b.t;
                        } ) / double(level.second.size());
                        body[2] += "("
                                   + to_string(level.first)
                                   + ","
                                   + to_string(averageT)
                                   + ")\n";
                    }

                }
            }
        }
    private:
        ResultType m_result;
    };
}

#endif //GEOMETRIC_SPANNERS_RESULT_H
