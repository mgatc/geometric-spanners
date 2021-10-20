#ifndef GEOMETRIC_SPANNERS_ANALYSIS_H
#define GEOMETRIC_SPANNERS_ANALYSIS_H

#include <fstream>
#include <string>
#include <vector>

#include "algorithms/BoundedDegreePlaneSpanners.h"

#include "printers/LatexPrinter.h"
#include "printers/PgfplotPrinter.h"
#include "printers/TablePrinter.h"

#include "tools/Results.h"
#include "tools/Utilities.h"

namespace spanners {

namespace analysis {

    double X_PLOT_SCALE = 1000.0;
    double Y_PLOT_SCALE = 1.0;

    using namespace std;

    typedef size_t level_t;
    typedef string spanner_t;
    typedef string distribution_t;
    typedef BoundedDegreeSpannerAnalysisResult result_t;

    typedef vector<result_t> ResultContainer;

    typedef map<level_t, ResultContainer> LevelResultMap;
    typedef map<level_t, result_t> ReducedLevelResultMap;

    typedef map<spanner_t, result_t> SpannerResultMap;
    typedef map<spanner_t, LevelResultMap> SpannerLevelResultMap;
    typedef map<distribution_t, SpannerLevelResultMap> DistributionSpannerLevelResultMap;

    typedef map<spanner_t, ReducedLevelResultMap> SpannerReducedLevelResultMap;
    typedef map<distribution_t, SpannerReducedLevelResultMap> DistributionSpannerReducedLevelResultMap;

    typedef map<distribution_t, LevelResultMap> DistributionLevelResultMap;
    typedef map<distribution_t, ReducedLevelResultMap> DistributionReducedLevelResultMap;
    typedef map<distribution_t, result_t> DistributionSummaryMap;

    const std::array<string,6> ANALYSIS_IVs = {
        "runtime","degree","degreeAvg","avgDegreePerPoint","stretchFactor","lightness"
    };

    vector<vector<string>>
    readFileIntoVector(const string &filename) {
        ifstream expIn;
        expIn.open(filename, ios_base::in);
        if (!expIn.is_open()) assert(!"Error opening file");

        vector<vector<string>> results;
        auto row = getNextLineAndSplitIntoTokens(expIn);
        while (!(row = getNextLineAndSplitIntoTokens(expIn)).empty()) {
            results.push_back(row);
        }
        expIn.close();

        return results;
    }

    DistributionSpannerLevelResultMap
    getDistributionSpannerLevelResults(const vector<vector<string>> &rawResults) {
        DistributionSpannerLevelResultMap results;
        auto headerRow = rawResults.front();
        bool first = true;
        for (auto row : rawResults) {
            if(first) {
                first = false;
                continue;
            }
            result_t result;
            result.loadRow(row);

            results[result.distribution]
                   [result.algorithm]
                   [result.n].push_back(result);
        }
        return results;
    }

    DistributionSpannerReducedLevelResultMap
    calculateDistributionSpannerSummary(const DistributionSpannerLevelResultMap &results) {
        auto samplesPerLevelSpanner = results.begin()->second.begin()->second.size();
        DistributionSpannerReducedLevelResultMap reduced;
        for (const auto &distribution : results) {
            auto distName = distribution.first;
            const auto &distributionResults = distribution.second;
            ReducedLevelResultMap levelSummary;
            for (const auto &spanner : distributionResults) {
                for (auto level : spanner.second) {
                    auto sum = std::accumulate(level.second.begin(), level.second.end(),
                                               result_t());
                    auto mean = sum / level.second.size();

                    levelSummary[level.first] = mean;
                }
                reduced[distName][spanner.first] = levelSummary;
            }
        }
        return reduced;
    }

    void
    plot(const DistributionSpannerReducedLevelResultMap &summary,
         LatexPrinter &document) {
        bool first = true;
        for (const auto &distribution : summary) {
            for(const auto &iv : ANALYSIS_IVs) {
                string caption = distribution.first;

                string plotName("document-plot-");
                plotName += removeSpaces(caption);
                plotName += "-";
                plotName += iv;

                PgfplotPrinter plot("./", plotName);
                plot.setCaption(caption);
                plot.plotAxis(iv, distribution.second, X_PLOT_SCALE, Y_PLOT_SCALE,first);
                document.addToDocument(plot);

                if(first) {
                    first = false;
                    string legendRefText = plot.getLegend();
                    document.addRawText(legendRefText);
                    document.addRawText("\n\n");
                }
            }
            document.clearpage();
        }
    }

    void
    tabulate(DistributionSpannerReducedLevelResultMap& spannerSummary,
             LatexPrinter& document) {
        for(auto distribution : spannerSummary) {
            for(auto iv : ANALYSIS_IVs) {
                vector<string> headers;
                vector<string> levels;
                const auto& frontRow = distribution.second.begin()->second;
                std::transform(frontRow.begin(), frontRow.end(), back_inserter(levels),
                               [](const auto& elem){
                                   return std::to_string(elem.first);
                               });

                vector<vector<string>> ivValues;
                for(auto spanner : distribution.second) {
                    headers.push_back(spanner.first);
                    ivValues.push_back({});
                    std::transform(spanner.second.begin(), spanner.second.end(), back_inserter(ivValues.back()),
                                   [iv](const auto& elem){
                                       return std::to_string(elem.second.template getIV<double>(iv));
                                   });
                }

                string caption(distribution.first + " x " + iv);
                string tableName("document-table-");
                tableName += removeSpaces(caption);
                tableName += "-" + iv;
                TablePrinter table("./", tableName);

                // add columns
                table.addColumn(N_SYMBOL,levels,0);

                for(size_t i=0; i<headers.size();++i){
                    table.addColumn( headers.at(i), ivValues.at(i),i+1);
                }

                table.tabulate();

                document.addRawText("\\subsection{" + caption + "}\n\n");
                document.addToDocument(table);
            }
        }
    }

    SpannerLevelResultMap
    getSpannerLevelResults(const vector<vector<string>>& rawResults) {
        SpannerLevelResultMap results;
        auto headerRow = rawResults.front();

        bool first = true;

        for(auto row : rawResults) {
            if(first) {
                first = false;
                continue;
            }
            result_t result;
            result.loadRow(row);

            results[result.algorithm]
                   [result.n].push_back(result);
        }
        return results;
    }

    SpannerReducedLevelResultMap
    calculateSpannerLevelSummary(const SpannerLevelResultMap& results) {
        SpannerReducedLevelResultMap summary;

        for(const auto& spanner : results) {
            for(const auto& level : spanner.second) {
                result_t sum = std::accumulate(level.second.begin(),level.second.end(),result_t());
                result_t avg = sum / level.second.size();

                summary[spanner.first]
                [level.first] = avg;
            }
        }
        return summary;
    }

    SpannerResultMap
    calculateSpannerSummary(const SpannerReducedLevelResultMap& results) {
        SpannerResultMap summary;
        for(auto spanner : results) {
            auto sum = accumulate(spanner.second.begin(),
                                  spanner.second.end(),
                                  result_t(),
                [](const result_t& sum, const auto& addend) {
                    return sum + addend.second;
                });
            auto avg = sum / spanner.second.size();
            summary[spanner.first] = avg;
        }
        return summary;
    }

    void
    plot(const SpannerReducedLevelResultMap& summary,
         LatexPrinter& document ) {
        for(const auto& iv : ANALYSIS_IVs) {
            string plotName("document-plot-summary-");
            plotName += iv;
            string caption = "Summary";

            PgfplotPrinter plot("./", plotName);
            plot.setCaption(caption);
            plot.plotAxis(iv, summary, X_PLOT_SCALE, Y_PLOT_SCALE, false);

            document.addToDocument(plot);
        }
    }

    void
    tabulate(const SpannerReducedLevelResultMap& summary,
             LatexPrinter& document) {
        for(auto iv : ANALYSIS_IVs) {

            string caption(iv);
            string tableName("document-table-summary");
            tableName += removeSpaces(caption);
            tableName += "-" + iv;
            TablePrinter table("./", tableName);

            // add columns
            vector<string> headers;
            vector<vector<string>> ivValues;

            headers.push_back(N_SYMBOL);
            ivValues.push_back({});
            for(auto level : summary.begin()->second) {
                ivValues.front().push_back(std::to_string(level.first));
            }
            for(auto spanner : summary) {
                headers.push_back("\\texttt{" + spanner.first + "}");
                ivValues.push_back({});
                for(auto level : spanner.second) {
                    ivValues.back().push_back(std::to_string(level.second.getIV<double>(iv)));
                }
            }
            for(size_t i=0; i<headers.size();++i){
                table.addColumn( headers.at(i), ivValues.at(i),i);
            }
            table.tabulate();

            document.addRawText("\\subsection{" + caption + "}\n\n");
            document.addToDocument(table);
        }
    }

    void
    tabulate(const SpannerResultMap& summary,
             LatexPrinter& document) {
        vector<string> header;
        header.push_back(ALGORITHM_SYMBOL);
        copy(next(ANALYSIS_IVs.begin()),ANALYSIS_IVs.end(),back_inserter(header));

        vector<vector<string>> body;
        body.push_back({});
        transform(ALGORITHM_NAMES.begin(),ALGORITHM_NAMES.end(),back_inserter(body.back()),
            [](const string& name) {
                return "\\texttt{" + name + "}";
            });

        for(auto iv : ANALYSIS_IVs) {
            if(iv == "runtime") continue;

            body.push_back({});
            for(auto spanner : ALGORITHM_NAMES) {
                body.back().push_back(std::to_string(summary.at(spanner).template getIV<double>(iv)));
            }
        }
        string caption = "Spanner Summary";
        string tableName("document-table-summary");
        TablePrinter table("./", tableName);

        for(size_t i=0; i<body.size(); i++) {
            table.addColumn(header[i], body[i],i);
        }
        table.tabulate();

        document.addRawText("\\subsection{" + caption + "}\n\n");
        document.addToDocument(table);
    }


} // analysis



void BoundedDegreePlaneSpannerAnalysis(const string filename) {
    using namespace spanners::analysis;

    LatexPrinter document("./", "document");
    vector<vector<string>> raw = readFileIntoVector(filename);


    cout<< "Finding summaries for each spanner algorithm within each distribution..."<<endl;
    DistributionSpannerLevelResultMap distributionSpannerLevelResults = getDistributionSpannerLevelResults(raw);
    DistributionSpannerReducedLevelResultMap distributionSpannerSummary = calculateDistributionSpannerSummary(distributionSpannerLevelResults);
    plot(distributionSpannerSummary, document);
    tabulate(distributionSpannerSummary, document);

    cout<<"Finding average per level for all distributions..."<<endl;
    document.addRawText("\\section{Overall Summary}\n\n");
    SpannerLevelResultMap spannerLevelResults = getSpannerLevelResults(raw);

//    auto numLevels = spannerLevelResults.size();
//    auto samplesPerLevel = spannerLevelResults.begin()->second.size();
//    auto numSamples = numLevels * samplesPerLevel;


    SpannerReducedLevelResultMap spannerLevelSummary = calculateSpannerLevelSummary(spannerLevelResults);
    plot(spannerLevelSummary, document);
    tabulate(spannerLevelSummary, document);

    SpannerResultMap spannerSummary = calculateSpannerSummary(spannerLevelSummary);
    tabulate(spannerSummary,document);


//    tabulateDistributionSummaries(distributionSpannerSummary, document);
    //plotDistributionSummaries();

//    document.clearpage();

    //tabulateOverallSummary(spannerLevelResults, document);

    document.display();
}


} // spanners



//    void
//    tabulateDistributionSummaries(const DistributionSpannerReducedLevelResultMap& summary,
//                             LatexPrinter& document) {
//        vector<string> headers;
//        headers.push_back("Distribution");
//        std::copy(ALGORITHM_NAMES.begin(),ALGORITHM_NAMES.end(),back_inserter(headers));
////        std::transform(summary.begin(), summary.end(), back_inserter(headers),
////            [](const auto& elem){
////                return elem.first;
////            });
//        for(auto iv : ANALYSIS_IVs) {
//
//            vector<vector<string>> ivValues;
//            ivValues.push_back(SYNTHETIC_DISTRIBUTION_NAMES);
//
//            for(auto spanner : ALGORITHM_NAMES) {
//                ivValues.push_back({});
//                for(auto distribution : summary) {
//                    ivValues.back().push_back(distribution.second.at(spanner));
//                }
//            }
//        }
//
//        string caption(iv);
//        string tableName("document-table-");
//        tableName += removeSpaces(caption);
//        TablePrinter table("./", tableName);
//
//        // add columns
//
//        for(size_t i=0; i<headers.size();++i){
//            table.addColumn( headers.at(i), ivValues.at(i),i);
//        }
//
//        table.tabulate();
//
//        document.addRawText("\\subsection{" + caption + "}\n\n");
//        document.addToDocument(table);
//        }
//    }

//    void
//    tabulateOverallSummary(const DistributionLevelResultMap& results,
//                           LatexPrinter& document) {
//        vector<string> headers = { "Algorithm", "Max degree",
//                                   "Avg degree", "Avg degree/point",
//                                   "Max s.f.", "Avg s.f." };
//        map<string,vector<string>> body;
//
//        for(const auto& spanner : results) {
//            body[headers[0]].push_back(spanner.first);
//            size_t maxDegree = std::accumulate(spanner.second.begin(),spanner.second.end(),SIZE_T_MAX,
//                [](const size_t maxFound, const LevelResultMap& levels){
//                    size_t maxInLevel = std::accumulate(levels.begin(), levels.end(), SIZE_T_MAX,
//                        [](const size_t maxFound, const ResultContainer& results ){
//                            size_t maxInResults = std::accumulate(results.begin(),results.end(),SIZE_T_MAX,
//                                [](const size_t maxFound, const BoundedDegreeSpannerResult& result) {
//                                    return std::max(maxFound,static_cast<size_t>(rint(std::get<double>(result.degree))));
//                                });
//                            return std::max(maxFound, maxInResults);
//                        });
//                    return std::max(maxFound, maxInLevel);
//                });
//            body[headers[1]].push_back(std::to_string(maxDegree));
//
//            // add max degreeMax, avg degreeMax, avg degreeAvg
//        }
//
//        TablePrinter table("./","document-table-overall");
//
//        for(auto header : headers) {
//            if(contains(body,header))
//                table.addColumn(header, body[header] );
//        }
//        table.tabulate();
//        document.addToDocument(table);
//    }


//    DistributionSpannerResultMap
//    calculateDistributionSummary(const SpannerLevelResultMap& results) {
//        SpannerReducedLevelResultMap summary;
//
//        for(const auto& spanner : results) {
//            for(const auto& level : spanner.second) {
//                result_t sum = std::accumulate(level.second.begin(),level.second.end(),result_t());
//                result_t avg = sum / level.second.size();
//
//                summary[spanner.first]
//                [level.first] = avg;
//            }
//        }
//        return summary;
//    }

#endif //GEOMETRIC_SPANNERS_ANALYSIS_H
