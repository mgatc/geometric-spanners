#ifndef GEOMETRIC_SPANNERS_ANALYSIS_H
#define GEOMETRIC_SPANNERS_ANALYSIS_H

#include <fstream>
#include <string>
#include <vector>

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
    calculateDistributionSpannerSummaries(const DistributionSpannerLevelResultMap &results) {
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
    plotDistributionSpannerReducedLevels(const DistributionSpannerReducedLevelResultMap &summary,
                                              LatexPrinter &document) {
        for (const auto &distribution : summary) {
            for(const auto &iv : ANALYSIS_IVs) {
                string caption = distribution.first;

                string plotName("plot-");
                plotName += removeSpaces(caption);
                plotName += "-";
                plotName += iv;

                PgfplotPrinter plot("./", plotName);
                plot.setCaption(caption);
                plot.plotAxis(iv, distribution.second, X_PLOT_SCALE, Y_PLOT_SCALE);
                document.addToDocument(plot);
            }
            document.clearpage();
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
    calculateSpannerSummary(const SpannerLevelResultMap& results) {
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

    void
    plotSpannerSummary(const SpannerReducedLevelResultMap& summary,
                LatexPrinter& document ) {
        for(const auto& iv : ANALYSIS_IVs) {
            string plotName("document-plot-summary-");
            plotName += iv;
            string caption = "Summary";

            PgfplotPrinter plot("./", plotName);
            plot.setCaption(caption);
            plot.plotAxis(iv, summary, X_PLOT_SCALE, Y_PLOT_SCALE);

            document.addToDocument(plot);
        }
    }

//    void
//    tabulateDistributionSummaries(const DistributionLevelResultMap& spannerLevel,
//                                  const DistributionReducedLevelResultMap& summary,
//                                  LatexPrinter& document) {
//        for(const auto& distribution : summary) {
//            string tableName("document-table-spanner-summary-");
//            tableName += removeSpaces(distribution.first);
//
//            TablePrinter table("./", tableName);
//
//            int i = 0;
//
//            vector<string> levels;
//            double nScale = 1.0;
//            for( auto level : distribution.second) {
//                levels.push_back(std::to_string(static_cast<int>(rint((double)level.first / nScale))));
//            }
//            table.addColumn("$n$}{", levels, i++);
//
//            map<string,vector<string>> ColumnMap;
//            for(auto name : IV_NAMES) {
//                for(auto level : distribution.second) {
//                    ColumnMap[name].push_back(spanners::to_string(level.second.IV[name]));
//                }
//                table.addColumn(name, ColumnMap[name], i++);
//            }
//
//            table.tabulate();
//            document.addToDocument(table);
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

} // analysis



void BoundedDegreePlaneSpannerAnalysis(const string filename) {
    using namespace spanners::analysis;

    LatexPrinter document("./", "document");
    vector<vector<string>> raw = readFileIntoVector(filename);


    cout<< "Finding summaries for each spanner algorithm within each distribution..."<<endl;
    DistributionSpannerLevelResultMap spannerResults = getDistributionSpannerLevelResults(raw);
    DistributionSpannerReducedLevelResultMap spannerSummary = calculateDistributionSpannerSummaries(spannerResults);
    plotDistributionSpannerReducedLevels(spannerSummary,document);
    //tabulateDistributionSpannerReducedLevels(spannerResults,spannerSummary,document);

    cout<<"Finding average per level for all distributions..."<<endl;
    document.addRawText("\\section{Overall Summary}\n\n");
    SpannerLevelResultMap results = getSpannerLevelResults(raw);

    auto numLevels = results.size();
    auto samplesPerLevel = results.begin()->second.size();
    auto numSamples = numLevels * samplesPerLevel;


    SpannerReducedLevelResultMap summary = calculateSpannerSummary(results);
    plotSpannerSummary(summary, document);
    //tabulateDistributionSummaries(results, summary, document);

//    document.clearpage();

    //tabulateOverallSummary(results, document);

    document.display();
}


} // spanners

#endif //GEOMETRIC_SPANNERS_ANALYSIS_H
