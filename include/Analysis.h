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

    double X_PLOT_SCALE = 1000000.0;
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

    const std::array<string,7> ANALYSIS_IVs = {
        "runtime","degree","degreeAvg","avgDegreePerPoint",
        "avgStretchFactor","maxStretchFactor","lightness"
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
        bool isFirst = true;
        for (auto row : rawResults) {
            if(isFirst) {
                isFirst = false;
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

        bool isFirst = true;
        for(const auto& distributionName : SYNTHETIC_DISTRIBUTION_NAMES ) {
            if(summary.find(distributionName)!=summary.end()){
                const auto &distribution = summary.at(distributionName);
                set<string> tdPlots = {"BGHP2010", "KPT2017"};
                set<string> linfPlots = {"BKPX2015"};
                set<string> l2Plots;
                copy_if(ALGORITHM_NAMES.begin(),ALGORITHM_NAMES.end(),
                        inserter(l2Plots,l2Plots.end()),
                        [&](const string& name){
                            return !(contains(tdPlots,name) || contains(linfPlots,name));
                        });

                set<string> nontdPlots;
                copy_if(ALGORITHM_NAMES.begin(), ALGORITHM_NAMES.end(),
                        inserter(nontdPlots, nontdPlots.end()),
                        [&tdPlots](const string &name) { // return true if name is not in tdPlots
                            return !contains(tdPlots, name);
                        });
                // split distribution.second into containers for td and non-td, l2 and nonl2
                //SpannerReducedLevelResultMap tdResults, nontdResults, linfResults, l2Results;

                // TD plots
                SpannerReducedLevelResultMap tdResults;
                copy_if(distribution.begin(), distribution.end(),
                        inserter(tdResults, tdResults.end()),
                        [&tdPlots](const auto &elem) {
                            return contains(tdPlots, elem.first);
                        });

                if( !tdResults.empty()) {
                    string caption = distributionName;
                    caption += "-";
                    caption += "TD";

                    string plotName("document-plot-");
                    plotName += removeSpaces(caption);
                    plotName += "-";
                    plotName += "runtime";

                    PgfplotPrinter tdPlot("./", plotName);
                    tdPlot.setCaption(caption);
                    tdPlot.plotAxis("runtime", tdResults, X_PLOT_SCALE, Y_PLOT_SCALE, false);
                    document.addToDocument(tdPlot);
                }

                // non td plots
                SpannerReducedLevelResultMap nontdResults;
                copy_if(distribution.begin(), distribution.end(),
                        inserter(nontdResults, nontdResults.end()),
                        [&nontdPlots](const auto &elem) {
                            return contains(nontdPlots, elem.first);
                        });

                if( !nontdResults.empty()) {
                    string caption = distributionName;
                    caption += "-";
                    caption += "nonTD";

                    string plotName = "document-plot-";
                    plotName += removeSpaces(caption);
                    plotName += "-";
                    plotName += "runtime";

                    PgfplotPrinter nontdPlot("./", plotName);
                    nontdPlot.setCaption(caption);
                    nontdPlot.plotAxis("runtime", nontdResults, X_PLOT_SCALE, Y_PLOT_SCALE, false);
                    document.addToDocument(nontdPlot);
                }

                // Linf plot
                SpannerReducedLevelResultMap linfResults;
                copy_if(distribution.begin(), distribution.end(),
                        inserter(linfResults, linfResults.end()),
                        [&linfPlots](const auto &elem) {
                            return contains(linfPlots, elem.first);
                        });

                if( !linfResults.empty()) {
                    string caption = distributionName;
                    caption += "-";

                    string plotName("document-plot-");
                    plotName += removeSpaces(caption);
                    plotName += "Linf";
                    plotName += "-";
                    plotName += "runtime";

                    caption += "$L_\\infty$";

                    PgfplotPrinter linfPlot("./", plotName);
                    linfPlot.setCaption(caption);
                    linfPlot.plotAxis("runtime", linfResults, X_PLOT_SCALE, Y_PLOT_SCALE, false);
                    document.addToDocument(linfPlot);
                }

                // L2 plot
                SpannerReducedLevelResultMap l2Results;
                copy_if(distribution.begin(), distribution.end(),
                        inserter(l2Results, l2Results.end()),
                        [&l2Plots](const auto &elem) {
                            return contains(l2Plots, elem.first);
                        });

                if( !l2Results.empty()) {
                    string caption = distributionName;
                    caption += "-";

                    string plotName("document-plot-");
                    plotName += removeSpaces(caption);
                    plotName += "L2";
                    plotName += "-";
                    plotName += "runtime";

                    caption += "$L_2$";

                    PgfplotPrinter l2Plot("./", plotName);
                    l2Plot.setCaption(caption);
                    l2Plot.plotAxis("runtime", l2Results, X_PLOT_SCALE, Y_PLOT_SCALE, false);
                    document.addToDocument(l2Plot);
                }

                for (const auto &iv: ANALYSIS_IVs) {
                    string caption = distributionName;

                    string plotName("document-plot-");
                    plotName += removeSpaces(caption);
                    plotName += "-";
                    plotName += iv;

                    PgfplotPrinter plot("./", plotName);
                    plot.setCaption(caption);
                    plot.plotAxis(iv, distribution, X_PLOT_SCALE, Y_PLOT_SCALE, isFirst);
                    document.addToDocument(plot);

                    if (isFirst) {
                        isFirst = false;
                        string legendRefText = plot.getLegend();
                        document.addRawText(legendRefText);
                        document.addRawText("\n\n");
                    }
                }
                document.clearpage();
        }
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
                                   return std::to_string(elem.first / X_PLOT_SCALE) + mathrm("M");
                               });

                vector<vector<string>> ivValues;
                for(auto spanner : distribution.second) {
                    headers.push_back(texttt(spanner.first));
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

                document.addRawText(subsection(caption) + "\n\n");
                document.addToDocument(table);
            }
        }
    }

    SpannerLevelResultMap
    getSpannerLevelResults(const vector<vector<string>>& rawResults) {
        SpannerLevelResultMap results;
        auto headerRow = rawResults.front();

        bool isFirst = true;

        for(auto row : rawResults) {
            if(isFirst) {
                isFirst = false;
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

            const auto& frontRow = summary.begin()->second;
            std::transform(frontRow.begin(), frontRow.end(),back_inserter(ivValues.front()),
                           [](const auto& elem){
                               return std::to_string(elem.first / X_PLOT_SCALE) + mathrm("M");
                           });
            for(auto spanner : summary) {
                headers.push_back(texttt(spanner.first));
                ivValues.push_back({});
                for(auto level : spanner.second) {
                    ivValues.back().push_back(std::to_string(level.second.getIV<double>(iv)));
                }
            }
            for(size_t i=0; i<headers.size();++i){
                table.addColumn( headers.at(i), ivValues.at(i),i);
            }
            table.tabulate();

            document.addRawText(subsection(caption) + "\n\n");
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
        copy_if(ALGORITHM_NAMES.begin(),ALGORITHM_NAMES.end(), back_inserter(body.back()),
            [&](const string& name){
                return contains(summary,name);
            });
        transform(body.back().begin(),body.back().end(),body.back().begin(),
            [](const string& name) {
                return "\\texttt{" + name + "}";
            });

        for(auto iv : ANALYSIS_IVs) {
            if(iv == "runtime") continue;

            body.push_back({});
            for(auto spanner : ALGORITHM_NAMES) {
                if(contains(summary,spanner))
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

        document.addRawText(subsection(caption) +  + "\n\n");
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
