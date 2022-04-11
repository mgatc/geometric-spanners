#ifndef GEOMETRIC_SPANNERS_ANALYSIS_H
#define GEOMETRIC_SPANNERS_ANALYSIS_H

#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "cpptex/cpptex.h"
#include "libspanner/BoundedDegreePlaneSpanners.h"
#include "libspanner/utilities.h"

#include "tools/Results.h"

namespace bdps_experiment {

namespace analysis {

    enum {
        SyntheticExperiment,
        RealExperiment,
        NoExperimentType,
    } EXPERIMENT_TYPE = NoExperimentType;

    std::string OUTPUT_DIRECTORY = "/tmp/";

    std::string X_PLOT_SCALE_UNIT = "";
    std::string X_PLOT_SCALE_SHORT_UNIT = "";
    double X_PLOT_SCALE = 1.0;
    double Y_PLOT_SCALE = 1.0;

    using namespace std;

    typedef size_t level_t;
    typedef string spanner_t;
    typedef string distribution_t;
    typedef GreedySpannerAnalysisResult result_t;

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

    std::vector<std::string> ANALYSIS_IVs = {
        "runtime","degree","degreeAvg","avgDegreePerPoint",
        "maxStretchFactor","avgStretchFactor","lightness"
    };



    std::vector<std::vector<std::string>>
    readFileIntoVector(const std::string &filename) {
        ifstream expIn;
        expIn.open(filename, ios_base::in);
        if (!expIn.is_open()) assert(!"Error opening file");

        vector<vector<string>> results;
        // skips first line! usually headers

        auto row = spanner::getNextLineAndSplitIntoTokens(expIn);
        try {
            stoi(row[1]); // second column should have a numeric value
            results.push_back(row);
        } catch( invalid_argument& ia) {
            cout<< "Throwing out header row";
        }
        while( !(row = spanner::getNextLineAndSplitIntoTokens(expIn)).empty()) {
            results.push_back(row);
        }
        expIn.close();

        return results;
    }



    DistributionSpannerLevelResultMap
    getDistributionSpannerLevelResults(const vector<vector<string>> &rawResults) {
        DistributionSpannerLevelResultMap results;
        for (auto row : rawResults) {
            result_t result;
            result.loadRow(row);

            results[result.distribution]
                   [result.algorithm]
                   [result.n].push_back(result);
        }
        return results;
    }

    SpannerLevelResultMap
    getSpannerLevelResults(const vector<vector<string>>& rawResults) {
        SpannerLevelResultMap results;

        for(auto row : rawResults) {
            result_t result;
            result.loadRow(row);

            results[result.algorithm]
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

    SpannerReducedLevelResultMap
    calculateSpannerLevelSummary(const SpannerLevelResultMap& results) {
        SpannerReducedLevelResultMap summary;

        for(const auto& spanner : results) {
            for(const auto& level : spanner.second) {
                result_t sum = std::accumulate(level.second.begin(),level.second.end(),result_t());
                size_t divisor = EXPERIMENT_TYPE == SyntheticExperiment ? level.second.size() : 1;
                result_t avg = sum / divisor;

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
    convertToResultMatrix(cpptex::PgfplotPrinter::ResultMatrix& out,
                          const std::string& iv,
                          const SpannerReducedLevelResultMap& results,
                          const double xScale = 1.0) {

        for( const auto& name : spanner::bdps::ALGORITHM_NAMES ) {
            if( spanner::contains(results,name) ) {
                out.emplace_back();
                for (const auto &level: results.at(name)) {
                    out.back().emplace_back(static_cast<double>(level.first) / xScale,
                                            level.second.getIV<double>(iv));
                }
            }
        }
    }

    void plot(const DistributionSpannerReducedLevelResultMap &summary,
              cpptex::LatexPrinter &document, bool lite = false) {

        set<string> tdPlots = {"BGHP2010", "KPT2017"};
        set<string> linfPlots = {"BKPX2015", "Degree3"};
        set<string> l2Plots;
        copy_if(spanner::bdps::ALGORITHM_NAMES.begin(),spanner::bdps::ALGORITHM_NAMES.end(),
                inserter(l2Plots,l2Plots.end()),
                [&](const string& name){
                    return !(spanner::contains(tdPlots,name) || spanner::contains(linfPlots,name));
                });

        set<string> nontdPlots;
        copy_if(spanner::bdps::ALGORITHM_NAMES.begin(), spanner::bdps::ALGORITHM_NAMES.end(),
                inserter(nontdPlots, nontdPlots.end()),
                [&tdPlots](const string &name) { // return true if name is not in tdPlots
                    return !spanner::contains(tdPlots, name);
                });

        int numCols = 2 + int(!lite);


        for(const auto& distributionName : spanner::SYNTHETIC_DISTRIBUTION_NAMES ) {
            if(summary.find(distributionName)!=summary.end()){
                const auto &distribution = summary.at(distributionName);

                string figureName("document-plot-");
                string caption = distributionName;
                caption += " Distribution Results";
                figureName += spanner::removeSpaces(caption);
                figureName += "-subfigure";

                cpptex::LatexPrinter parentFigure(OUTPUT_DIRECTORY + figureName);
                parentFigure.setCaption(caption);
                // split distribution.second into containers for td and non-td, l2 and nonl2

                // all results
                SpannerReducedLevelResultMap allResults;
                copy(distribution.begin(), distribution.end(),
                     inserter(allResults, allResults.end()));

                if( !lite && !allResults.empty()) {
                    string caption = distributionName;
                    caption += "-";
                    caption += "all";

                    string plotName("document-plot-");
                    plotName += spanner::removeSpaces(caption);
                    plotName += "-";
                    plotName += "runtime";

                    cpptex::PgfplotPrinter tdPlot(OUTPUT_DIRECTORY + plotName);
                    tdPlot.setCaption(caption);
                    cpptex::PgfplotPrinter::ResultMatrix resMatrix;
                    convertToResultMatrix(resMatrix, "runtime", allResults, X_PLOT_SCALE);
                    tdPlot.plotAxis(resMatrix, spanner::bdps::ALGORITHM_NAMES, "$n$", "Average execution time", caption);
                    parentFigure.addToDocumentAsSubfigure(tdPlot, numCols);
                }

                // non td plots
                SpannerReducedLevelResultMap nontdResults;
                copy_if(distribution.begin(), distribution.end(),
                        inserter(nontdResults, nontdResults.end()),
                        [&nontdPlots](const auto &elem) {
                            return spanner::contains(nontdPlots, elem.first);
                        });

                if( !nontdResults.empty()) {
                    string caption = distributionName;
                    caption += "-";
                    caption += "nonTD";

                    string plotName = "document-plot-";
                    plotName += spanner::removeSpaces(caption);
                    plotName += "-";
                    plotName += "runtime";

                    cpptex::PgfplotPrinter nontdPlot(OUTPUT_DIRECTORY + plotName);
                    nontdPlot.setCaption(caption);
                    cpptex::PgfplotPrinter::ResultMatrix resMatrix;
                    convertToResultMatrix(resMatrix, "runtime", nontdResults, X_PLOT_SCALE);
                    nontdPlot.plotAxis(resMatrix, spanner::bdps::ALGORITHM_NAMES, "$n$", "Average execution time", caption);
                    parentFigure.addToDocumentAsSubfigure( nontdPlot, numCols);
                }

                // Linf plot
                SpannerReducedLevelResultMap linfResults;
                copy_if(distribution.begin(), distribution.end(),
                        inserter(linfResults, linfResults.end()),
                        [&linfPlots](const auto &elem) {
                            return spanner::contains(linfPlots, elem.first);
                        });

                if( !linfResults.empty()) {
                    string caption = distributionName;
                    caption += "-";

                    string plotName("document-plot-");
                    plotName += spanner::removeSpaces(caption);
                    plotName += "Linf";
                    plotName += "-";
                    plotName += "runtime";

                    caption += "$L_\\infty$";

                    cpptex::PgfplotPrinter linfPlot(OUTPUT_DIRECTORY + plotName);
                    linfPlot.setCaption(caption);
                    cpptex::PgfplotPrinter::ResultMatrix resMatrix;
                    convertToResultMatrix(resMatrix, "runtime", linfResults, X_PLOT_SCALE);
                    linfPlot.plotAxis(resMatrix, spanner::bdps::ALGORITHM_NAMES, "$n$", "Average execution time", caption);
                    //parentFigure.addToDocumentAsSubfigure(linfPlot);
                }

                // L2 plot
                SpannerReducedLevelResultMap l2Results;
                copy_if(distribution.begin(), distribution.end(),
                        inserter(l2Results, l2Results.end()),
                        [&l2Plots](const auto &elem) {
                            return spanner::contains(l2Plots, elem.first);
                        });

                if( !l2Results.empty()) {
                    string caption = distributionName;
                    caption += "-";

                    string plotName("document-plot-");
                    plotName += spanner::removeSpaces(caption);
                    plotName += "L2";
                    plotName += "-";
                    plotName += "runtime";

                    caption += "$L_2$";

                    cpptex::PgfplotPrinter l2Plot(OUTPUT_DIRECTORY + plotName);
                    l2Plot.setCaption(caption);
                    cpptex::PgfplotPrinter::ResultMatrix resMatrix;
                    convertToResultMatrix(resMatrix, "runtime", l2Results, X_PLOT_SCALE);
                    l2Plot.plotAxis(resMatrix, spanner::bdps::ALGORITHM_NAMES, "$n$", "Average execution time", caption);
                    parentFigure.addToDocumentAsSubfigure(l2Plot, numCols);
                }

                if(!lite) {
                    for (const auto &iv: ANALYSIS_IVs) {
                        string caption = distributionName;

                        if (iv == "runtime") continue;

                        string plotName("document-plot-");
                        plotName += spanner::removeSpaces(caption);
                        plotName += "-";
                        plotName += iv;

                        cpptex::PgfplotPrinter plot(OUTPUT_DIRECTORY + plotName);
                        plot.setCaption(caption);
                        cpptex::PgfplotPrinter::ResultMatrix resMatrix;
                        convertToResultMatrix(resMatrix, iv, distribution, X_PLOT_SCALE);
                        plot.plotAxis(resMatrix, spanner::bdps::ALGORITHM_NAMES, "$n$", iv, caption);
                        parentFigure.addToDocumentAsSubfigure(plot, numCols);
                    }
                }
                document.addToDocumentAsFigure(parentFigure);
            }
        }
    }

    void plot(const SpannerReducedLevelResultMap& summary,
              cpptex::LatexPrinter& document, bool lite = false ) {
        bool isFirst = true;

        int numCols = 2 + int(!lite);

        set<string> tdPlots = {"BGHP2010", "KPT2017"};
        set<string> linfPlots = {"BKPX2015", "Degree3"};
        set<string> l2Plots;
        copy_if(spanner::bdps::ALGORITHM_NAMES.begin(),spanner::bdps::ALGORITHM_NAMES.end(),
                inserter(l2Plots,l2Plots.end()),
                [&](const string& name){
                    return !(spanner::contains(tdPlots,name) || spanner::contains(linfPlots,name));
                });

        set<string> nontdPlots;
        copy_if(spanner::bdps::ALGORITHM_NAMES.begin(), spanner::bdps::ALGORITHM_NAMES.end(),
                inserter(nontdPlots, nontdPlots.end()),
                [&tdPlots](const string &name) { // return true if name is not in tdPlots
                    return !spanner::contains(tdPlots, name);
                });


        string caption("Summary Results");
        string filename("document-plot-");
        filename += spanner::removeSpaces(caption) + "-subfigure";

        cpptex::LatexPrinter parentFigure(OUTPUT_DIRECTORY + filename);
        parentFigure.setCaption(caption);

        // all plots
        SpannerReducedLevelResultMap allResults;
        copy(summary.begin(), summary.end(),
             inserter(allResults, allResults.end()));

        if( !lite && !allResults.empty()) {
            string caption = "summary";
            caption += "-";
            caption += "all";

            string plotName("document-plot-");
            plotName += spanner::removeSpaces(caption);
            plotName += "-";
            plotName += "runtime";

            cpptex::PgfplotPrinter tdPlot(OUTPUT_DIRECTORY + plotName);
            tdPlot.setCaption(caption);
            cpptex::PgfplotPrinter::ResultMatrix resMatrix;
            convertToResultMatrix(resMatrix, "runtime", allResults, X_PLOT_SCALE);
            tdPlot.plotAxis(resMatrix, spanner::bdps::ALGORITHM_NAMES, "$n$", "Average execution time", caption);
            parentFigure.addToDocumentAsSubfigure(tdPlot, numCols);
        }

        // non td plots
        SpannerReducedLevelResultMap nontdResults;
        copy_if(summary.begin(), summary.end(),
                inserter(nontdResults, nontdResults.end()),
                [&nontdPlots](const auto &elem) {
                    return spanner::contains(nontdPlots, elem.first);
                });

        if( !nontdResults.empty()) {
            string caption = "summary";
            caption += "-";
            caption += "nonTD";

            string plotName = "document-plot-";
            plotName += spanner::removeSpaces(caption);
            plotName += "-";
            plotName += "runtime";

            cpptex::PgfplotPrinter nontdPlot(OUTPUT_DIRECTORY + plotName);
            nontdPlot.setCaption(caption);
            cpptex::PgfplotPrinter::ResultMatrix resMatrix;
            convertToResultMatrix(resMatrix, "runtime", nontdResults, X_PLOT_SCALE);
            nontdPlot.plotAxis(resMatrix, spanner::bdps::ALGORITHM_NAMES, "$n$", "Average execution time", caption);
            parentFigure.addToDocumentAsSubfigure( nontdPlot, numCols);
        }

        // Linf plot
        SpannerReducedLevelResultMap linfResults;
        copy_if(summary.begin(), summary.end(),
                inserter(linfResults, linfResults.end()),
                [&linfPlots](const auto &elem) {
                    return spanner::contains(linfPlots, elem.first);
                });

        if( !linfResults.empty()) {
            string caption = "summary";
            caption += "-";

            string plotName("document-plot-");
            plotName += spanner::removeSpaces(caption);
            plotName += "Linf";
            plotName += "-";
            plotName += "runtime";

            caption += "$L_\\infty$";

            cpptex::PgfplotPrinter linfPlot(OUTPUT_DIRECTORY + plotName);
            linfPlot.setCaption(caption);
            cpptex::PgfplotPrinter::ResultMatrix resMatrix;
            convertToResultMatrix(resMatrix, "runtime", linfResults, X_PLOT_SCALE);
            linfPlot.plotAxis(resMatrix, spanner::bdps::ALGORITHM_NAMES, "$n$", "Average execution time", caption);
//            parentFigure.addToDocumentAsSubfigure(linfPlot);
        }

        // L2 plot
        SpannerReducedLevelResultMap l2Results;
        copy_if(summary.begin(), summary.end(),
                inserter(l2Results, l2Results.end()),
                [&l2Plots](const auto &elem) {
                    return spanner::contains(l2Plots, elem.first);
                });

        if( !l2Results.empty()) {
            string caption = "summary";
            caption += "-";

            string plotName("document-plot-");
            plotName += spanner::removeSpaces(caption);
            plotName += "L2";
            plotName += "-";
            plotName += "runtime";

            caption += "$L_2$";

            cpptex::PgfplotPrinter l2Plot(OUTPUT_DIRECTORY + plotName);
            l2Plot.setCaption(caption);
            cpptex::PgfplotPrinter::ResultMatrix resMatrix;
            convertToResultMatrix(resMatrix, "runtime", l2Results, X_PLOT_SCALE);
            l2Plot.plotAxis(resMatrix, spanner::bdps::ALGORITHM_NAMES, "$n$", "Average execution time", caption);
            parentFigure.addToDocumentAsSubfigure(l2Plot, numCols);
        }

        if(!lite) {
            string plotLegendText;
            for (const auto &iv: ANALYSIS_IVs) {

                if (iv == "runtime")
                    continue;

                string plotName("document-plot-summary-");
                plotName += iv;
                string caption = "Summary";

                cpptex::PgfplotPrinter plot(OUTPUT_DIRECTORY + plotName);
                plot.setCaption(caption);
                cpptex::PgfplotPrinter::ResultMatrix resMatrix;
                convertToResultMatrix(resMatrix, iv, summary, X_PLOT_SCALE);
                plot.plotAxis(resMatrix, spanner::bdps::ALGORITHM_NAMES, "$n$", iv, caption);

                parentFigure.addToDocumentAsSubfigure(plot, numCols);

                if (isFirst) {
                    isFirst = false;
                    plotLegendText = plot.getLegend();
                }
            }

            cpptex::LatexPrinter legend(OUTPUT_DIRECTORY +"document-plot-legend");
            legend.addRawText(plotLegendText);
            document.addToDocumentAsFigure(legend);
        }
        document.addToDocumentAsFigure(parentFigure);
    }



    void tabulate(DistributionSpannerReducedLevelResultMap& spannerSummary,
                  cpptex::LatexPrinter& document, bool lite = false) {


        for(auto distribution : spannerSummary) {

            string caption = "Spanner Summary";
            string tableName("document-table-summary-");
            string distNoSpace = spanner::removeSpaces(distribution.first);
            tableName += distNoSpace;


                vector<string> degreeHeader;
                degreeHeader.push_back(spanner::bdps::ALGORITHM_SYMBOL);
                auto firstDegreeIv = find_if(ANALYSIS_IVs.begin(), ANALYSIS_IVs.end(), [](const auto& iv) {
                    return iv.find(string("egree")) != std::string::npos;
                });
                auto lastDegreeIv = find_if(next(firstDegreeIv), ANALYSIS_IVs.end(), [](const auto& iv) {
                   return iv.find("egree") == std::string::npos;
                });
                transform(firstDegreeIv, lastDegreeIv, back_inserter(degreeHeader),
                          [](const auto &iv) {
                              return cpptex::TablePrinter::m_ivNiceNames.find(iv) ==
                                     cpptex::TablePrinter::m_ivNiceNames.end() ?
                                     iv : cpptex::TablePrinter::m_ivNiceNames.at(iv);
                          });

                vector<vector<string>> degreeBody;
                degreeBody.push_back({});
                copy_if(spanner::bdps::ALGORITHM_NAMES.begin(), spanner::bdps::ALGORITHM_NAMES.end(),
                        back_inserter(degreeBody.back()),
                        [&](const string &name) {
                            return spanner::contains(distribution.second, name);
                        });
                transform(degreeBody.back().begin(), degreeBody.back().end(), degreeBody.back().begin(),
                          [](const string &name) {
                              return "\\texttt{" + name + "}";
                          });

                for (auto ivit=firstDegreeIv; ivit!=lastDegreeIv; ++ivit) {
                    auto iv = *ivit;
                    if (iv == "runtime") continue;


                    degreeBody.push_back({});

                    for (auto name: spanner::bdps::ALGORITHM_NAMES) {
                        if (spanner::contains(distribution.second, name)) {
                            auto sum = accumulate(distribution.second.at(name).begin(),
                                                  distribution.second.at(name).end(),
                                                  result_t(),
                                                  [](const result_t &sum, const auto &addend) {
                                                      return sum + addend.second;
                                                  });
                            size_t total = distribution.second.at(name).size();
                            degreeBody.back().push_back(std::to_string((sum / total).template getIV<double>(iv)));
                        }
                    }


                }

                auto degreeTableName = tableName + "-degree";
                cpptex::TablePrinter degreeTable(OUTPUT_DIRECTORY + degreeTableName);

                for (size_t i = 0; i < degreeBody.size(); i++) {
                    int precision = degreeHeader[i] == cpptex::TablePrinter::m_ivNiceNames["degree"] ? 0 : 3;
                    degreeTable.addColumn(degreeHeader[i], degreeBody[i], precision, i);
                }
                degreeTable.tabulate();

                document.addToDocument(degreeTable);



            vector<string> sfHeader;
            sfHeader.push_back(spanner::bdps::ALGORITHM_SYMBOL);
            auto firstSfIv = lastDegreeIv;
            auto lastSfIv = find_if(next(firstSfIv), ANALYSIS_IVs.end(), [](const auto& iv) {
                return iv.find("lightness") != std::string::npos;
            });
            transform(firstSfIv, lastSfIv, back_inserter(sfHeader),
                      [](const auto &iv) {
                          return cpptex::TablePrinter::m_ivNiceNames.find(iv) ==
                                 cpptex::TablePrinter::m_ivNiceNames.end() ?
                                 iv : cpptex::TablePrinter::m_ivNiceNames.at(iv);
                      });

            vector<vector<string>> sfBody;
            sfBody.push_back({});
            copy_if(spanner::bdps::ALGORITHM_NAMES.begin(), spanner::bdps::ALGORITHM_NAMES.end(),
                    back_inserter(sfBody.back()),
                    [&](const string &name) {
                        return spanner::contains(distribution.second, name);
                    });
            transform(sfBody.back().begin(), sfBody.back().end(), sfBody.back().begin(),
                      [](const string &name) {
                          return "\\texttt{" + name + "}";
                      });

            for (auto ivit=firstSfIv; ivit!=lastSfIv; ++ivit) {
                auto iv = *ivit;
                if (iv == "runtime") continue;


                sfBody.push_back({});

                for (auto name: spanner::bdps::ALGORITHM_NAMES) {
                    if (spanner::contains(distribution.second, name)) {
                        auto sum = accumulate(distribution.second.at(name).begin(),
                                              distribution.second.at(name).end(),
                                              result_t(),
                                              [](const result_t &sum, const auto &addend) {
                                                  return sum + addend.second;
                                              });
                        size_t total = distribution.second.at(name).size();
                        sfBody.back().push_back(std::to_string((sum / total).template getIV<double>(iv)));
                    }
                }


            }

            auto sfTableName = tableName + "-sf";
            cpptex::TablePrinter sfTable(OUTPUT_DIRECTORY + sfTableName);

            for (size_t i = 0; i < sfBody.size(); i++) {
                int precision = sfHeader[i] == cpptex::TablePrinter::m_ivNiceNames["degree"] ? 0 : 3;
                sfTable.addColumn(sfHeader[i], sfBody[i], precision, i);
            }
            sfTable.tabulate();

            document.addToDocument(sfTable);





            vector<string> lightnessHeader;
            lightnessHeader.push_back(spanner::bdps::ALGORITHM_SYMBOL);
            auto firstLightnessIv = lastSfIv;
            auto lastLightnessIv = ANALYSIS_IVs.end();
            transform(firstLightnessIv, lastLightnessIv, back_inserter(lightnessHeader),
                      [](const auto &iv) {
                          return cpptex::TablePrinter::m_ivNiceNames.find(iv) ==
                                 cpptex::TablePrinter::m_ivNiceNames.end() ?
                                 iv : cpptex::TablePrinter::m_ivNiceNames.at(iv);
                      });

            vector<vector<string>> lightnessBody;
            lightnessBody.push_back({});
            copy_if(spanner::bdps::ALGORITHM_NAMES.begin(), spanner::bdps::ALGORITHM_NAMES.end(),
                    back_inserter(lightnessBody.back()),
                    [&](const string &name) {
                        return spanner::contains(distribution.second, name);
                    });
            transform(lightnessBody.back().begin(), lightnessBody.back().end(), lightnessBody.back().begin(),
                      [](const string &name) {
                          return "\\texttt{" + name + "}";
                      });

            for (auto ivit=firstLightnessIv; ivit!=lastLightnessIv; ++ivit) {
                auto iv = *ivit;
                if (iv == "runtime") continue;


                lightnessBody.push_back({});

                for (auto name: spanner::bdps::ALGORITHM_NAMES) {
                    if (spanner::contains(distribution.second, name)) {
                        auto sum = accumulate(distribution.second.at(name).begin(),
                                              distribution.second.at(name).end(),
                                              result_t(),
                                              [](const result_t &sum, const auto &addend) {
                                                  return sum + addend.second;
                                              });
                        size_t total = distribution.second.at(name).size();
                        lightnessBody.back().push_back(std::to_string((sum / total).template getIV<double>(iv)));
                    }
                }


            }

            auto lightnessTableName = tableName + "-lightness";
            cpptex::TablePrinter lightnessTable(OUTPUT_DIRECTORY + lightnessTableName);

            for (size_t i = 0; i < lightnessBody.size(); i++) {
                int precision = lightnessHeader[i] == cpptex::TablePrinter::m_ivNiceNames["degree"] ? 0 : 3;
                lightnessTable.addColumn(lightnessHeader[i], lightnessBody[i], precision, i);
            }
            lightnessTable.tabulate();

            document.addToDocument(lightnessTable);



            document.clearpage();
        }
    }

    void tabulate(const SpannerReducedLevelResultMap& summary,
                  cpptex::LatexPrinter& document, bool lite = false) {
        int tableCount = 0;
        for(auto iv : ANALYSIS_IVs) {

            string caption(iv);
            string tableName("document-table-summary");
            tableName += spanner::removeSpaces(caption);
            cpptex::TablePrinter table(OUTPUT_DIRECTORY + tableName);

            // add columns
            vector<string> headers;
            vector<vector<string>> ivValues;

            headers.push_back(N_SYMBOL);
            ivValues.push_back({});

            const auto& frontRow = summary.begin()->second;
            cpptex::SetPrecision precisionSetter{int(X_PLOT_SCALE==1000000)};
            transform(frontRow.begin(), frontRow.end(),back_inserter(ivValues.front()),
                [precisionSetter](const auto& elem){
                    return precisionSetter(std::to_string(elem.first / X_PLOT_SCALE)) + cpptex::mathrm(X_PLOT_SCALE_SHORT_UNIT);
                });

            for(auto name : spanner::bdps::ALGORITHM_NAMES) {
                if( spanner::contains(summary, name)) {
                    auto& spanner = summary.at(name);
                    headers.push_back(cpptex::texttt(name));
                    ivValues.push_back({});
                    for(auto level : spanner) {
                        ivValues.back().push_back(std::to_string(level.second.getIV<double>(iv)));
                    }
                }
            }


            for(size_t i=0; i<headers.size(); ++i){
                int precision = i==0 ? -1 :
                        iv == "degree" ? 0 : 3;
                table.addColumn( headers.at(i), ivValues.at(i),precision,i);
            }
            table.tabulate(true);
            table.setCaption(caption);

            document.addToDocument(table);
        }
    }

    void tabulate(const SpannerResultMap& summary,
                  cpptex::LatexPrinter& document, bool lite = false) {
        if(lite)
            return;

        vector<string> header;
        header.push_back(spanner::bdps::ALGORITHM_SYMBOL);
        transform(next(ANALYSIS_IVs.begin()),ANALYSIS_IVs.end(),back_inserter(header),
            [](const auto& iv){
                return cpptex::TablePrinter::m_ivNiceNames.find(iv) == cpptex::TablePrinter::m_ivNiceNames.end() ?
                        iv : cpptex::TablePrinter::m_ivNiceNames.at(iv);
            });

        vector<vector<string>> body;
        body.push_back({});
        copy_if(spanner::bdps::ALGORITHM_NAMES.begin(),spanner::bdps::ALGORITHM_NAMES.end(), back_inserter(body.back()),
            [&](const string& name){
                return spanner::contains(summary,name);
            });
        transform(body.back().begin(),body.back().end(),body.back().begin(),
            [](const string& name) {
                return "\\texttt{" + name + "}";
            });

        for(auto iv : ANALYSIS_IVs) {
            if(iv == "runtime") continue;

            body.push_back({});
            for(auto spanner : spanner::bdps::ALGORITHM_NAMES) {
                if(spanner::contains(summary,spanner))
                    body.back().push_back(std::to_string(summary.at(spanner).template getIV<double>(iv)));
            }
        }
        string caption = "Spanner Summary";
        string tableName("document-table-summary");
        cpptex::TablePrinter table(OUTPUT_DIRECTORY + tableName);

        for(size_t i=0; i<body.size(); i++) {
            int precision = header[i] == cpptex::TablePrinter::m_ivNiceNames["degree"] ? 0 : 3;
            table.addColumn(header[i], body[i],precision,i);
        }
        table.tabulate();

        document.addRawText(cpptex::subsection(caption) +  + "\n\n");
        document.addToDocument(table);
    }





void BoundedDegreePlaneSpannerAnalysisSynthetic(const string& filename) {
    using namespace bdps_experiment::analysis;

    cout << "Starting synthetic analysis..."<<endl;


    cpptex::LatexPrinter document("/tmp/bdps_analysis");
    vector<vector<string>> raw = readFileIntoVector(filename);

    int firstLevelInResults = stod(*next(raw.front().begin()));

    bool isLargeExperiment = firstLevelInResults >= 1e5;

    if(isLargeExperiment) {
        X_PLOT_SCALE = 1e6;
        X_PLOT_SCALE_UNIT = "millions";
        X_PLOT_SCALE_SHORT_UNIT = "M";
        ANALYSIS_IVs.resize(1);
    } else {
        X_PLOT_SCALE = 1e3;
        X_PLOT_SCALE_UNIT = "thousands";
        X_PLOT_SCALE_SHORT_UNIT = "K";
    }

    DistributionSpannerLevelResultMap distributionSpannerLevelResults = getDistributionSpannerLevelResults(raw);
    DistributionSpannerReducedLevelResultMap distributionSpannerSummary = calculateDistributionSpannerSummary(distributionSpannerLevelResults);

    SpannerLevelResultMap spannerLevelResults = getSpannerLevelResults(raw);

//    auto numLevels = spannerLevelResults.size();
//    auto samplesPerLevel = spannerLevelResults.begin()->second.size();
//    auto numSamples = numLevels * samplesPerLevel;


    cout<<"Finding average per level for all distributions..."<<endl;
    document.addRawText("\\section{Overall Summary}\n\n");

    SpannerReducedLevelResultMap spannerLevelSummary = calculateSpannerLevelSummary(spannerLevelResults);
    SpannerResultMap spannerSummary = calculateSpannerSummary(spannerLevelSummary);

    tabulate(spannerSummary, document, isLargeExperiment);
//    plot(spannerLevelSummary, document, isLargeExperiment);
//    tabulate(spannerLevelSummary, document, isLargeExperiment);

    if(!isLargeExperiment) document.clearpage();

    cout<< "Finding summaries for each spanner algorithm within each distribution..."<<endl;
    document.addRawText("\\section{Distribution Summaries}\n\n");

//    plot(distributionSpannerSummary, document, isLargeExperiment);
    tabulate(distributionSpannerSummary, document, isLargeExperiment);




//    document.clearpage();

    //tabulateOverallSummary(spannerLevelResults, document);
    document.save();
    //document.display();
}

void BoundedDegreePlaneSpannerAnalysisReal(const string& filename) {
    using namespace bdps_experiment::analysis;

    cout << "Starting real analysis..."<<endl;

    cpptex::LatexPrinter document(OUTPUT_DIRECTORY + "document");
    vector<vector<string>> raw = readFileIntoVector(filename);

//    cout<< "Finding summaries for each spanner algorithm within each distribution..."<<endl;
//    DistributionSpannerLevelResultMap distributionSpannerLevelResults = getDistributionSpannerLevelResults(raw);
//    DistributionSpannerReducedLevelResultMap distributionSpannerSummary = calculateDistributionSpannerSummary(distributionSpannerLevelResults);
//    plot(distributionSpannerSummary, document);
//    tabulate(distributionSpannerSummary, document);

    cout<<"Finding average per level for all distributions..."<<endl;
    document.addRawText("\\section{Overall Summary}\n\n");
    SpannerLevelResultMap spannerLevelResults = getSpannerLevelResults(raw);

//    auto numLevels = spannerLevelResults.size();
//    auto samplesPerLevel = spannerLevelResults.begin()->second.size();
//    auto numSamples = numLevels * samplesPerLevel;


    SpannerReducedLevelResultMap spannerLevelSummary = calculateSpannerLevelSummary(spannerLevelResults);
    //plot(spannerLevelSummary, document);
    tabulate(spannerLevelSummary, document);

    SpannerResultMap spannerSummary = calculateSpannerSummary(spannerLevelSummary);
    tabulate(spannerSummary,document);


//    tabulateDistributionSummaries(distributionSpannerSummary, document);
    //plotDistributionSummaries();

//    document.clearpage();

    //tabulateOverallSummary(spannerLevelResults, document);

    document.save();
}

} // analysis



void BoundedDegreePlaneSpannerAnalysis(const std::string& filename) {
    using namespace analysis;
    cout << "Data analysis started... type: ";

    // get first token of filename, separated by '-'
    const char slash = '/';
    const auto firstAlphanumeric = filename.find_last_of(slash)+1;
    const char dash = '-';
    const auto dashPosition = filename.find_first_of(dash);
    const size_t expTypeStringLength = dashPosition - firstAlphanumeric;
    string experimentType = filename.substr(firstAlphanumeric, expTypeStringLength);
    cout << experimentType << endl;

    if(experimentType == "synthetic") {
        EXPERIMENT_TYPE = SyntheticExperiment;
        BoundedDegreePlaneSpannerAnalysisSynthetic(filename);
    } else if(experimentType == "real") {
        EXPERIMENT_TYPE = RealExperiment;
        BoundedDegreePlaneSpannerAnalysisReal(filename);
    } else {
        cout<< "Invalid experiment type. Exiting...";
    }
}


} // bdps_experiment



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
//        TablePrinter table(OUTPUT_DIRECTORY + tableName);
//
//        // add columns
//
//        for(size_t i=0; i<headers.size();++i){
//            table.addColumn( headers.at(i), ivValues.at(i),6,i);
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
//            // add max degreeMax, avg degreeMax, avg avgDegreePerPoint
//        }
//
//        TablePrinter table(OUTPUT_DIRECTORY +"document-table-overall");
//
//        for(auto header : headers) {
//            if(contains(body,header))
//                table.addColumn(header, body[header],6);
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
