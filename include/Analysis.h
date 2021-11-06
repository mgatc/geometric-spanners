#ifndef GEOMETRIC_SPANNERS_ANALYSIS_H
#define GEOMETRIC_SPANNERS_ANALYSIS_H

#include <fstream>
#include <iterator>
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

    enum {
        SyntheticExperiment,
        RealExperiment,
        NoExperimentType,
    } EXPERIMENT_TYPE = NoExperimentType;

    string OUTPUT_DIRECTORY = "/tmp/spanners/";

    string X_PLOT_SCALE_UNIT = "";
    string X_PLOT_SCALE_SHORT_UNIT = "";
    double X_PLOT_SCALE = 1.0;
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

    std::vector<string> ANALYSIS_IVs = {
        "runtime","degree","degreeAvg","avgDegreePerPoint",
        "maxStretchFactor","avgStretchFactor","lightness"
    };



    vector<vector<string>>
    readFileIntoVector(const string &filename) {
        ifstream expIn;
        expIn.open(filename, ios_base::in);
        if (!expIn.is_open()) assert(!"Error opening file");

        vector<vector<string>> results;
        // skips first line! usually headers
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



    void plot(const DistributionSpannerReducedLevelResultMap &summary,
              LatexPrinter &document, bool lite = false) {

        set<string> tdPlots = {"BGHP2010", "KPT2017"};
        set<string> linfPlots = {"BKPX2015", "Degree3"};
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


        for(const auto& distributionName : SYNTHETIC_DISTRIBUTION_NAMES ) {
            if(summary.find(distributionName)!=summary.end()){
                const auto &distribution = summary.at(distributionName);

                string figureName("document-plot-");
                string caption = distributionName + " Distribution Results";
                figureName += removeSpaces(caption);
                figureName += "-subfigure";

                LatexPrinter parentFigure(OUTPUT_DIRECTORY, figureName);
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
                    plotName += removeSpaces(caption);
                    plotName += "-";
                    plotName += "runtime";

                    PgfplotPrinter tdPlot(OUTPUT_DIRECTORY, plotName);
                    tdPlot.setCaption(caption);
                    tdPlot.plotAxis("runtime", allResults, X_PLOT_SCALE, X_PLOT_SCALE_UNIT, false);
                    parentFigure.addToDocumentAsSubfigure(tdPlot);
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

                    PgfplotPrinter nontdPlot(OUTPUT_DIRECTORY, plotName);
                    nontdPlot.setCaption(caption);
                    nontdPlot.plotAxis("runtime", nontdResults, X_PLOT_SCALE, X_PLOT_SCALE_UNIT, false);
                    parentFigure.addToDocumentAsSubfigure(nontdPlot);
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

                    PgfplotPrinter linfPlot(OUTPUT_DIRECTORY, plotName);
                    linfPlot.setCaption(caption);
                    linfPlot.plotAxis("runtime", linfResults, X_PLOT_SCALE, X_PLOT_SCALE_UNIT, false);
                    //parentFigure.addToDocumentAsSubfigure(linfPlot);
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

                    PgfplotPrinter l2Plot(OUTPUT_DIRECTORY, plotName);
                    l2Plot.setCaption(caption);
                    l2Plot.plotAxis("runtime", l2Results, X_PLOT_SCALE, X_PLOT_SCALE_UNIT, false);
                    parentFigure.addToDocumentAsSubfigure(l2Plot);
                }

                if(!lite) {
                    for (const auto &iv: ANALYSIS_IVs) {
                        string caption = distributionName;

                        if (iv == "runtime") continue;

                        string plotName("document-plot-");
                        plotName += removeSpaces(caption);
                        plotName += "-";
                        plotName += iv;

                        PgfplotPrinter plot(OUTPUT_DIRECTORY, plotName);
                        plot.setCaption(caption);
                        plot.plotAxis(iv, distribution, X_PLOT_SCALE, X_PLOT_SCALE_UNIT);
                        parentFigure.addToDocumentAsSubfigure(plot);
                    }
                }
                document.addToDocumentAsFigure(parentFigure);
            }
        }
    }

    void plot(const SpannerReducedLevelResultMap& summary,
         LatexPrinter& document, bool lite = false ) {
        bool isFirst = true;

        set<string> tdPlots = {"BGHP2010", "KPT2017"};
        set<string> linfPlots = {"BKPX2015", "Degree3"};
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


        string caption("Summary Results");
        string filename("document-plot-");
        filename += removeSpaces(caption) + "-subfigure";

        LatexPrinter parentFigure(OUTPUT_DIRECTORY, filename);
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
            plotName += removeSpaces(caption);
            plotName += "-";
            plotName += "runtime";

            PgfplotPrinter tdPlot(OUTPUT_DIRECTORY, plotName);
            tdPlot.setCaption(caption);
            tdPlot.plotAxis("runtime", allResults, X_PLOT_SCALE, X_PLOT_SCALE_UNIT, false);
            parentFigure.addToDocumentAsSubfigure(tdPlot);
        }

        // non td plots
        SpannerReducedLevelResultMap nontdResults;
        copy_if(summary.begin(), summary.end(),
                inserter(nontdResults, nontdResults.end()),
                [&nontdPlots](const auto &elem) {
                    return contains(nontdPlots, elem.first);
                });

        if( !nontdResults.empty()) {
            string caption = "summary";
            caption += "-";
            caption += "nonTD";

            string plotName = "document-plot-";
            plotName += removeSpaces(caption);
            plotName += "-";
            plotName += "runtime";

            PgfplotPrinter nontdPlot(OUTPUT_DIRECTORY, plotName);
            nontdPlot.setCaption(caption);
            nontdPlot.plotAxis("runtime", nontdResults, X_PLOT_SCALE, X_PLOT_SCALE_UNIT, false);
            parentFigure.addToDocumentAsSubfigure(nontdPlot);
        }

        // Linf plot
        SpannerReducedLevelResultMap linfResults;
        copy_if(summary.begin(), summary.end(),
                inserter(linfResults, linfResults.end()),
                [&linfPlots](const auto &elem) {
                    return contains(linfPlots, elem.first);
                });

        if( !linfResults.empty()) {
            string caption = "summary";
            caption += "-";

            string plotName("document-plot-");
            plotName += removeSpaces(caption);
            plotName += "Linf";
            plotName += "-";
            plotName += "runtime";

            caption += "$L_\\infty$";

            PgfplotPrinter linfPlot(OUTPUT_DIRECTORY, plotName);
            linfPlot.setCaption(caption);
            linfPlot.plotAxis("runtime", linfResults, X_PLOT_SCALE, X_PLOT_SCALE_UNIT, false);
            //parentFigure.addToDocumentAsSubfigure(linfPlot);
        }

        // L2 plot
        SpannerReducedLevelResultMap l2Results;
        copy_if(summary.begin(), summary.end(),
                inserter(l2Results, l2Results.end()),
                [&l2Plots](const auto &elem) {
                    return contains(l2Plots, elem.first);
                });

        if( !l2Results.empty()) {
            string caption = "summary";
            caption += "-";

            string plotName("document-plot-");
            plotName += removeSpaces(caption);
            plotName += "L2";
            plotName += "-";
            plotName += "runtime";

            caption += "$L_2$";

            PgfplotPrinter l2Plot(OUTPUT_DIRECTORY, plotName);
            l2Plot.setCaption(caption);
            l2Plot.plotAxis("runtime", l2Results, X_PLOT_SCALE, X_PLOT_SCALE_UNIT, false);
            parentFigure.addToDocumentAsSubfigure(l2Plot);
        }

        if(!lite) {
            string plotLegendText;
            for (const auto &iv: ANALYSIS_IVs) {

                if (iv == "runtime")
                    continue;

                string plotName("document-plot-summary-");
                plotName += iv;
                string caption = "Summary";

                PgfplotPrinter plot(OUTPUT_DIRECTORY, plotName);
                plot.setCaption(caption);
                plot.plotAxis(iv, summary, X_PLOT_SCALE, X_PLOT_SCALE_UNIT, isFirst);

                parentFigure.addToDocumentAsSubfigure(plot);

                if (isFirst) {
                    isFirst = false;
                    plotLegendText = plot.getLegend();
                }
            }
            document.addToDocumentAsFigure(parentFigure);

            LatexPrinter legend(OUTPUT_DIRECTORY,"document-plot-legend");
            legend.addRawText(plotLegendText);
            document.addToDocumentAsFigure(legend);
        }
    }



    void tabulate(DistributionSpannerReducedLevelResultMap& spannerSummary,
             LatexPrinter& document, bool lite = false) {

        for(auto distribution : spannerSummary) {
            for(auto iv : ANALYSIS_IVs) {
                vector<string> levels;
                const auto& frontRow = distribution.second.begin()->second;
                SetPrecision precisionSetter{int(X_PLOT_SCALE==1000000)};
                transform(frontRow.begin(), frontRow.end(), back_inserter(levels),
                          [precisionSetter](const auto& elem){
                              return precisionSetter(std::to_string(elem.first / X_PLOT_SCALE)) + mathrm(X_PLOT_SCALE_SHORT_UNIT);
                          });

                vector<string> headers;
                vector<vector<string>> ivValues;
                for(auto name : ALGORITHM_NAMES) {
                    if( contains(distribution.second, name)) {
                        auto& spanner = distribution.second[name];
                        headers.push_back(texttt(name));
                        ivValues.push_back({});
                        transform(spanner.begin(), spanner.end(), back_inserter(ivValues.back()),
                                  [iv](const auto& elem){
                                      return std::to_string(elem.second.template getIV<double>(iv));
                                  });
                    }
                }

                string caption(distribution.first + " x " + iv);
                string tableName("document-table-");
                tableName += removeSpaces(caption);
                TablePrinter table(OUTPUT_DIRECTORY, tableName);

                // add columns
                table.addColumn(N_SYMBOL,levels,-1,0);

                for(size_t i=0; i<headers.size();++i){
                    int precision = iv == "degree" ? 0 : 3;
                    table.addColumn( headers.at(i), ivValues.at(i),precision,i+1);
                }

                table.tabulate();

                document.addRawText(subsection(caption) + "\n\n");
                document.addToDocument(table);
            }
        }
    }

    void tabulate(const SpannerReducedLevelResultMap& summary,
             LatexPrinter& document, bool lite = false) {
        for(auto iv : ANALYSIS_IVs) {

            string caption(iv);
            string tableName("document-table-summary");
            tableName += removeSpaces(caption);
            TablePrinter table(OUTPUT_DIRECTORY, tableName);

            // add columns
            vector<string> headers;
            vector<vector<string>> ivValues;

            headers.push_back(N_SYMBOL);
            ivValues.push_back({});

            const auto& frontRow = summary.begin()->second;
            SetPrecision precisionSetter{int(X_PLOT_SCALE==1000000)};
            transform(frontRow.begin(), frontRow.end(),back_inserter(ivValues.front()),
                [precisionSetter](const auto& elem){
                    return precisionSetter(std::to_string(elem.first / X_PLOT_SCALE)) + mathrm(X_PLOT_SCALE_SHORT_UNIT);
                });

            for(auto name : ALGORITHM_NAMES) {
                if( contains(summary, name)) {
                    auto& spanner = summary.at(name);
                    headers.push_back(texttt(name));
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
            table.tabulate();

            document.addRawText(subsection(caption) + "\n\n");
            document.addToDocument(table);
        }
    }

    void tabulate(const SpannerResultMap& summary,
             LatexPrinter& document, bool lite = false) {
        if(lite)
            return;

        vector<string> header;
        header.push_back(ALGORITHM_SYMBOL);
        transform(next(ANALYSIS_IVs.begin()),ANALYSIS_IVs.end(),back_inserter(header),
            [](const auto& iv){
                return TablePrinter::m_ivNiceNames.find(iv) == TablePrinter::m_ivNiceNames.end() ?
                        iv : TablePrinter::m_ivNiceNames.at(iv);
            });

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
        TablePrinter table(OUTPUT_DIRECTORY, tableName);

        for(size_t i=0; i<body.size(); i++) {
            int precision = header[i] == TablePrinter::m_ivNiceNames["degree"] ? 0 : 3;
            table.addColumn(header[i], body[i],precision,i);
        }
        table.tabulate();

        document.addRawText(subsection(caption) +  + "\n\n");
        document.addToDocument(table);
    }





void BoundedDegreePlaneSpannerAnalysisSynthetic(const string& filename) {
    using namespace spanners::analysis;

    cout << "Starting synthetic analysis..."<<endl;


    LatexPrinter document(OUTPUT_DIRECTORY, "document");
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
    plot(spannerLevelSummary, document, isLargeExperiment);
    tabulate(spannerLevelSummary, document, isLargeExperiment);

    if(!isLargeExperiment) document.clearpage();

    cout<< "Finding summaries for each spanner algorithm within each distribution..."<<endl;
    document.addRawText("\\section{Distribution Summaries}\n\n");

    plot(distributionSpannerSummary, document, isLargeExperiment);
    tabulate(distributionSpannerSummary, document, isLargeExperiment);



//    tabulateDistributionSummaries(distributionSpannerSummary, document);
    //plotDistributionSummaries();

//    document.clearpage();

    //tabulateOverallSummary(spannerLevelResults, document);
    document.save();
    //document.display();
}

void BoundedDegreePlaneSpannerAnalysisReal(const string& filename) {
    using namespace spanners::analysis;

    cout << "Starting real analysis..."<<endl;

    LatexPrinter document(OUTPUT_DIRECTORY, "document");
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



void BoundedDegreePlaneSpannerAnalysis(const string& filename) {
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
//        TablePrinter table(OUTPUT_DIRECTORY, tableName);
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
//            // add max degreeMax, avg degreeMax, avg degreeAvg
//        }
//
//        TablePrinter table(OUTPUT_DIRECTORY,"document-table-overall");
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
