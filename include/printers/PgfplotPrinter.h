#ifndef SPANNERS_PGFPLOTSPRINTER_H
#define SPANNERS_PGFPLOTSPRINTER_H

#include <map>
#include <iostream>
#include <optional>
#include <vector>

#include "printers/GraphPrinter.h"
#include "Names.h"
#include "tools/Results.h"
#include "tools/Utilities.h"

namespace spanners {

    using namespace std;

    const size_t SCALING_FACTOR=1000;

//    void plotResults(const BoundedDegreeSpannerResultSet &results, LatexPrinter* addToPrinter);

    class PgfplotPrinter : public TikzPrinter {
    public:
        // Palette generated using https://coolors.co/
        vector<string> Colors = {
                "f94144","f3722c","f8961e",
                "f9844a","f9c74f","90be6d",
                "43aa8b","4d908e","577590",
                "277da1","360568"
        };

        std::vector<std::string> Marks = {
                //"otimes*",
                "*", "triangle*", "square*",  "pentagon*", "diamond*"
        };

        PgfplotPrinter(string directory, string filename, string documentType = "standalone")
                : TikzPrinter(directory,filename,documentType){
            m_body = Body{getTikzHeader(),
                          "",
                          getTikzFooter() };
            for(auto color : Colors){
                defineColor(color);
            }
        }





//        void plotAxis(unsigned iv, const BoundedDegreeSpannerResultSet &results, bool legend) {
//
//            // build axis header
//            string allPlotsOfAxis = getAxisHeader(iv, results, legend);
//
//            for( const auto& alg : results.m_reducedSamples ) {
//                // build the plot
//                string plot = getPlotHeader(alg.first);
//
//                for( const auto& level : alg.second ) {
//                    // add the level, measurement pair
//                    string entry = "("
//                                   + std::to_string(static_cast<double>(level.first) / SCALING_FACTOR)
//                                   + ","
//                                   + spanners::to_string(level.second.IV.at(IV_NAMES[iv]))
//                                   + ")\n";
//                    plot += entry;
//                }
//                plot += getPlotFooter(iv);
//
////                if(includeLegend)
////                    plot += addLegendEntry(alg.first);
//
//                allPlotsOfAxis += plot;
//            }
//            allPlotsOfAxis += getAxisFooter();
//
////            if(legend)
////                allPlotsOfAxis += getLegend();
//
//            m_body.content = allPlotsOfAxis;
//        }
//        void plotResults(const BoundedDegreeSpannerResultSet &results) {
//            spanners::plotResults(results, this);
//        }
//        string addLegendEntry(Algorithm alg) {
//            string legendEntry = "\\addlegendentry{\\textsc{\\tiny{"
//                                 + ALGORITHM_NAMES[alg]
//                                 + "}}}";
//            return legendEntry;
//        }
        string getLegendText( const vector<string>& names) {
            //TikzPrinter tikz("templegend");
            string legend = string("")
                            +"$\\matrix[\n"
                            +"matrix of nodes,\n"
                            +"draw,\n"
                            +"anchor=north,\n"
                            +"column 1/.style={nodes={anchor=center}},\n"
                            +"column 2/.style={nodes={anchor=west}}\n"
                            +"] at([yshift=-2.5em]current axis.south)\n"
                            +"{\n\n";
            for(auto name : names) {
                legend += "\\ref{plots:"
                          + name
                          + "} & $\\textbf{"
                          + name
                          + "}\\\\\n";
            }
            //    +"\\ref{plots:gamma}&$\\textbf{Gamma}~\\alpha=1.7435~ \\beta=21.263$\\\\"
            legend += "}$";

            return legend;
        }
        string getLegend() {
            return "\\ref{planespannerlegend}\n\n";
        }
        string getPlotHeader(Algorithm alg) {
            string plotHeader = "\n\n\\addplot+";
            plotHeader += string("[solid")
                    + ",color="
                    + Colors[alg%Colors.size()];
            plotHeader += string(",")
                          + "mark="
                          + Marks[alg % Marks.size()];
            plotHeader += "]";
            plotHeader += " coordinates {\n";
            return plotHeader;
        }
        string getPlotFooter(unsigned iv) {
            string plotFooter("}node [pos=1.15, above left] {};\n\n");
//            plotFooter += "\\label{plots:"
//                        + ALGORITHM_NAMES[iv]
//                        + "}";
            return plotFooter;
        }
//        string getAxisHeader(unsigned iv,
//                             const BoundedDegreeSpannerResultSet &results,
//                             bool legend ) {
//
//            string axisHeader = string("")
//                                + "\\begin{axis}[";
//            if(!m_caption.empty()) axisHeader += "title={" + m_caption + "},";
//            axisHeader += string("scaled ticks=false,grid=major,xlabel={$n$ (in thousands)}, ")
//                        + "xtick=";//{";
//            string xTicks = "data";
//            //assert(!results.m_reducedSamples.empty());
////            for( auto level : results.m_reducedSamples.begin()->second ){
////                xTicks += std::to_string(static_cast<double>(level.first) / SCALING_FACTOR);
////                xTicks += ",";
////            }
////            xTicks = xTicks.substr( 0, xTicks.size()-1 );
//
//            axisHeader += xTicks
//                          //+ "\n}"
//                          + ", ylabel near ticks,ylabel={"
//                          + IV_NICE_NAMES.at(iv);
//            if(!IV_UNITS.at(iv).empty())
//                axisHeader += " (in "
//                          + IV_UNITS.at(iv)
//                          + ")";
//
//            axisHeader += "}";
//
//            if(legend)
//                axisHeader += getLegendAxisAttributes();
//
//            axisHeader += "]"; // close axis environment attributes
//            return axisHeader;
//        }

        string getLegendAxisAttributes() {

            string legend = string(",\nlegend columns=3,\n")
                + "legend to name=planespannerlegend,\n"
                + "legend entries={";
            for(auto name : ALGORITHM_NAMES)
                legend += name + ",";
            legend += "}";

            return legend;
        }
        string getAxisFooter() {
            string axisFooter("\n\n\\end{axis}\n\n");
            return axisFooter;
        }

    private:

    }; // class PgfplotsPrinter


//
//
//    void plotResults(const string& dist, const BoundedDegreeSpannerResultSet &results, LatexPrinter* addToPrinter) {
//
//        // Create plot names
//        vector<string> plotNames;
//        transform(PGFPLOT_NAMES.begin(),
//                  PGFPLOT_NAMES.end(),
//                  back_inserter(plotNames),
//                  [&dist](const string& str) {
//                      string plotName = str + " (" + dist + ")";
//                      //cout<<plotName;
//                      return plotName;
//                  });
//
//        const double scale = 1;
//        auto plotNameIt = plotNames.begin();
//        bool first = true;
//        for( unsigned iv=0;iv<IV_NAMES.size();++iv) {
//        //for( const auto& iv : IV_NAMES ) {
//            //assert(plotNameIt != plotNames.end());
//            string caption = *plotNameIt++;
//            string filename = string("exp-plot-") + caption;
//            boost::erase_all(filename, " ");
//            boost::erase_all(filename, ")");
//            boost::erase_all(filename, "(");
//
//            PgfplotPrinter singlePlotter(addToPrinter->m_directory, filename);
//            //singlePlotter.setCaption(caption);
//            singlePlotter.plotAxis(iv,results,first);
//            first = false;
//
//            addToPrinter->addToDocument(singlePlotter,PRECOMPILE_SUBDOCUMENTS);
//        }
//    }

} // namespace spanners


#endif //SPANNERS_PGFPLOTSPRINTER_H
