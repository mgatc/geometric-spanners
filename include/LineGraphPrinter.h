#ifndef GEOMETRIC_SPANNERS_PGFPLOTSPRINTER_H
#define GEOMETRIC_SPANNERS_PGFPLOTSPRINTER_H

#include <map>
#include <iostream>
#include <optional>
#include <vector>

#include "GeometricSpannerPrinter.h"
#include "names.h"
#include "Result.h"
#include "utilities.h"

namespace unf_spanners {

    using namespace std;

    void plotResults(const BoundedDegreeSpannerResultSet &results, LatexPrinter* addToPrinter);

    class PgfplotsPrinter : public TikzPrinter {
    public:
        // Palette generated using https://coolors.co/
        vector<string> Colors = {
                "f94144","f3722c","f8961e",
                "f9844a","f9c74f","90be6d",
                "43aa8b","4d908e","577590",
                "277da1","360568"
        };


        PgfplotsPrinter(string filename, string documentType = "standalone")
                : TikzPrinter(filename,documentType){
            m_body = Body{getTikzHeader(),
                          "",
                          getTikzFooter() };
            for(auto color : Colors){
                defineColor(color);
            }
        }




//        string getLegendAxisAttributes() {
////            string legend = string(",\nlegend columns=3,\n")
////                + "legend to name=planespannerlegend\n";
//
//            return legend;
//        }
//        string getLegend() {
//            TikzPrinter tikz("templegend");
//            string legend = "\\begin{axis}[%\n"
//                            "            hide axis,\n"
//                            "                    xmin=10,\n"
//                            "                    xmax=50,\n"
//                            "                    ymin=0,\n"
//                            "                    ymax=0.4,\n"
//                            "                    legend style={draw=white!15!black,legend cell align=left},\n"
//                            "                    legend columns=3\n"
//                            "            ]\n\n";
//            for(string name : ALGORITHM_NAMES) {
//                legend += "\\addlegendentry{"
//                        + name
//                        +  "}";
//            }
//            legend += "\\end{axis}";
//
//            tikz.addRawText(legend);
//
//        }
        void plotAxis(unsigned iv, const BoundedDegreeSpannerResultSet &results) {

            // build axis header
            string allPlotsOfAxis = getAxisHeader(iv, results);

            for( const auto& alg : results.m_reducedSamples ) {
                // build the plot
                string plot = getPlotHeader(alg.first);

                for( const auto& level : alg.second ) {
                    // add the level, measurement pair
                    string entry = "("
                                   + std::to_string(level.first)
                                   + ","
                                   + unf_spanners::to_string(level.second.IV.at(IV_NAMES[iv]))
                                   + ")\n";
                    plot += entry;
                }
                plot += getPlotFooter(iv);

//                if(includeLegend)
//                    plot += addLegendEntry(alg.first);

                allPlotsOfAxis += plot;
            }
            allPlotsOfAxis += getAxisFooter();

            m_body.content = allPlotsOfAxis;
        }
        void plotResults(const BoundedDegreeSpannerResultSet &results) {
            unf_spanners::plotResults(results,this);
        }
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
        void addLegend() {
            m_body.footer += getLegendText(ALGORITHM_NAMES);
        }
        string getPlotHeader(Algorithm alg) {
            string plotHeader = "\n\n\\addplot";
            plotHeader += "[thick,color="
                    + Colors[alg]
                    + ",mark=*]";
            plotHeader += " coordinates {\n";
            return plotHeader;
        }
        string getPlotFooter(unsigned iv) {
            string plotFooter("}node [pos=1.15, above left] {};\n\n");
            plotFooter += "\\label{plots:"
                        + ALGORITHM_NAMES[iv]
                        + "}";
            return plotFooter;
        }
        string getAxisHeader(unsigned iv,
                             const BoundedDegreeSpannerResultSet &results) {
            string axisHeader = string("")
                                + "\\begin{axis}["
                                + "title={"
                                + PGFPLOT_NAMES[iv]
                                + "},"
                                + "grid=major,xlabel=$n$, "
                                + "xtick={";
            string xTicks;
            assert(!results.m_reducedSamples.empty());
            for( auto level : results.m_reducedSamples.begin()->second ){
                xTicks += std::to_string(level.first);
                xTicks += ",";
            }
            xTicks = xTicks.substr( 0, xTicks.size()-1 );

            axisHeader += xTicks
                          + "\n}, ylabel near ticks,ylabel={"
                          + IV_NICE_NAMES.at(iv)
                          + "}";

//            if(first)
//                axisHeader += getLegendAxisAttributes();

            axisHeader += "]"; // close axis environment attributes
            return axisHeader;
        }
        string getAxisFooter() {
            string axisFooter("\n\n\\end{axis}\n\n");
            return axisFooter;
        }

    private:

    }; // class PgfplotsPrinter




    void plotResults(const BoundedDegreeSpannerResultSet &results, LatexPrinter* addToPrinter, vector<string> plotNames = {}) {


        const double scale = 1;
        auto plotNameIt = plotNames.begin();
        bool first = false;
        for( unsigned iv=0;iv<IV_NAMES.size();++iv) {
        //for( const auto& iv : IV_NAMES ) {
            assert(plotNameIt != plotNames.end());
            string outputname = string("RPIS-")
                                + IV_NAMES[iv];
            PgfplotsPrinter singlePlotter(outputname);
            singlePlotter.plotAxis(iv,results);
            singlePlotter.setCaption(*plotNameIt++);
            if(first)singlePlotter.addLegend();
            first = false;

            addToPrinter->addToDocument(singlePlotter);
        }
    }

} // namespace unf_spanners


#endif //GEOMETRIC_SPANNERS_PGFPLOTSPRINTER_H
