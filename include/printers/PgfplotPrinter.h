#ifndef SPANNERS_PGFPLOTSPRINTER_H
#define SPANNERS_PGFPLOTSPRINTER_H

#include <map>
#include <iostream>
#include <optional>
#include <vector>

#include "printers/GraphPrinter.h"
#include "tools/Results.h"
#include "tools/Utilities.h"

namespace spanners {

    using namespace std;

    class PgfplotPrinter : public TikzPrinter {
    public:
        PgfplotPrinter(string directory, string filename, string documentType = "standalone")
                : TikzPrinter(directory,filename,documentType){
            m_body = Body{getTikzHeader("scale=0.55"),
                          "",
                          getTikzFooter() };
            for(auto color : Colors){
                defineColor(color);
            }
        }
        template<class PlotMap>
        void plotAxis(const string& iv, const PlotMap& results, const double xScale = 1.0, const string xScaleUnit= "", const bool isFirst = true) {

            assert(abs(xScale) > 0);

            if(m_algorithmMarkers.empty())
                for(const auto& spannerName : ALGORITHM_NAMES) // provides the ordering we want
                    addNewMarker(spannerName);


            // build axis header
            string allPlotsOfAxis = getAxisHeader(iv, results, xScale, xScaleUnit, isFirst);

            // build the plot
            for( const auto& name : ALGORITHM_NAMES ) {
                if( contains(results,name) ) {
                    const auto &spanner = results.at(name);
                    string ivPlot = getPlotHeader(name);//
                    for (const auto &level: spanner) {
                        auto xValue = static_cast<double>(level.first) / xScale;
                        auto yValue = level.second.template getIV<double>(iv);
                        string entry = "("
                                       + std::to_string(xValue)
                                       + ","
                                       + std::to_string(yValue)
                                       + ")\n";
                        ivPlot += entry;
                    }
                    ivPlot += getPlotFooter();
                    if (isFirst) {
                        ivPlot += getLegendEntry(name);
                    }
                    allPlotsOfAxis += ivPlot;
                }
            }
            allPlotsOfAxis += getAxisFooter(m_filename);
            m_body.content = allPlotsOfAxis;
        }
        string getLegend() {
            string refOpener("\\ref{");
            string refText = removeSpaces(m_caption) + "-legend";
            string refCloser("}");

            return refOpener + refText + refCloser;
        }
        string getLegendEntry(const string& legendText) {
            string legendEntry = "\\addlegendentry{\\texttt{"
                                 + legendText
                                 + "}}";
            return legendEntry;
        }
        string getPlotHeader(const string& algorithm) {
            string plotHeader = "\n\n\\addplot[";

            auto color = getColor();

            plotHeader += string("solid,")
                        + getMarkerText(algorithm);

            plotHeader += "]";
            plotHeader += " coordinates {\n";
            return plotHeader;
        }
        string getPlotFooter() {
            string plotFooter("}node [pos=1.15, above left] {};\n\n");
            return plotFooter;
        }
        template<class ResultMap>
        string getAxisHeader(const string& iv, const ResultMap& results, const double xScale = 1.0,
                             const string xScaleUnit = "", bool first = true) {

            string axisHeader = string("")
                                + "\\begin{axis}[";
//            if(!m_caption.empty()) axisHeader += "title={" + m_caption + "},";
            string legendText = "legend to name=" + removeSpaces(m_caption)
                    + "-legend, legend columns={3}, ";
            axisHeader += "yticklabel style={rotate=90,anchor=base,yshift=0.2cm}, ";
            axisHeader += string("scaled ticks=false,grid=major,xlabel={$n$")
                          + (xScaleUnit.empty() ? "" : string(" (in ") + xScaleUnit + ")}, ")
                          + (first ? legendText : "")
                          + "xtick={";
            string xTicks = "";

            for( auto level : results.begin()->second ){
                xTicks += std::to_string(static_cast<double>(level.first) / xScale);
                xTicks += ",";
            }
            xTicks = xTicks.substr( 0, xTicks.size()-1 );


            auto firstTick = results.begin()->second.begin();
            auto secondTick = std::next(firstTick);
            auto lastTick = std::prev(results.begin()->second.end());

            auto firstTickScaled = static_cast<double>(firstTick->first) / xScale;
            auto secondTickScaled = static_cast<double>(secondTick->first) / xScale;
            auto lastTickScaled = static_cast<double>(lastTick->first) / xScale;

            double margin = std::abs(firstTickScaled - secondTickScaled);
            double xMin = firstTickScaled - margin;
            double xMax = lastTickScaled + margin;

            axisHeader += xTicks
                          + "\n}"
                          + ", ylabel near ticks,ylabel={"
                          + m_ivNiceNames.at(iv)
                          + "}, legend pos=north west,"
                          + "xmin="
                          + std::to_string(xMin)
                          + ",xmax="
                          + std::to_string(xMax);

            axisHeader += "]"; // close axis environment attributes
            return axisHeader;
        }
        string getAxisFooter(const string& plotLabel = "") {
            string axisFooter("\n\n\\end{axis}\n\n");
            axisFooter += "\\label{plots:"
                          + plotLabel
                          + "}";
            return axisFooter;
        }
    private:
        static vector<string> Marks;
        static size_t m_markIndex; // a valid index of Marks
        static string getMark() {
            auto currentMark = m_markIndex;
            m_markIndex = (m_markIndex+1) % Marks.size();
            return Marks.at(currentMark);
        }

        static vector<string> Colors;
        static size_t m_colorIndex;
        static string getColor() {
            auto currentColor = m_colorIndex;
            m_colorIndex = (m_colorIndex+1) % Colors.size();
            return Colors.at(currentColor);
        }
        static void addNewMarker(const string& algorithm) {
            string color = getColor();
            string markerText;
            markerText += "mark=" + getMark() + ","
                    + "mark color=" + color + ","
                    + "color=" + color;
            m_algorithmMarkers[algorithm] = markerText;
        }
        static map<string,string> m_algorithmMarkers;
        static string getMarkerText(const string& algorithm) {
            if(m_algorithmMarkers.find(algorithm) == m_algorithmMarkers.end())
                addNewMarker(algorithm);
            return m_algorithmMarkers.at(algorithm);
        }


        static map<string,string> m_ivNiceNames;
    }; // class PgfplotsPrinter

    map<string,string> PgfplotPrinter::m_algorithmMarkers;
    vector<string> PgfplotPrinter::Marks = {
            //"otimes*",
            "o", "triangle", "pentagon", "square",  "diamond"
    };
    // Palette generated using https://coolors.co/
    vector<string> PgfplotPrinter::Colors = {
            "22dddd",
            "2222dd",
            "dd22dd",
//            "264653",
//            "1789a6",
//            "e76f51",
//            "33c1b1",
            //"287271",
//            "2a9d8f",
//            "e9c46a",
            //"f4a261",

    };
    size_t PgfplotPrinter::m_markIndex = 0;
    size_t PgfplotPrinter::m_colorIndex = 1;

    map<string,string> PgfplotPrinter::m_ivNiceNames = {
            {"runtime",             "Average execution time (s)"},
            {"degree",              "Maximum observed degree"},
            {"degreeAvg",           "Average observed degree per spanner"},
            {"avgDegreePerPoint",   "Average observed degree per point"},
            {"maxStretchFactor",    "Maximum observed stretch-factor"},
            {"avgStretchFactor",    "Average observed stretch-factor"},
            {"lightness",           "Average lightness"},
    };


} // namespace spanners


#endif //SPANNERS_PGFPLOTSPRINTER_H
