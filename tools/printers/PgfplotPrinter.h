#ifndef SPANNERS_PGFPLOTSPRINTER_H
#define SPANNERS_PGFPLOTSPRINTER_H

#include <map>
#include <iostream>
#include <optional>
#include <vector>

#include "GraphPrinter.h"
#include "../Results.h"
#include "../Utilities.h"

namespace bdps_experiment {

    class PgfplotPrinter : public TikzPrinter {
    public:
        PgfplotPrinter(std::string directory, std::string filename, std::string documentType = "standalone")
                : TikzPrinter(directory,filename,documentType){
            m_body = Body{getTikzHeader("scale=0.55"),
                          "",
                          getTikzFooter() };
            for(auto color : Colors){
                defineColor(color);
            }
        }
        template<class PlotMap>
        void plotAxis(const std::string& iv, const PlotMap& results, const double xScale = 1.0, const std::string xScaleUnit= "", const bool isFirst = true) {

            assert(abs(xScale) > 0);

            if(m_algorithmMarkers.empty())
                for(const auto& spannerName : spanner::bdps::ALGORITHM_NAMES) // provides the ordering we want
                    addNewMarker(spannerName);


            // build axis header
            std::string allPlotsOfAxis = getAxisHeader(iv, results, xScale, xScaleUnit, isFirst);

            // build the plot
            for( const auto& name : spanner::bdps::ALGORITHM_NAMES ) {
                if( contains(results,name) ) {
                    const auto &spanner = results.at(name);
                    std::string ivPlot = getPlotHeader(name);//
                    for (const auto &level: spanner) {
                        auto xValue = static_cast<double>(level.first) / xScale;
                        auto yValue = level.second.template getIV<double>(iv);
                        std::string entry = "("
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
        std::string getLegend() {
            std::string refOpener("\\ref{");
            std::string refText = removeSpaces(m_caption) + "-legend";
            std::string refCloser("}");

            return refOpener + refText + refCloser;
        }
        std::string getLegendEntry(const std::string& legendText) {
            std::string legendEntry = "\\addlegendentry{\\texttt{"
                                 + legendText
                                 + "}}";
            return legendEntry;
        }
        std::string getPlotHeader(const std::string& algorithm) {
            std::string plotHeader = "\n\n\\addplot[";

            auto color = getColor();

            plotHeader += string("solid,")
                        + getMarkerText(algorithm);

            plotHeader += "]";
            plotHeader += " coordinates {\n";
            return plotHeader;
        }
        std::string getPlotFooter() {
            std::string plotFooter("}node [pos=1.15, above left] {};\n\n");
            return plotFooter;
        }
        template<class ResultMap>
        std::string getAxisHeader(const std::string& iv, const ResultMap& results, const double xScale = 1.0,
                             const std::string xScaleUnit = "", bool first = true) {

            std::string axisHeader = std::string("")
                                + "\\begin{axis}[";
//            if(!m_caption.empty()) axisHeader += "title={" + m_caption + "},";
            std::string legendText = "legend to name=" + removeSpaces(m_caption)
                    + "-legend, legend columns={3}, ";
            axisHeader += "yticklabel style={rotate=90,anchor=base,yshift=0.2cm}, ";
            axisHeader += std::string("scaled ticks=false,grid=none,xlabel={$n$")
                          + (xScaleUnit.empty() ? "" : std::string(" (in ") + xScaleUnit + ")}, ")
                          + (first ? legendText : "")
                          + "xtick={";
            std::string xTicks = "";

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
        std::string getAxisFooter(const std::string& plotLabel = "") {
            std::string axisFooter("\n\n\\end{axis}\n\n");
            axisFooter += "\\label{plots:"
                          + plotLabel
                          + "}";
            return axisFooter;
        }
    private:
        static std::vector<std::string> MarkStyles;
        static std::vector<std::string> Marks;
        static size_t m_markIndex; // a valid index of Marks
        static std::string getMark() {
            auto currentMark = m_markIndex;
            m_markIndex = (m_markIndex+1) % Marks.size();
            return Marks.at(currentMark);
        }

        static std::vector<std::string> Colors;
        static size_t m_colorIndex;
        static std::string getColor() {
            auto currentColor = m_colorIndex;
            m_colorIndex = (m_colorIndex+1) % Colors.size();
            return Colors.at(currentColor);
        }
        static int m_markStyleCurrent;
        static std::string getMarkStyle() {
            std::string currentStyle = MarkStyles[m_markStyleCurrent];
            m_markStyleCurrent = (m_markStyleCurrent+1) % MarkStyles.size();
            return currentStyle;
        }
        static void addNewMarker(const std::string& algorithm) {
//            string color = getColor();
            std::string markerText;
//            markerText += "mark=" + getMark() + ","
//                    + "mark color=" + color + ","
//                    + "color=" + color;
            markerText = getMarkStyle();
            m_algorithmMarkers[algorithm] = markerText;
        }
        static std::map<std::string,std::string> m_algorithmMarkers;
        static std::string getMarkerText(const std::string& algorithm) {
            if(m_algorithmMarkers.find(algorithm) == m_algorithmMarkers.end())
                addNewMarker(algorithm);
            return m_algorithmMarkers.at(algorithm);
        }


        static std::map<std::string,std::string> m_ivNiceNames;
    }; // class PgfplotsPrinter

    int PgfplotPrinter::m_markStyleCurrent = 0;
    std::map<std::string,std::string> PgfplotPrinter::m_algorithmMarkers;
    std::vector<std::string> PgfplotPrinter::MarkStyles = {
            "color=black,mark options={fill=black},mark=square*",
            "color=black,mark options={fill=black},mark=pentagon*",
            "color=black,mark options={fill=black},mark=diamond*",
            "color=black,mark options={fill=black},mark=*",
            "color=black,mark options={fill=black},mark=triangle*",
            "color=black,mark=square",
            "color=black,mark=pentagon",
            "color=black,mark=diamond",
            "color=black,mark=o",
            "color=black,mark=otimes",
            "color=black,mark=triangle",
            "color=black,mark=oplus",
    };
    std::vector<std::string> PgfplotPrinter::Marks = {
            "otimes*",
            "oplus*",
            "o", "triangle", "pentagon", "square",  "diamond"
    };
    // Palette generated using https://coolors.co/
    std::vector<std::string> PgfplotPrinter::Colors = {
            "000000",
            "2288DD"
//            "22dddd",
//            "2222dd",
//            "dd22dd",
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

    std::map<std::string,std::string> PgfplotPrinter::m_ivNiceNames = {
            {"runtime",             "Average execution time (s)"},
            {"degree",              "Maximum observed degree"},
            {"degreeAvg",           "Average observed degree per spanner"},
            {"avgDegreePerPoint",   "Average observed degree per point"},
            {"maxStretchFactor",    "Maximum observed stretch-factor"},
            {"avgStretchFactor",    "Average observed stretch-factor"},
            {"lightness",           "Average lightness"},
    };


} // namespace bdps_experiment


#endif //SPANNERS_PGFPLOTSPRINTER_H
