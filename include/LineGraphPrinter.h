#ifndef GEOMETRIC_SPANNERS_PGFPLOTSPRINTER_H
#define GEOMETRIC_SPANNERS_PGFPLOTSPRINTER_H

/*
 * TODO: Separate statistics calculations into a separate
 *  class that will be accepted by PgfplotsPrinter and
 *  TablePrinter to print appropriate results without the
 *  need to perform redundant measurements on the data
 */

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

    class PgfplotsPrinter : public TikzPrinter {
    public:
        //typedef map< Algorithm, map< size_t, vector<BoundedDegreeSpannerResult>>> ResultMap;

        // Palette generated using https://coolors.co/
        vector<string> Colors = {
                "f94144","f3722c","f8961e",
                "f9844a","f9c74f","90be6d",
                "43aa8b","4d908e","577590",
                "277da1","360568"
        };


        PgfplotsPrinter(string filename, string documentType = "standalone")
                : TikzPrinter(filename,documentType){
            m_body = Body(getTikzHeader(), "", getTikzFooter());
        }




        string getLegendAxisAttributes() {
            string legend = string("legend columns=1,\n")
                + "legend to name=planespannerlegend\n";

            return legend;
        }
        string getLegend() {
            return "\\ref{planespannerlegend}\n\n";
        }
        void plotResults(const BoundedDegreeSpannerResultSet::ReducedResultMap &results) {
            bool addLegend = true;
            const double scale = 1;
            string plotHeader = string("")
                    + "\\begin{tikzpicture}\n\n"
                    + "\\begin{axis}[grid=major,xlabel=$n$, "
                    + "xtick={";

            string xTicks("");
            for( const auto& level : results.begin()->second ){
                xTicks += to_string(level.first);
                xTicks += ",";
            }
            xTicks = xTicks.substr( 0, xTicks.size()-1 );

            plotHeader += xTicks
                      + "\n}, ylabel near ticks,ylabel={IV}";


            vector<string> body(3, plotHeader);
            size_t colorIndex = 0;

            body[0] += string(",") + getLegendAxisAttributes();

            for( auto& graph : body ) {
                graph += "]";
            }

            for( auto alg : results ) {
                for( auto& graph : body) {
                    graph += "\n\n\\addplot";
                    graph += "[thick,color="
                            + Colors[colorIndex]
                            + ",mark=*]";
                    graph += " coordinates {\n";
                }
                for( auto level : alg.second ) {
                    body[0] += "("
                              + to_string(level.first)
                              + ","
                              + to_string(level.second.runtime)
                              + ")\n";

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
                for( auto& graph : body ) {
                    graph += string("}node [pos=1.15, above left] {};\n\n");
                }
                body[0] += "\\addlegendentry{\\textsc{\\tiny{"
                    + Names[alg.first]
                    + "}}}";
                ++colorIndex;
            }
            for( auto& graph : body ) {
                graph += string("\n\n\\end{axis}\n\n");
                graph += "\\end{tikzpicture}\n\n";
            }
            m_body.content += std::accumulate( body.begin(), body.end(), string("") );

        }

    private:
        //ResultMap _results;

    }; // class PgfplotsPrinter



} // namespace unf_spanners


#endif //GEOMETRIC_SPANNERS_PGFPLOTSPRINTER_H
