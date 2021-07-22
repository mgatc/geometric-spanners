//
// Created by matt on 7/22/21.
//

#ifndef GEOMETRIC_SPANNERS_TABLEPRINTER_H
#define GEOMETRIC_SPANNERS_TABLEPRINTER_H

#include "LatexPrinter.h"

namespace unf_spanners {

    class TablePrinter : public LatexPrinter {
    public:

        TablePrinter(string filename, string documentType = "standalone")
        : LatexPrinter(filename,documentType){
            m_body = Body{getTableHeader(IV_NICE_NAMES), "", getTableFooter()};
        }
        void tabulateResults(const BoundedDegreeSpannerResultSet &results) {
            string table;
            for( const auto& alg : results.m_reducedLevels ) {
                table += "\\textsc{"
                        + ALGORITHM_NAMES[alg.first]
                        + "}";
                for( auto iv : IV_NAMES){
                    table += " & $" + unf_spanners::to_string(alg.second.IV.at(iv)) + "$";
                }
//                for( auto iv : alg.second.IV ) {
//                    table += " & $" + unf_spanners::to_string(iv.second) + "$";
//                }
                table += " \\\\\n\n";
            }
            m_body.content = table;
        }
        string getTableHeader(const vector<string>& colHeaders) {
            std::string tableHeader = "\\begin{tabular}{|";
            unsigned tableCols = colHeaders.size() + 1;
            for( unsigned i=0; i<tableCols; ++i ){
                tableHeader += "c|";
            }
            tableHeader += "}\n";
            tableHeader += "\\hline\n\n";
            tableHeader += "\\textbf{Algorithm} & ";

            for( auto attr : colHeaders ){
                tableHeader += "\\textbf{" + attr + "} & ";
            }
            tableHeader.erase( std::prev(tableHeader.end(),2) );
            tableHeader += "\\\\";

            tableHeader += "\n\n\\hline\n";
            return tableHeader;
        }
        string getTableFooter() {
            std::string tableFooter = std::string("")
                                      + "\\hline\n"
                                      + "\\end{tabular}\n\n";
            return tableFooter;
        }
    private:

    };

}

#endif //GEOMETRIC_SPANNERS_TABLEPRINTER_H
