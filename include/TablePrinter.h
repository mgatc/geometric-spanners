//
// Created by matt on 7/22/21.
//

#ifndef GEOMETRIC_SPANNERS_TABLEPRINTER_H
#define GEOMETRIC_SPANNERS_TABLEPRINTER_H

#include <iomanip>
#include <sstream>

#include "LatexPrinter.h"
#include "utilities.h"

namespace unf_spanners {

    class TablePrinter : public LatexPrinter {
    public:

        TablePrinter(string filename, string documentType = "standalone")
        : LatexPrinter(filename,documentType){
            m_body = Body{getTableHeader(IV_SYMBOLS), "", getTableFooter()};
        }
        void ignoreIV(unsigned iv) {
            m_ignore.insert(iv);
            m_body.header = getTableHeader(IV_SYMBOLS);
        }
        void addColumn(string header, const vector<string>& values) {
            m_added.emplace(header,values);
            m_body.header = getTableHeader(IV_SYMBOLS);
        }
        void tabulateResults(const BoundedDegreeSpannerResultSet &results) {
            string table;
            for( const auto& alg : results.m_reducedLevels ) {
                table += "\\textsc{"
                        + ALGORITHM_NAMES[alg.first]
                        + "}";
                for( auto column : m_added ) {
                    table += " & $"
                            + column.second.at(alg.first)
                            + "$";
                }
                for( unsigned iv=0; iv<IV_SYMBOLS.size();++iv ) {//auto attr : colHeaders ){
                    if (!contains(m_ignore, iv)) {
                        ostringstream stream;
                        stream << fixed;
                        stream << setprecision(IV_PRECISION[iv]);
                        stream << alg.second.IV.at(IV_NAMES[iv]); // TODO: this is a variant, whose << overload uses to_string before setprecision can affect the output
                        string formattedNum = stream.str();
                        table += " & $" + formattedNum + "$";
                    }
                }
                table += " \\\\\n\n";
            }
            m_body.content = table;
        }
        string getTableHeader(const vector<string>& colHeaders) {
            std::string tableHeader = "\\begin{tabular}{|";
            unsigned tableCols = colHeaders.size() + 1;
            tableCols -= m_ignore.size();
            tableCols += m_added.size();
            for( unsigned i=0; i<tableCols; ++i ){
                tableHeader += "c|";
            }
            tableHeader += "}\n";
            tableHeader += "\\hline\n\n";
            tableHeader += "\\textbf{Algorithm}";

            for( auto column : m_added ) {
                tableHeader += " & \\textbf{" + column.first + "}";
            }

            for( unsigned iv=0; iv<colHeaders.size();++iv ) //auto attr : colHeaders ){
                if( !contains(m_ignore,iv) )
                    tableHeader += " & \\textbf{" + colHeaders[iv] + "}";

            //tableHeader.erase( std::prev(tableHeader.end(),2) );
            tableHeader += "\\\\";

            tableHeader += "\n\n\\hline\n";
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
        set<unsigned> m_ignore;
        map<string,vector<string>> m_added;
    };

}

#endif //GEOMETRIC_SPANNERS_TABLEPRINTER_H
