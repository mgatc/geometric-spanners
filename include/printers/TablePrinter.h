#ifndef SPANNERS_TABLEPRINTER_H
#define SPANNERS_TABLEPRINTER_H

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "printers/LatexPrinter.h"
#include "tools/Utilities.h"

namespace spanners {

    using namespace std;

    const string TABLE_COLOR_1 = "ffffff";
    const string TABLE_COLOR_2 = "cccccc";


    class TablePrinter : public LatexPrinter {
    public:

        enum CellHighlightStyle {
            None,
            MaxInColumn,
            MaxInRow
        };

        TablePrinter(string directory, string filename, string documentType = "standalone")
        : LatexPrinter(directory,filename,documentType){
            defineColor(TABLE_COLOR_1);
            defineColor(TABLE_COLOR_2);
        }
        void ignoreIV(unsigned iv) {
            m_ignore.insert(iv);
        }
        void addColumn(string header, const vector<string>& values, int precision = -1, int priority = -1) {
            auto valuesCopy(values);
            SetPrecision precisionSetter{precision};
            transform(valuesCopy.begin(),valuesCopy.end(),valuesCopy.begin(),precisionSetter);
            if( priority < 0 ) {
                priority = 0;
                while (contains(m_added, priority))
                    ++priority;
            }
            m_added.emplace(priority,make_pair(header,valuesCopy));
        }
        void tabulate(TablePrinter::CellHighlightStyle highlightStyle = TablePrinter::CellHighlightStyle::None) {
            m_body = Body{getTableHeader(),
                          getTableBody(highlightStyle),
                          getTableFooter() };
        }
        string getTableBody(TablePrinter::CellHighlightStyle highlightStyle) {
            assert(!m_added.empty());
            typedef TablePrinter::CellHighlightStyle CellHighlightStyle;

            string table;
            size_t numRows = m_added.begin()->second.second.size();

            vector<string> highlightCols;
            if(highlightStyle==CellHighlightStyle::MaxInColumn) {
                // first column has algorithm names, skip
                transform(next(m_added.begin()),m_added.end(),back_inserter(highlightCols),
                    [](const auto& column) -> string {

                        return *min_element(column.second.second.cbegin(),column.second.second.cend(),
                            []( const auto& lhs, const auto& rhs ) -> bool {
                                return stod(lhs) < stod(rhs);
                        });
                });
            }
            for(size_t row=0; row<numRows; ++row) {
                optional<string> highlightValue;
                if(highlightStyle==CellHighlightStyle::MaxInRow) {
                    // first column has n levels, skip
                    highlightValue = make_optional(min_element(next(m_added.begin()),m_added.end(),
                  [row](auto &column1, auto &column2){
                            return stod(column1.second.second.at(row)) < stod(column2.second.second.at(row));
                        })->second.second.at(row));
                }
                int columnNum = -1; // 0 is the first index of the column highlight-data but we want to skip the first col (names)
                for(const auto& cell : m_added ) { //auto attr : colHeaders )
                    string cellValue = cell.second.second.at(row);
                    string prefix = "$";
                    string suffix = "$ &";
                    if( (highlightValue && *highlightValue == cellValue )
                      ||(!highlightCols.empty() && highlightCols[columnNum] == cellValue )) {
                        prefix += "\\textbf{";
                        suffix = "}" + suffix;
                    }
                    table += prefix + cellValue + suffix;

                    ++columnNum;
                }
                table = table.substr(0, table.size() - 1 ); // remove trailing character (&)
                table += "\\\\\n\n";
            }
            return table;
        }
        string getTableHeader() {

            std::string tableHeader = "\\rowcolors{1}{"
                                      + TABLE_COLOR_1 + "}{"
                                      + TABLE_COLOR_2 + "}\n\n"
                                      + "\\begin{tabular}{|";
            unsigned tableCols = m_added.size();
            tableCols -= m_ignore.size();
            assert(m_added.size()>m_ignore.size());
            for( unsigned i=0; i<tableCols; ++i ){
                tableHeader += "c|";
            }
            tableHeader += "}\n";
            tableHeader += "\\hline\n\n";


            //for(int cell = 0; cell < (int) m_added.size(); ++cell) { //auto attr : colHeaders )
            for( auto column : m_added ){
                tableHeader += column.second.first + " & ";
            }

            tableHeader = tableHeader.substr( 0, tableHeader.size()-2 ); // remove trailing character (&)
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
        static map<string,string> m_ivNiceNames;
    protected:
        set<unsigned> m_ignore;
        map<int,pair<string,vector<string>>> m_added;

    };

    map<string,string> TablePrinter::m_ivNiceNames = {
            {"runtime",             "$\\mathrm{runtime (s)}$"},
            {"degree",              "$\\Delta$"},
            {"degreeAvg",           "$\\Delta_\\mathrm{avg}$"},
            {"avgDegreePerPoint",   "$\\Delta_\\mathrm{point}$"},
            {"maxStretchFactor",    "$t_\\mathrm{max}$"},
            {"avgStretchFactor",    "$t_\\mathrm{avg}$"},
            {"lightness",           "$\\ell$"},
    };

}

#endif //SPANNERS_TABLEPRINTER_H
