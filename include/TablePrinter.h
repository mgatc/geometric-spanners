#ifndef PLANESPANNERS_TABLEPRINTER_H
#define PLANESPANNERS_TABLEPRINTER_H

#include <iomanip>
#include <sstream>
#include <string>

#include "LatexPrinter.h"
#include "utilities.h"

namespace planespanners {

    const string TABLE_COLOR_1 = "ffffff";
    const string TABLE_COLOR_2 = "cccccc";


    class TablePrinter : public LatexPrinter {
    public:

        enum CellHighlightStyle {
            None,
            MaxInColumn,
            MaxInRow
        };

        TablePrinter(string filename, string documentType = "standalone")
        : LatexPrinter(filename,documentType){
            defineColor(TABLE_COLOR_1);
            defineColor(TABLE_COLOR_2);
        }
        void ignoreIV(unsigned iv) {
            m_ignore.insert(iv);
        }
        void addColumn(string header, const vector<string>& values, int priority = -1) {
            if( priority < 0 ) {
                priority = 0;
                while (contains(m_added, priority))
                    ++priority;
            }
            m_added.emplace(priority,make_pair(header,values));
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
            optional<size_t> highlightRow, highlightCol;
            optional<vector<string>> highlightCols;
            if(highlightStyle==CellHighlightStyle::MaxInColumn) {
                transform(m_added.begin(),m_added.end(),back_inserter(*highlightCols),
                    [](const auto& column) {
                        return *min_element(column.second.second.begin(),column.second.second.end(),
                            []( const auto& lhs, const auto& rhs ) {
                                return stod(lhs) < stod(lhs);
                        });
                });
            }
            for(size_t row=0; row<numRows; ++row) {
                optional<string> highlightValue;
                if(highlightStyle==CellHighlightStyle::MaxInRow) {
                    highlightValue = make_optional(min_element(next(m_added.begin()),m_added.end(),
                  [row](auto &column1, auto &column2){
                            return stod(column1.second.second.at(row)) < stod(column2.second.second.at(row));
                        })->second.second.at(row));
                }
                size_t columnNum = 0;
                for(const auto& cell : m_added ) { //auto attr : colHeaders )
                    string cellValue = cell.second.second.at(row);
                    string prefix = "$";
                    string suffix = "$ &";
                    if( (highlightValue && *highlightValue == cellValue )
                      ||(highlightCols && (*highlightCols)[columnNum] == cellValue )) {
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
                tableHeader += "$\\textbf{" + column.second.first + "}$ &";
            }

            tableHeader = tableHeader.substr( 0, tableHeader.size()-1 ); // remove trailing character (&)
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
    protected:
        set<unsigned> m_ignore;
        map<int,pair<string,vector<string>>> m_added;
    };

    void tabulateSummaryResults(const DistributionType dist, const BoundedDegreeSpannerResultSet &results, LatexPrinter* addToPrinter) {

        string distributionNameNoSpaces = DISTRIBUTION_NAMES.at(dist);
        boost::erase_all(distributionNameNoSpaces, " ");
        string tableFilename = string("exp-table")
                               + distributionNameNoSpaces;

//                TablePrinter table(tableFilename);
//
        string tableCaption = DISTRIBUTION_NAMES.at(dist) + " Summary";
//                table.setCaption(tableCaption);
        string filename = string("exp-table_summary-") + tableCaption;
        boost::erase_all(filename, " ");
        boost::erase_all(filename, ")");
        boost::erase_all(filename, "(");

        TablePrinter table(filename);
        table.addColumn( ALGORITHM_SYMBOL, ALGORITHM_NAMES, 0 );
        table.addColumn(DEGREE_BOUND_SYMBOL, DEGREE_BOUND_PER_ALGORITHM, 1);
        table.addColumn(STRETCH_FACTOR_BOUND_SYMBOL, STRETCH_FACTOR_BOUND_PER_ALGORITHM, 4);

        size_t col = 0;

        // skip runtime for summary table by starting at 1
        for( int iv=1; iv<(int)IV_SYMBOLS.size(); ++iv) { //auto attr : colHeaders ){
            //if( iv == 0 ) continue;

            vector<string> singleColumn;
            for(auto row : results.m_reducedLevels) {
                ostringstream stream;
                stream << fixed;
                stream << setprecision(IV_PRECISION.at(iv));
                stream << row.second.IV.at(IV_NAMES.at(iv));
                singleColumn.push_back(stream.str());
            }
            table.addColumn( IV_SYMBOLS.at(iv),singleColumn);
        }
        table.tabulate(TablePrinter::CellHighlightStyle::MaxInColumn);
        addToPrinter->addToDocument(table, PRECOMPILE_SUBDOCUMENTS);
    }


    void getSingleIVColumn(size_t iv, Algorithm alg, const BoundedDegreeSpannerResultSet &results,  vector<string>& singleColumn) {
        for( const auto& level : results.m_reducedSamples.at(alg) ) {
            ostringstream stream;
            stream << fixed;
            stream << setprecision(IV_PRECISION.at(iv));
            stream << level.second.IV.at(IV_NAMES.at(iv));
            singleColumn.push_back(stream.str());
        }
    }

    void tabulateIVs(const DistributionType dist, const BoundedDegreeSpannerResultSet &results, LatexPrinter* addToPrinter) {

        // Create table names
        vector<string> tableNames;
        transform(PGFPLOT_NAMES.begin(),
                  PGFPLOT_NAMES.end(),
                  back_inserter(tableNames),
                  [&dist](const string& str) {
                      string tableName = str + " (" + DISTRIBUTION_NAMES.at(dist) + ")";
                      //cout<<tableName;
                      return tableName;
                  });

        // get levels of n
        vector<string> nLevels;
        //for( const auto& alg :  ) {
            for( const auto& level : results.m_reducedSamples.begin()->second ) {
                nLevels.push_back(std::to_string(level.first));
            }
        //}

        auto tableNameIt = tableNames.begin();
        //bool first = true;
        for( unsigned iv=0;iv<IV_NAMES.size();++iv) {
            //for( const auto& iv : IV_NAMES ) {
            //assert(tableNameIt != plotNames.end());
            string caption = *tableNameIt++;
            string filename = string("exp-table_iv-") + caption;
            boost::erase_all(filename, " ");
            boost::erase_all(filename, ")");
            boost::erase_all(filename, "(");

            TablePrinter singleTabulater(filename);
            //singleTabulater.setCaption(caption);
            size_t i=0;
            singleTabulater.addColumn(N_SYMBOL, nLevels, i++);

            for(int alg=Algorithm::AlgorithmFirst;
                alg!=Algorithm::AlgorithmLast; ++alg ) {
                vector<string> singleColumn;
                getSingleIVColumn( iv, Algorithm(alg), results, singleColumn);
                singleTabulater.addColumn(ALGORITHM_NAMES[alg], singleColumn, i++);
            }
            singleTabulater.tabulate(TablePrinter::CellHighlightStyle::MaxInRow);
            addToPrinter->addToDocument(singleTabulater, PRECOMPILE_SUBDOCUMENTS);
        }
    }

}

#endif //PLANESPANNERS_TABLEPRINTER_H
