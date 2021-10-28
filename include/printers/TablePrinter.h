#ifndef SPANNERS_TABLEPRINTER_H
#define SPANNERS_TABLEPRINTER_H

#include <iomanip>
#include <sstream>
#include <string>

#include "printers/LatexPrinter.h"
#include "tools/Utilities.h"

namespace spanners {

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





//    void tabulateSummaryResults(const string& dist, const BoundedDegreeSpannerResultSet &results, LatexPrinter* addToPrinter) {
//
//        string distributionNameNoSpaces(dist);
//        boost::erase_all(distributionNameNoSpaces, " ");
//        string tableFilename = string("exp-table")
//                               + distributionNameNoSpaces;
//
////                TablePrinter table(tableFilename);
////
//        string tableCaption = dist + " Summary";
////                table.setCaption(tableCaption);
//        string filename = string("exp-table_summary-") + tableCaption;
//        boost::erase_all(filename, " ");
//        boost::erase_all(filename, ")");
//        boost::erase_all(filename, "(");
//
//        TablePrinter table(addToPrinter->m_directory,filename);
//        table.addColumn( ALGORITHM_SYMBOL, ALGORITHM_NAMES, -1,0 );
//        table.addColumn(DEGREE_BOUND_SYMBOL, DEGREE_BOUND_PER_ALGORITHM,0,1);
//        table.addColumn(STRETCH_FACTOR_BOUND_SYMBOL, STRETCH_FACTOR_BOUND_PER_ALGORITHM,-1,4);
//
//        size_t col = 0;
//
//        // skip runtime for summary table by starting at 1
//        for( int iv=1; iv<(int)IV_SYMBOLS.size(); ++iv) { //auto attr : colHeaders ){
//            //if( iv == 0 ) continue;
//
//            vector<string> singleColumn;
//            for(auto row : results.m_reducedLevels) {
//                ostringstream stream;
//                stream << fixed;
//                stream << setprecision(IV_PRECISION.at(iv));
//                stream << row.second.IV.at(IV_NAMES.at(iv));
//                singleColumn.push_back(stream.str());
//            }
//            table.addColumn( IV_SYMBOLS.at(iv),singleColumn,6);
//        }
//        table.tabulate(TablePrinter::CellHighlightStyle::MaxInColumn);
//        addToPrinter->addToDocument(table, PRECOMPILE_SUBDOCUMENTS);
//    }


//    void getSingleIVColumn(size_t iv, BoundedDegreePlaneSpannerAlgorithm alg, const BoundedDegreeSpannerResultSet &results,  vector<string>& singleColumn) {
//        for( const auto& level : results.m_reducedSamples.at(alg) ) {
//            ostringstream stream;
//            stream << fixed;
//            stream << setprecision(IV_PRECISION.at(iv));
//            stream << level.second.IV.at(IV_NAMES.at(iv));
//            singleColumn.push_back(stream.str());
//        }
//    }
//
//    void tabulateIVsFromConfigExperiment(const string& dist, const BoundedDegreeSpannerResultSet &results, const map<index_t,string>& REAL_POINTSET_NAMES, LatexPrinter* addToPrinter) {
//        // Create table names
//        vector<string> tableNames;
//        transform(PGFPLOT_NAMES.begin(),
//                  PGFPLOT_NAMES.end(),
//                  back_inserter(tableNames),
//                  [&dist](const string& str) {
//                      string tableName = str + " (" + dist + ")";
//                      //cout<<tableName;
//                      return tableName;
//                  });
//
//        // get levels of n
//        vector<string> nLevels;
//        //for( const auto& alg :  ) {
//        for( const auto& level : results.m_reducedSamples.begin()->second ) {
//            nLevels.push_back(std::to_string(level.first));
//        }
//        //}
//        vector<string> names(REAL_POINTSET_NAMES.size());
//        transform(REAL_POINTSET_NAMES.begin(),REAL_POINTSET_NAMES.end(),names.begin(),
//            [](const auto& entry) -> string {
//                return "\\text{" + entry.second + "}";
//        });
//
//        auto tableNameIt = tableNames.begin();
//        //bool first = true;
//        for( unsigned iv=0;iv<IV_NAMES.size();++iv) {
//            //for( const auto& iv : IV_NAMES ) {
//            //assert(tableNameIt != plotNames.end());
//            string caption = *tableNameIt++;
//            string filename = string("exp-table_iv-") + caption;
//            boost::erase_all(filename, " ");
//            boost::erase_all(filename, ")");
//            boost::erase_all(filename, "(");
//
//            TablePrinter singleTabulater(addToPrinter->m_directory, filename);
//            //singleTabulater.setCaption(caption);
//            size_t i=0;
//            singleTabulater.addColumn("Point set", names, -1,i++);
//            singleTabulater.addColumn(N_SYMBOL, nLevels, 0,i++);
//
//            vector<string> AlgorithmNamesSmallText(ALGORITHM_NAMES);
////            transform(ALGORITHM_NAMES.begin(),ALGORITHM_NAMES.end(),back_inserter(AlgorithmNamesSmallText),
////                [](const auto& text){
////                    return "\\tiny{" + text + "}";
////            });
//
//            for(int alg=BoundedDegreePlaneSpannerAlgorithm::AlgorithmFirst;
//                alg!=BoundedDegreePlaneSpannerAlgorithm::AlgorithmLast; ++alg ) {
//                vector<string> singleColumn;
//                getSingleIVColumn( iv, BoundedDegreePlaneSpannerAlgorithm(alg), results, singleColumn);
//                singleTabulater.addColumn(AlgorithmNamesSmallText[alg], singleColumn, -1,i++);
//            }
//            singleTabulater.tabulate(TablePrinter::CellHighlightStyle::MaxInRow);
//            addToPrinter->addToDocument(singleTabulater, PRECOMPILE_SUBDOCUMENTS);
//        }
//    }
//
//
//
//
//    void tabulateIVs(const string& dist, const BoundedDegreeSpannerResultSet &results, LatexPrinter* addToPrinter) {
//
//        // Create table names
//        vector<string> tableNames;
//        transform(PGFPLOT_NAMES.begin(),
//                  PGFPLOT_NAMES.end(),
//                  back_inserter(tableNames),
//                  [&dist](const string& str) {
//                      string tableName = str + " (" + dist + ")";
//                      //cout<<tableName;
//                      return tableName;
//                  });
//
//        // get levels of n
//        vector<string> nLevels;
//        //for( const auto& alg :  ) {
//            for( const auto& level : results.m_reducedSamples.begin()->second ) {
//                nLevels.push_back(std::to_string(level.first));
//            }
//        //}
//
//        auto tableNameIt = tableNames.begin();
//        //bool first = true;
//        for( unsigned iv=0;iv<IV_NAMES.size();++iv) {
//            //for( const auto& iv : IV_NAMES ) {
//            //assert(tableNameIt != plotNames.end());
//            string caption = *tableNameIt++;
//            string filename = string("exp-table_iv-") + caption;
//            boost::erase_all(filename, " ");
//            boost::erase_all(filename, ")");
//            boost::erase_all(filename, "(");
//
//            TablePrinter singleTabulater(addToPrinter->m_directory, filename);
//            //singleTabulater.setCaption(caption);
//            size_t i=0;
//            singleTabulater.addColumn(N_SYMBOL, nLevels,0, i++);
//
//            vector<string> AlgorithmNamesSmallText(ALGORITHM_NAMES);
////            transform(ALGORITHM_NAMES.begin(),ALGORITHM_NAMES.end(),back_inserter(AlgorithmNamesSmallText),
////                [](const auto& text){
////                    return "\\tiny{" + text + "}";
////            });
//
//            for(int alg=BoundedDegreePlaneSpannerAlgorithm::AlgorithmFirst;
//                alg!=BoundedDegreePlaneSpannerAlgorithm::AlgorithmLast; ++alg ) {
//                vector<string> singleColumn;
//                getSingleIVColumn( iv, BoundedDegreePlaneSpannerAlgorithm(alg), results, singleColumn);
//                singleTabulater.addColumn(AlgorithmNamesSmallText[alg], singleColumn, -1,i++);
//            }
//            singleTabulater.tabulate(TablePrinter::CellHighlightStyle::MaxInRow);
//            addToPrinter->addToDocument(singleTabulater, PRECOMPILE_SUBDOCUMENTS);
//        }
//    }

}

#endif //SPANNERS_TABLEPRINTER_H
