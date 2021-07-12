#ifndef GEOMETRIC_SPANNERS_LINEGRAPHPRINTER_H
#define GEOMETRIC_SPANNERS_LINEGRAPHPRINTER_H


#include <map>
#include <iostream>
#include <optional>
#include <vector>

#include "utilities.h"

namespace unf_planespanners {



    using namespace std;

    class LineGraphPrinter {
    public:
        typedef map< Algorithm, map< size_t, vector<Result>>> ResultMap;
        typedef pair<string,optional<string>> Option;
        typedef vector<Option> OptionsList;

        unordered_set<string> _colors;

        // Palette generated using https://coolors.co/
        vector<string> Colors = {
                "f94144","f3722c","f8961e",
                "f9844a","f9c74f","90be6d",
                "43aa8b","4d908e","577590",
                "277da1","360568"
        };
        vector<string> Names = {
                "BGS2005",
                "LW2004",
                "BSX2009",
                "KPX2010",
                "KX2012",
                "BCC2012-7",
                "BCC2012-6",
                "BHS2017",
                "BGHP2010",
                "KPT2017",
                "BKPX2015"
        };


        LineGraphPrinter() {
        }

        void clear() { _document = ""; }


        string defineColor( const string& hex ) {
            // parse the hex value
            vector<size_t> color = parseHexRGB( hex );
            // add color to document
            string def("\\definecolor{"
                         +  hex + "}{RGB}{ "
                         +  to_string(color.at(0)) + ", "
                         +  to_string(color.at(1)) + ", "
                         +  to_string(color.at(2))
                         + " }\n");
            // add color to colormap
            _colors.insert( hex );
            return def;
        }

        static string expandOptions( const OptionsList& options ) {
            string optionsString;
            for( auto& o : options ) {
                optionsString += o.first
                                 + ( o.second?("=" + *o.second):"") // include second param if given
                                 + ",";
            }
            // remove trailing comma
            return optionsString.substr( 0, optionsString.size()-1 );
        }
        static vector<size_t> parseHexRGB( const string& hex_str ) {
            // the hex string should contain 6 digits
            // three 2-digit hex numbers
            vector<size_t> rgb(3, 0);
            // parse each 2-digit number and convert to base 10
            for( size_t i=0; i<3&&i<hex_str.size()/2; ++i ) {
                rgb[i] = stoi( hex_str.substr(2*i, 2), 0, 16 );
            }
            return rgb;
        }
        void print( string fName ) {
            fName = _outputFilePrefix + fName;
            boost::erase_all(fName, ".");
            FILE *fp = fopen( fName.c_str(), "w" );
            std::cout<<"Opening file...\n\n";
            if(!fp) {
                std::cout<<"Error opening file!"<<std::endl;
                return;
            }
            fprintf( fp, "%s", getHeader().c_str() );
            fprintf( fp, "%s", getBody().c_str() );
            fprintf( fp, "%s", _footer.c_str() );

            fclose(fp);
            //cout<<fName<<endl;
            //cout<< _header<<_document<<_footer;

//        cout << "\nOutput PDF generation started...\n";
            string command = "pdflatex " + fName + " > /dev/null";
            ignore = system(command.c_str());
//        cout << "PDF generation terminated...\n";

            command = "evince " + fName + ".pdf &";
            ignore = system(command.c_str());
        }

        void registerResult(const Result& result ) {
            auto el1 = _results.begin();
            bool inserted = false;
            tie(el1, inserted) = _results.emplace(result.algorithm,map<size_t, vector<Result>>());

            auto el2 = el1->second.begin();
            tie(el2, inserted) = el1->second.emplace(result.n, vector<Result>());

            // Skip first run
            if(!inserted) {
                el2->second.push_back(result);
            }
        }

        string getHeader() {

            string header = string("")
                         + "\\documentclass[tikz]{standalone}\n"
                         + "\\usepackage{pgfplots}\n"
                         + "\\pgfplotsset{compat=1.15}\n\n"
                         + "\\begin{document}\n\n";

            // define colors in the document
            for( const auto& color : Colors ) {
                header += defineColor(color);
            }
            return header;
        }
        string getBody() {
            const double scale = 1;
            string plotHeader = string("")
                    + "\\begin{tikzpicture}\n\n"
                    + "\\begin{axis}[legend pos=north west,xlabel=$n$, "
                    + "xtick={";

            string xTicks("");
            for( const auto& level : _results.begin()->second ){
                xTicks += to_string(level.first);
                xTicks += ",";
            }
            xTicks = xTicks.substr( 0, xTicks.size()-1 );

            plotHeader += xTicks
                      + "\n}, ylabel near ticks,ylabel={IV}]\n\n";

            vector<string> body(3, plotHeader);
            size_t colorIndex = 0;
            for( auto alg : _results ) {
                for( auto& graph : body) {
                    graph += "\\addplot[thick,color="
                            + Colors[colorIndex]
                            + ",mark=*] coordinates {\n";
                }
                for( auto level : alg.second ) {
                    double averageRuntime = std::accumulate(level.second.begin(), level.second.end(), 0.0,[&]( const auto& a, const auto& b ) {
                        return a + double(b.runtime) / scale;
                    } ) / level.second.size();
                    body[0] += "("
                              + to_string(level.first)
                              + ","
                              + to_string(averageRuntime)
                              + ")\n";

                    double averageDegree = std::accumulate(level.second.begin(), level.second.end(), 0,[]( const auto& a, const auto& b ) {
                        return a + b.degree;
                    } ) / double(level.second.size());
                    body[1] += "("
                               + to_string(level.first)
                               + ","
                               + to_string(averageDegree)
                               + ")\n";

                    if(level.second.front().t) {
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
                    graph += string("}node [pos=1.15, above left] {};\n")
                            + "\\addlegendentry{\\textsc{\\tiny "
                            + Names.at(alg.first)
                            + "}};"
                            + "\n\n";
                }
                ++colorIndex;
            }
            for( auto& graph : body ) {
                graph += string("\n\n\\end{axis}\n\n")
                    + "\\end{tikzpicture}\n\n";
            }
            return std::accumulate( body.begin(), body.end(), string("") );
        }

    private:
        ResultMap _results;
        string _outputFilePrefix = "result-";
        string _document;
        string _footer = string("")
                         + "\\end{document}\n";

    }; // class LineGraphPrinter


    LineGraphPrinter lineGraph;

} // namespace unf_planespanners


#endif //GEOMETRIC_SPANNERS_LINEGRAPHPRINTER_H
