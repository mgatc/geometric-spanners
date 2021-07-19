//
// Created by matt on 7/19/21.
//

#ifndef GEOMETRIC_SPANNERS_LATEXPRINTER_H
#define GEOMETRIC_SPANNERS_LATEXPRINTER_H

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>

#include "names.h"

namespace unf_planespanners {

    using namespace std;

    class LatexPrinter {
    public:
        typedef pair<string,string> Option;
        typedef vector<Option> OptionsList;

        LatexPrinter(string filename, string documentType = "article")
            : m_filename(filename), m_documentType(documentType) {}

        // Document-level getters
        string getName() const {
            return m_filename + m_bodySuffix;
        }
        string getFullDocumentText() const {
            return getDocumentHeader()
                   + getBodyText()
                   + getDocumentFooter();
        }
        string getBodyText() const {
            return m_body.header + m_body.content + m_body.footer;
        }
        string getDocumentHeader() const {
            string header = "\\documentclass{"
                    + m_documentType
                    + "}\n\n"
                    + "\\usepackage{tikz,pgfplots}\n"
                    + "\\usetikzlibrary{shapes}\n"
                    + "\\pgfplotsset{compat=1.15}\n\n"
                    + getColorDefinitions()
                    + "\n\n\n\n"
                    + "\\begin{document}\n\n";
            return header;
        }
        string getDocumentFooter() const {
            return "\\end{document}";
        }
        // Figure getters
        string getFigureHeader(string options = "") const {
            string header = "\\begin{figure}";
            if(!options.empty()){
                header += "[" + options + "]";
            }
            header += "\n\\centering\n\n";
            return header;
        }
        string getFigureFooter() const {
            string footer;
            if( !m_caption.empty() ) {
                footer += string("")
                    + "\\caption{"
                    + m_caption
                    + "}\n";
            }
            footer += "\\end{figure}\n\n";
            return footer;
        }
        // Figure captions
        string getCaption() const {
            return "\\caption{"
                + m_caption
                + "}\n";
        }
        void setCaption(const string& caption ) {
            m_caption = caption;
            addLatexComment(caption);
        }
        void setCaption( const Result& result ) {
            string caption = string("\\textsc{")
                             + Names.at(result.algorithm)
                             + "}: "
                             + "$\\Delta = "
                             + to_string(result.degree);

            if( result.t ) {
                caption += ",\\ t = "
                           + to_string(*result.t);
            }

            caption +="$";
            setCaption(caption);
        }
        void clearCaptionFile() const {
            string captionFilename = "captions.txt";

            FILE *fp = fopen(captionFilename.c_str(), "w");
            fprintf(fp, "%s", "");
            fclose(fp);
        }
        // Adding to the document's body
        /** Add the contents of a LatexPrinter object to this document*/
        void addToDocument(const LatexPrinter& printer) {
            printer.saveBody();
            addInput(printer.getName());

            for( auto hex : printer.m_colors ) {
                defineColor(hex);
            }
        }
        /** Add the contents of a LatexPrinter object to this document as a figure. */
        void addToDocumentAsFigure(const LatexPrinter& printer, const bool captionAbove = true) {

            addRawText(getFigureHeader());

            string caption = printer.getCaption();
            if(!caption.empty() && captionAbove)
                addRawText(caption);

            addToDocument(printer);

            if(!caption.empty() && !captionAbove)
                addRawText(caption);

            addRawText(getFigureFooter());
        }
        void addInput(string name) {
            m_body.content += "\\input{" + name + "}\n\n";
        }
        void addRawText(const string& text) {
            m_body.content += text;
        }
        void addLatexComment( const string& comment ) {
            addRawText("% " + comment + "\n");
        }

        void save() const {
            string texFilename = getTexFilename();
            cout<<"Saving file "<<texFilename<<"...";

            FILE *fileOut = fopen(texFilename.c_str(), "w");

            fprintf(fileOut, "%s", getFullDocumentText().c_str());

            fclose(fileOut);
            cout<<"done.\n";
        }
        void saveBody() const {
            string texFilename = getTexFilenameForBody();
            cout<<"Saving file "<<texFilename<<"...";
            FILE *fileOut = fopen(texFilename.c_str(), "w");

            fprintf(fileOut, "%s", getBodyText().c_str());

            fclose(fileOut);
            cout<<"done.\n";
        }
        void compile() const {
            save();

            cout<<"Compiling "<< getTexFilename()<<"...";
            string command = "pdflatex " + getTexFilename() + " > /dev/null";
            ignore = system(command.c_str());
            cout<<"done.\n";
        }
        void display() const {
            compile();

            cout<<"Opening "<<getPdfFilename()<<" for viewing.\n";
            string command = "evince " + getPdfFilename() + " &";
            ignore = system(command.c_str());
        }
        /** Define a hex rgb color for use in the document. */
        void defineColor( const string& hex ) {
            // add color to colormap
            m_colors.insert( hex );
        }

    protected:
        struct Body {
            string header;
            string content;
            string footer;
        };
        Body m_body;
        string m_filename;
        string m_documentType;
        string m_caption;
        unordered_set<string> m_colors;
        string m_bodySuffix = "_body";

        string getTexFilename() const {
            return m_filename + ".tex";
        }
        string getTexFilenameForBody() const {
            return getName() + ".tex";
        }
        string getPdfFilename() const {
            return m_filename + ".pdf";
        }
        static string expandOptions( const OptionsList& options ) {
            string optionsString;
            for( auto& o : options ) {
                optionsString += o.first
                                 + ( !o.second.empty()?("=" + o.second):"") // include second param if given
                                 + ",";
            }
            // remove trailing comma
            return optionsString.substr( 0, optionsString.size()-1 );
        }

        string getColorDefinitions() const {
            string definitions;
            for( auto hex : m_colors ) {
                // parse the hex value
                vector<size_t> color = parseHexRGB( hex );
                // add color to document
                definitions += "\\definecolor{"
                               +  hex + "}{RGB}{ "
                               +  to_string(color.at(0)) + ", "
                               +  to_string(color.at(1)) + ", "
                               +  to_string(color.at(2))
                               + " }\n";
            }
            definitions += "\n";
            return definitions;
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
    };

} // unf_planespanners

#endif //GEOMETRIC_SPANNERS_LATEXPRINTER_H
