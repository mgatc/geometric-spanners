#ifndef GSNUNF_GRAPHPRINTER_H
#define GSNUNF_GRAPHPRINTER_H


#include <cctype>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <list>
#include <optional>
#include <sstream>
#include <string>
#include <tuple> // ignore
#include <utility> // pair

#include <boost/algorithm/string.hpp>

#include "DelaunayGraph.h"

namespace gsnunf {

using namespace std;

class GraphPrinter {
  public:
    typedef pair<string,optional<string>> Option;
    typedef vector<Option> OptionsList;

    double _scaleFactor;
    unordered_set<string> _colors;

    string activeEdgeColor =     "e98570";
    string inactiveEdgeColor =   "95acc3";
    string activeVertexColor =   "174681";
    string inactiveVertexColor = "95acc3";
    string backgroundColor =     "FEFEF6";
    string textColor =           "111116";
    double vertexRadius = 17;
    //double vertexBorderWidth = 0.63;
    double activeEdgeWidth = 2.0;
    double inactiveEdgeWidth = 1.36;

    GraphPrinter( double scale = 1.0 )
      : _scaleFactor(scale) {
        // define colors in the document
        defineColor(activeEdgeColor);
        defineColor(inactiveEdgeColor);
        defineColor(activeVertexColor);
        defineColor(inactiveVertexColor);
        defineColor(backgroundColor);
        defineColor(textColor);
    }

    void clear() { _document = ""; }

    template< typename RandomAccessIterator >
    void drawEdges( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd, const OptionsList& options = {} ) {
        for( auto e=edgesBegin; e!=edgesEnd; ++e ) {
            drawLine( e->first.x(), e->first.y(), e->second.x(), e->second.y(), options );
        }
    }

    template< typename Triangulation >
    void drawEdges( const Triangulation& T, const OptionsList& options = {} ) {
        for( auto eit = T.finite_edges_begin(); eit != T.finite_edges_end(); ++eit ) {
            auto e = *eit;
            double x1 = e.first->vertex( (e.second+1)%3 )->point().x();
            double y1 = e.first->vertex( (e.second+1)%3 )->point().y();
            double x2 = e.first->vertex( (e.second+2)%3 )->point().x();
            double y2 = e.first->vertex( (e.second+2)%3 )->point().y();
            drawLine( x1,y1,x2,y2,options );
        }
        _document += "\n";
    }

    template< typename T >
    void drawVertices( const T &Triangulation, const OptionsList& options = {}, const OptionsList& borderOptions = {} ) {
        for( typename T::Finite_vertices_iterator it = Triangulation.finite_vertices_begin(); it != Triangulation.finite_vertices_end(); ++it )
            drawVertexWithLabel( it->point().x(), it->point().y(), "", options, borderOptions );
        _document += "\n";
    }

    template< typename T >
    void drawVerticesWithInfo( const T &Triangulation, const OptionsList& options = {}, const OptionsList& borderOptions = {} ) {
        for( typename T::Finite_vertices_iterator it = Triangulation.finite_vertices_begin(); it != Triangulation.finite_vertices_end(); ++it )
            drawVertexWithLabel( it->point().x(), it->point().y(), to_string(it->info()), options, borderOptions );
        _document += "\n";
    }
    template< typename T >
    void drawVertexPair( const pair<typename T::Vertex_handle,typename T::Vertex_handle>& vertices, const OptionsList& options = {} ) {
        drawVertex( vertices.first->point().x(), vertices.first->point().y(), options );
        drawVertex( vertices.second->point().x(), vertices.second->point().y(), options );
    }

    void drawVertex( double x, double y, const OptionsList& options = {} ) {
        drawVertexWithLabel( x, y, "", options );
    }

    void drawVertexWithLabel( double x, double y, string label, const OptionsList& options = {}, const OptionsList& borderOptions = {} ) {

        if( label == "" ){
            stringstream stream;
            stream << fixed << setprecision(1);
            stream << x << "," << y;
            label = stream.str();
        }
        _document += "\\node (vertex" + label + ") [fill,"
                + expandOptions( options )
            + "] at ("
                + to_string( x*_scaleFactor ) + ","
                + to_string( y*_scaleFactor )
            + ") {"
                + label
            + "};\n"
            + "\\node () [draw,"
                + expandOptions( borderOptions )
            + "] at ("
                + to_string( x*_scaleFactor ) + ","
                + to_string( y*_scaleFactor )
            + ") {};\n\n";

    }

    void drawEdges( const DelaunayGraph& DG, const OptionsList& options = {} ) {
        for( auto el : DG._E ) {
            for( auto v : el.second ) {
                //std::cout << el.first->point() << " " << v->point() << "\n";
                double x1 = el.first->point().x();
                double y1 = el.first->point().y();
                double x2 = v->point().x();
                double y2 = v->point().y();
                drawLine( x1,y1,x2,y2,options );
            }
        }
        _document += "\n";
    }

    void drawLine( double x1, double y1, double x2, double y2, const OptionsList& options = {} ) {
        _document += "\\draw [" + expandOptions( options ) + "] ("
            + to_string(x1*_scaleFactor) + "," + to_string(y1*_scaleFactor) + ") -- ("
            + to_string(x2*_scaleFactor) + "," + to_string(y2*_scaleFactor) + ");\n";
    }
    void defineColor( string hex ) {
        // parse the hex value
        vector<size_t> color = parseHexRGB( hex );
        // add color to document
        _document += "\\definecolor{"
                  +  hex + "}{RGB}{ "
                      +  to_string(color.at(0)) + ", "
                      +  to_string(color.at(1)) + ", "
                      +  to_string(color.at(2))
                  + " }\n";
        // add color to colormap
        _colors.insert( hex );
    }

    string expandOptions( const OptionsList& options ) {
        string optionsString("");
        for( auto& o : options ) {
            optionsString += o.first
                + ( o.second?("=" + *o.second):"") // include second param if given
                + ",";
        }
        // remove trailing comma
        return optionsString.substr( 0, optionsString.size()-1 );
    }
    vector<size_t> parseHexRGB( string hex_str ) {
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
        boost::erase_all(fName, ".");
        FILE *fp = fopen( fName.c_str(), "w" );

        fprintf( fp, "%s", _header.c_str() );
        fprintf( fp, "%s", _document.c_str() );
        fprintf( fp, "%s", _footer.c_str() );

        fclose(fp);
        //cout<< _header<<_document<<_footer;

//        cout << "\nOutput PDF generation started...\n";
        string command = "pdflatex " + fName + " > /dev/null";
        ignore = system(command.c_str());
//        cout << "PDF generation terminated...\n";

        command = "evince " + fName + ".pdf &";
        ignore = system(command.c_str());
    }

  private:
    string _document = "";
    string _header = string("")
        + "\\documentclass[margin=3mm]{standalone} \n"
        + "\\usepackage{tikz}  \n "
        + "\\usetikzlibrary{backgrounds}\n\n"
        + "\\usetikzlibrary{fit}\n\n"
        + "\\begin{document}\n\n"
        + "\\begin{tikzpicture}[ background rectangle/.style={fill="
            + backgroundColor
        + "}, show background rectangle, "
        + "vertex/.style = {circle, fill, minimum size=#1, inner sep=0pt, outer sep=0pt}, "
        + "vertex/.default = 6pt, "
        + "border/.style = {circle, draw, minimum size=#1, inner sep=0pt, outer sep=0pt}, "
        + "border/.default = 6pt ]\n\n";
    string _footer = string("")
        + "\n\n\\end{tikzpicture}\n\n"
        + "\\end{document}";

}; // class GraphPrinter

} // namespace gsnunf

#endif // GSNUNF_GRAPHPRINTER_H
