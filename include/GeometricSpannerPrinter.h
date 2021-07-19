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
#include "names.h"
#include "utilities.h"

namespace unf_planespanners {

using namespace std;

class GraphPrinter {
  public:
    typedef pair<string,optional<string>> Option;
    typedef vector<Option> OptionsList;
    unordered_set<string> _colors;

    string activeEdgeColor =     "000000";
    string inactiveEdgeColor =   "bbbbbb";
    string worstPathEdgeColor =  "ff0000";
    string activeVertexColor =   "0000ff";
    string inactiveVertexColor = "f7b267";
    string backgroundColor =     "FEFEF6";
    string textColor =           "111116";
    double vertexRadius = 4.0;
    //double vertexBorderWidth = 0.63;
    double activeEdgeWidth =0.6;
    double inactiveEdgeWidth = 1.0;

    explicit GraphPrinter( double scale = 1.0 )
      : _scaleFactor(scale), _resizeFactor(1) {

    }
    // begin and end are iterators over the point set that will be printed
    // this is for proper scaling of the points to work with the latex document
    // the points will not be printed unless you call drawVertices() on them
    template< class InputIterator>
    void autoscale( InputIterator pointsBegin, InputIterator pointsEnd, double documentSize = 10 )
    {
        double minX = pointsBegin->x(),
               maxX = minX,
               minY = pointsBegin->y(),
               maxY = minY;

        for( auto p=pointsBegin; p!=pointsEnd; ++p )
        {
            minX = CGAL::min( p->x(), minX );
            maxX = CGAL::max( p->x(), maxX );
            minY = CGAL::min( p->y(), minY );
            maxY = CGAL::max( p->y(), maxY );
        }

        double deltaX = abs( minX - maxX ),
               deltaY = abs( minY - maxY ),
               delta  = CGAL::max( deltaX, deltaY );

        _scaleFactor = documentSize / delta;
    }

    void clear() { _document = ""; }

    template< typename RandomAccessIterator, typename PointContainer >
    void drawEdges( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd, const PointContainer &P, const OptionsList& options = {} ) {
        for( auto e=edgesBegin; e!=edgesEnd; ++e ) {
            drawLine( P[e->first].x(),
                      P[e->first].y(),
                      P[e->second].x(),
                      P[e->second].y(), options );
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

    template< typename Triangulation >
    void drawEdgesOfSDG( const Triangulation& T, const OptionsList& options = {} ) {
        for( auto eit = T.finite_edges_begin(); eit != T.finite_edges_end(); ++eit ) {
            auto e = *eit;
            double x1 = e.first->vertex( (e.second+1)%3 )->site().point().x();
            double y1 = e.first->vertex( (e.second+1)%3 )->site().point().y();
            double x2 = e.first->vertex( (e.second+2)%3 )->site().point().x();
            double y2 = e.first->vertex( (e.second+2)%3 )->site().point().y();
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

    template< typename InputIterator >
    void drawVertices( const InputIterator &pointsStart, const InputIterator &pointsEnd, const OptionsList& options = {}, const OptionsList& borderOptions = {} ) {
        size_t id = 0;
        for( auto it=pointsStart; it!=pointsEnd; ++it )
            drawVertexWithLabel( it->x(), it->y(), "", options, borderOptions );
        _document += "\n";
    }

    template< typename T >
    void drawVerticesWithInfo( const T &Triangulation, const OptionsList& options = {}, const OptionsList& borderOptions = {} ) {
        for( typename T::Finite_vertices_iterator it = Triangulation.finite_vertices_begin(); it != Triangulation.finite_vertices_end(); ++it )
            drawVertexWithLabel( it->point().x(), it->point().y(), to_string(it->info()), options, borderOptions );
        _document += "\n";
    }

    template< typename InputIterator >
    void drawVerticesWithInfo( const InputIterator &pointsStart, const InputIterator &pointsEnd, const OptionsList& options = {}, const OptionsList& borderOptions = {} ) {
        size_t id = 0;
        for( auto it=pointsStart; it!=pointsEnd; ++it )
            drawVertexWithLabel( it->x(), it->y(), to_string(id++), options, borderOptions );
        _document += "\n";
    }

    template< typename T >
    void drawVerticesWithInfoSDG( const T &Triangulation, const OptionsList& options = {}, const OptionsList& borderOptions = {} ) {
        for( typename T::Finite_vertices_iterator it = Triangulation.finite_vertices_begin(); it != Triangulation.finite_vertices_end(); ++it )
            drawVertexWithLabel( it->site().point().x(), it->site().point().y(), to_string(it->storage_site().info()), options, borderOptions );
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

    void drawVertexWithLabel( double x, double y, const string &label, const OptionsList& options = {}, const OptionsList& borderOptions = {} ) {

        if( label.empty() ){
            stringstream stream;
            stream << fixed << setprecision(1);
            stream << x << "," << y;
            //label = stream.str();
        }
        _document += "\\node (vertex" + label + ") [fill,"
                + expandOptions( options )
            + "] at ("
                + to_string( x*_scaleFactor ) + ","
                + to_string( y*_scaleFactor )
            + ") {"
                + label
            + "};\n";
        if( !borderOptions.empty() ) {
            _document += "\\node () [draw,"
                + expandOptions( borderOptions )
            + "] at ("
                + to_string( x*_scaleFactor ) + ","
                + to_string( y*_scaleFactor )
            + ") {};\n\n";
        }

    }

    void drawEdges( const DelaunayGraph& DG, const OptionsList& options = {} ) {
        for( const auto &el : DG.m_E ) {
            for( const auto& v : el.second ) {
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

    template< typename Triangulation >
    void drawEdgesOfHalfTheta( const Triangulation& T, const OptionsList& options = {} ) {
                using CGAL::to_double;
        //typedef typename Triangulation::VertexHandle VertexHandle;
        for( auto eit = T.edges_begin(); eit != T.edges_end(); ++eit ) {
            auto e = *eit;

            typename Triangulation::Point_2
                vp = T._G[source(e,T._G)],
                vq = T._G[target(e,T._G)];

            double x1 = to_double(vp.x());
            double y1 = to_double(vp.y());
            double x2 = to_double(vq.x());
            double y2 = to_double(vq.y());
            drawLine( x1,y1,x2,y2,options );
        }
        _document += "\n";
    }

    void drawLine( double x1, double y1, double x2, double y2, const OptionsList& options = {} ) {
        _document += "\\draw [" + expandOptions( options ) + "] ("
            + to_string(x1*_scaleFactor) + "," + to_string(y1*_scaleFactor) + ") -- ("
            + to_string(x2*_scaleFactor) + "," + to_string(y2*_scaleFactor) + ");\n";
    }
    string defineColor( const string& hex ) {
        // parse the hex value
        vector<size_t> color = parseHexRGB( hex );
        // add color to document
        string definition = string("")
                    + "\\definecolor{"
                    +  hex + "}{RGB}{ "
                        +  to_string(color.at(0)) + ", "
                        +  to_string(color.at(1)) + ", "
                        +  to_string(color.at(2))
                    + " }\n";
        // add color to colormap
        _colors.insert( hex );
        return definition;
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
    void addLatexComment( const string& comment ) {
        _document += "% " + comment + "\n";
    }
    void print( string fName ) {
        fName = _outputFilePrefix + fName;
        string texFilename = fName + ".tex",
               pdfFilename = fName + ".pdf",
               captionFilename = _outputFilePrefix + "captions.txt"
                                         "";
        boost::erase_all(fName, ".");

        FILE *fp = fopen( texFilename.c_str(), "w" );

        fprintf( fp, "%s", getDocumentHeader().c_str() );
        fprintf( fp, "%s", _document.c_str() );
        fprintf( fp, "%s", getDocumentFooter().c_str() );

        fclose(fp);

        fp = fopen( captionFilename.c_str(), "a" );

        auto captionWithNewLine = _caption + "\n";
        fprintf( fp, "%s", captionWithNewLine.c_str() );

        fclose(fp);
        //cout<<fName<<endl;
        //cout<< _header<<_document<<_footer;

//        cout << "\nOutput PDF generation started...\n";
        string command = "pdflatex " + texFilename + " > /dev/null";
        ignore = system(command.c_str());
//        cout << "PDF generation terminated...\n";

        command = "evince " + pdfFilename + " &";
        ignore = system(command.c_str());
    }
    void clearCaptionFile() {
        string captionFilename = _outputFilePrefix + "captions.txt";

        FILE *fp = fopen( captionFilename.c_str(), "w" );
        fprintf( fp, "%s", "" );
        fclose(fp);
    }
    void setCaption(const string& caption ) {
        _caption = caption;
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
    string getDocumentHeader() {
        string documentHeader = string("")
                                 + "\\documentclass{standalone}\n"
                                 + "\\usepackage{tikz}\n"
                                 //+ "\\usetikzlibrary{backgrounds}\n\n"
                                 + "\\usetikzlibrary{shapes}\n\n";

        // define colors in the document
        string colors = defineColor(activeEdgeColor)
            + defineColor(inactiveEdgeColor)
            + defineColor(worstPathEdgeColor)
            + defineColor(activeVertexColor)
            + defineColor(inactiveVertexColor)
            + defineColor(backgroundColor)
            + defineColor(textColor);

        documentHeader += colors
                + "\\begin{document}\n\n";

        return documentHeader;
    }
    string getDocumentFooter() {
        string footer = "\\end{document}";
        return footer;
    }
    string getFigureFooter() {
        string footer = string("")
                        + "\n\n\\end{tikzpicture}\n";
                        //+ "}\n\n";
        //footer += string("")
                 // + "\\caption{"
                 // + _caption
                 // +"}\n";
                 // + "\\end{figure}\n\n";
        return footer;
    }
    void beginFigure() {
        _document += _figureHeader;
    }
    void endFigure() {
        _document += getFigureFooter();
    }


private:
    //Result _result;
    double _scaleFactor;
    double _resizeFactor;
    string _caption;
    string _outputFilePrefix = "Demo-";
    string _document;
    string _figureHeader = string("")
                           //+ "\\begin{figure}\n\n"
                           //+ "\\resizebox{\\textwidth}{!}{\n"
                           + "\\begin{tikzpicture}[ "
                           //+ "background rectangle/.style={fill="
                           //    + backgroundColor
                           //+ "}, show background rectangle, "
                           + "vertex/.style = {fill, minimum size=#1, inner sep=0pt, outer sep=0pt}, "
                           + "vertex/.default = 6pt, "
                           + "border/.style = {draw, minimum size=#1, inner sep=0pt, outer sep=0pt}, "
                           + "border/.default = 6pt ]\n\n"
                           + "\\draw[step=1.0,black,thin,dotted] (-5.5,-5.5) grid (5.5,5.5);\n\n"
                           + "\\tiny\n\n";

}; // class GraphPrinter

} // namespace unf_planespanners

#endif // GSNUNF_GRAPHPRINTER_H
