#ifndef GSNUNF_GEOMETRICSPANNERPRINTER_H
#define GSNUNF_GEOMETRICSPANNERPRINTER_H

#include <cstdlib>
#include <iostream>
#include <list>
#include <optional>
#include <string>
#include <tuple> // ignore
#include <utility> // pair

#include "DelaunayGraph.h"

namespace gsnunf {

using namespace std;

class GeometricSpannerPrinter {
  public:
    using tikzOption = pair< string,optional<string> >;
    double radiusOfPoints;

    GeometricSpannerPrinter()
        : _header( "\\documentclass{standalone} \n\\usepackage{tikz} \n \n\n\\begin{document}\n\n\n\n\\begin{tikzpicture}[scale=1.0, x=0.10cm, y=0.10cm]\n\n" ),
          _footer( "\n\n\\end{tikzpicture}\n\n\\end{document}") {

    }

    void drawEdges( const DelaunayGraph::Delaunay_triangulation_2 &Triangulation, const vector<tikzOption>& options = {} ) {

        for( auto eit = Triangulation.finite_edges_begin(); eit != Triangulation.finite_edges_end(); ++eit ) {
            auto e = *eit;
            double x1 = e.first->vertex( (e.second+1)%3 )->point().x();
            double y1 = e.first->vertex( (e.second+1)%3 )->point().y();
            double x2 = e.first->vertex( (e.second+2)%3 )->point().x();
            double y2 = e.first->vertex( (e.second+2)%3 )->point().y();
            drawLine( x1,y1,x2,y2,options );
        }
        _document += "\n";
    }

    void drawEdges( const DelaunayGraph& DG, const vector<tikzOption>& options = {} ) {

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

    template< typename T >
    void drawVertices( const T &Triangulation, double size, const vector<tikzOption>& options = {} ) {
        for( typename T::Finite_vertices_iterator it = Triangulation.finite_vertices_begin(); it != Triangulation.finite_vertices_end(); ++it )
            _document += "\\filldraw [" + parseOptions( options ) + "] ("
                + to_string(it->point().x()) + ","+ to_string(it->point().y()) +") circle ("+to_string(size)+");\n";

        _document += "\n";
    }

    void drawLine( double x1, double y1, double x2, double y2, const vector<tikzOption>& options = {} ) {
        _document += "\\draw [" + parseOptions( options ) + "] ("
            + to_string(x1) + "," + to_string(y1) + ") -- ("
            + to_string(x2) + "," + to_string(y2) + ");\n";
    }

    string parseOptions( const vector<tikzOption>& options ) {
        if( options.empty() ) return "";

        string optionsString("");
        for( auto& o : options ) {
            //optionsString += o.second + ",";
            optionsString += o.first;
            optionsString += o.second ? ("=" + *o.second) : "";
            optionsString += ",";
        }
        // remove trailing comma
        return optionsString.substr( 0, optionsString.size()-1 );
    }

    void print( string fName ) {
        FILE *fp = fopen( fName.c_str(), "w" );

        fprintf( fp, "%s", _header.c_str() );
        fprintf( fp, "%s", _document.c_str() );
        fprintf( fp, "%s", _footer.c_str() );

        fclose(fp);
        //cout<< _header<<_document<<_footer;

        cout << "\nOutput PDF generation started...\n";
        string command = "pdflatex " + fName + " > /dev/null";
        int pdflatexResult = system(command.c_str());
        cout << "PDF generation terminated with status "
             << pdflatexResult << "\n";

        command = "evince " + fName + ".pdf &";
        ignore = system(command.c_str());
    }

  private:
    string _document;
    string _header;
    string _footer;

}; // class TriangulationPrinter

} // namespace gsnunf

#endif // GSNUNF_GEOMETRICSPANNERPRINTER_H
