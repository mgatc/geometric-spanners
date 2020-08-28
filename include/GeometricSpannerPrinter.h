#ifndef GSNUNF_GEOMETRICSPANNERPRINTER_H
#define GSNUNF_GEOMETRICSPANNERPRINTER_H

#include <cstdlib>
#include <iostream>
#include <list>
#include <string>
#include <tuple> // ignore
#include <utility> // pair

#include "DelaunayGraph.h"

namespace gsnunf {

using namespace std;

class GeometricSpannerPrinter {
  public:

    double radiusOfPoints;

    GeometricSpannerPrinter()
        : _header( "\\documentclass{standalone} \n\\usepackage{tikz} \n \n\n\\begin{document}\n\n\n\n\\begin{tikzpicture}\n\n" ),
          _footer( "\n\n\\end{tikzpicture}\n\n\\end{document}") {

    }

    template< typename RandomAccessIterator >
    void drawEdges( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd, const vector<pair<string,string>>& options = {} ) {
        for( auto e=edgesBegin; e!=edgesEnd; ++e ) {
            drawLine( e->first.x(), e->first.y(), e->second.x(), e->second.y(), options );
        }
    }

    void drawEdges( const DelaunayGraph::Delaunay_triangulation_2 &Triangulation, const vector<pair<string,string>>& options = {} ) {
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

    template< typename T >
    void drawVertices( const T &Triangulation, const vector<pair<string,string>>& options = {} ) {
        double radius = 0.09;
        for( typename T::Finite_vertices_iterator it = Triangulation.finite_vertices_begin(); it != Triangulation.finite_vertices_end(); ++it )
            _document += "\\draw [" + parseOptions( options ) + "] ("
                + to_string(it->point().x()) + ","+ to_string(it->point().y()) +") circle [radius="+to_string(radius)+"];\n";

        _document += "\n";
    }
    void drawVertexPair( const pair<Vertex_handle,Vertex_handle>& vertices, const vector<pair<string,string>>& options = {} ) {
        string optionsStr = parseOptions( options );
        double radius = 3;
        _document += "\\draw [" + optionsStr + "] ("
            + to_string(vertices.first->point().x()) + ","+ to_string(vertices.first->point().y()) +") circle [radius="+to_string(radius)+"];\n";
        _document += "\\draw [" + optionsStr + "] ("
            + to_string(vertices.second->point().x()) + ","+ to_string(vertices.second->point().y()) +") circle [radius="+to_string(radius)+"];\n";
    }

    void drawEdges( const DelaunayGraph& DG, const vector<pair<string,string>>& options = {} ) {
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

    void drawLine( double x1, double y1, double x2, double y2, const vector<pair<string,string>>& options = {} ) {
        _document += "\\draw [" + parseOptions( options ) + "] ("
            + to_string(x1) + "," + to_string(y1) + ") -- ("
            + to_string(x2) + "," + to_string(y2) + ");\n";
    }

    string parseOptions( const vector<pair<string,string>>& options ) {
        string optionsString("");
        for( auto& o : options ) {
            optionsString += o.second + ",";//o.first + "=" + o.second + ",";
        }
        return optionsString;
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
        ignore = system(command.c_str());
        cout << "PDF generation terminated...\n";

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
