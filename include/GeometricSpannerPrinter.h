#ifndef GSNUNF_GEOMETRICSPANNERPRINTER_H
#define GSNUNF_GEOMETRICSPANNERPRINTER_H

#include "CGALComponents.h"
#include "DelaunayGraph.h"

namespace gsnunf {

class GeometricSpannerPrinter {
  public:

    double radiusOfPoints;

    GeometricSpannerPrinter( double radiusOfPoints = 0.09 ) : radiusOfPoints(radiusOfPoints) {}

    template< typename T >
    void draw( const T &Triangulation, std::string fName ) {
        FILE *fp = fopen(  fName.c_str() ,"w");
        fprintf(fp,"\\documentclass{standalone} \n\\usepackage{tikz} \n \n\n\\begin{document}\n");
        fprintf(fp,"\n\n\n\\begin{tikzpicture}\n\n");

        for( typename T::Finite_edges_iterator eit = Triangulation.finite_edges_begin(); eit != Triangulation.finite_edges_end(); ++eit ) {
           typename T::Edge e = *eit;
           double x1 = e.first->vertex( (e.second+1)%3 )->point().x();
           double y1 = e.first->vertex( (e.second+1)%3 )->point().y();
           double x2 = e.first->vertex( (e.second+2)%3 )->point().x();
           double y2 = e.first->vertex( (e.second+2)%3 )->point().y();
           fprintf(fp, "\\draw [thick] (%f,%f) -- (%f,%f);\n", x1,y1, x2,y2);
        }

        for( typename T::Finite_vertices_iterator it = Triangulation.finite_vertices_begin(); it != Triangulation.finite_vertices_end(); ++it )
             fprintf(fp,"\\draw [fill=red,stroke=red] (%f,%f) circle [radius=%f];\n",it->point().x(),it->point().y(),radiusOfPoints);


        fprintf(fp,"\n\n\\end{tikzpicture}");
        fprintf(fp,"\n\n\\end{document}");
        fclose(fp);

        std::cout << "\nOutput PDF generation started...\n";
        std::string command = "pdflatex " + fName + " > /dev/null";
        system(command.c_str());
        std::cout << "PDF generation terminated...\n";

        command = "evince " + fName + ".pdf &";
        system(command.c_str());

    }

    template< typename T >
    void draw( const DelaunayGraph<T>& DG, std::string fName ) {
        FILE *fp = fopen(  fName.c_str() ,"w");
        fprintf(fp,"\\documentclass{standalone} \n\\usepackage{tikz} \n \n\n\\begin{document}\n");
        fprintf(fp,"\n\n\n\\begin{tikzpicture}\n\n");

        // Draw vertices
        for( typename T::Finite_vertices_iterator it = DG._DT.finite_vertices_begin(); it != DG._DT.finite_vertices_end(); ++it)
            fprintf(fp,"\\draw [fill=red,stroke=red] (%f,%f) circle [radius=%f];\n",it->point().x(),it->point().y(),radiusOfPoints);

        // draw spanning graph
        for( auto el : DG._E ) {
            for( auto v : el.second ) {
                //std::cout << el.first->point() << " " << v->point() << "\n";
                double x1 = el.first->point().x();
                double y1 = el.first->point().y();
                double x2 = v->point().x();
                double y2 = v->point().y();
                fprintf(fp, "\\draw [thick,stroke=red] (%f,%f) -- (%f,%f);\n", x1,y1, x2,y2);
            }
        }

        fprintf(fp,"\n\n\\end{tikzpicture}");
        fprintf(fp,"\n\n\\end{document}");
        fclose(fp);

        std::cout << "\nOutput PDF generation started...\n";
        std::string command = "pdflatex " + fName + " > /dev/null";
        system(command.c_str());
        std::cout << "PDF generation terminated...\n";

        command = "evince " + fName + ".pdf &";
        system(command.c_str());

    }

}; // class TriangulationPrinter

} // namespace gsnunf

#endif // GSNUNF_GEOMETRICSPANNERPRINTER_H
