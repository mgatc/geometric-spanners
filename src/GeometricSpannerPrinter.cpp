#include "GeometricSpannerPrinter.h"

namespace gsnunf{

GeometricSpannerPrinter::GeometricSpannerPrinter(double radiusOfPoints) {
    this->radiusOfPoints = radiusOfPoints;
}

void GeometricSpannerPrinter::drawTriangulation(const DelaunayTriangulation &T, std::string fName) {
    FILE *fp = fopen(  fName.c_str() ,"w");
    fprintf(fp,"\\documentclass{standalone} \n\\usepackage{tikz} \n \n\n\\begin{document}\n");
    fprintf(fp,"\n\n\n\\begin{tikzpicture}\n\n");

    for(Finite_edges_iterator eit = T.finite_edges_begin(); eit != T.finite_edges_end(); ++eit) {
       DelaunayTriangulation::Edge e = *eit;
       double x1 = e.first->vertex( (e.second+1)%3 )->point().x();
       double y1 = e.first->vertex( (e.second+1)%3 )->point().y();
       double x2 = e.first->vertex( (e.second+2)%3 )->point().x();
       double y2 = e.first->vertex( (e.second+2)%3 )->point().y();
       fprintf(fp, "\\draw [thick] (%f,%f) -- (%f,%f);\n", x1,y1, x2,y2);
    }

    for(Finite_vertices_iterator it = T.finite_vertices_begin(); it != T.finite_vertices_end(); ++it)
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
void GeometricSpannerPrinter::drawGraph(Graph &graph, std::string fName) {
    FILE *fp = fopen(  fName.c_str() ,"w");
    fprintf(fp,"\\documentclass{standalone} \n\\usepackage{tikz} \n \n\n\\begin{document}\n");
    fprintf(fp,"\n\n\n\\begin{tikzpicture}\n\n");

    // Draw vertices
    for(Finite_vertices_iterator it = graph._DT.finite_vertices_begin(); it != graph._DT.finite_vertices_end(); ++it)
        fprintf(fp,"\\draw [fill=red,stroke=red] (%f,%f) circle [radius=%f];\n",it->point().x(),it->point().y(),radiusOfPoints);

    // draw spanning graph
    for( auto el : graph._E ) {
        for( auto v : el.second ) {
            std::cout << el.first->point() << " " << v->point() << "\n";
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
}
