#ifndef GSNUNF_CGALCOMPONENTS_H
#define GSNUNF_CGALCOMPONENTS_H



/**
 *  Global CGAL includes
 */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // CGAL
#include <CGAL/Triangulation_vertex_base_with_info_2.h>         // Triangulations
#include <CGAL/Delaunay_triangulation_2.h>                      // Triangulations
#include <CGAL/point_generators_2.h>                            // Random point generation, testing
#include <CGAL/circulator.h>                                    // Vertex neighbor searches
#include <CGAL/algorithm.h>

/**
 *  Global stl includes
 */
#include <string>
#include <unordered_set>
#include <unordered_map>

/**
 *  Other global includes
 */
#include "Vertex_info.h"

using namespace std;

/**
 *  Global types
 */
typedef CGAL::Exact_predicates_inexact_constructions_kernel                        K;
    typedef K::Point_2                                                         Point;

typedef CGAL::Triangulation_vertex_base_with_info_2< gsnunf::Vertex_info, K >     Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                                 Tds;

typedef CGAL::Delaunay_triangulation_2< K, Tds >               DelaunayTriangulation;
    typedef DelaunayTriangulation::Finite_vertices_iterator Finite_vertices_iterator;
    typedef DelaunayTriangulation::Finite_edges_iterator       Finite_edges_iterator;
    typedef DelaunayTriangulation::Face_handle                           Face_handle;
    typedef DelaunayTriangulation::Vertex_handle                       Vertex_handle;
    typedef DelaunayTriangulation::Vertex_circulator               Vertex_circulator;

typedef CGAL::Creator_uniform_2<double,Point>                                Creator;
typedef CGAL::Container_from_circulator<Vertex_circulator>          Vertex_container;



// TODO: After the TriangulationPrinter class is separated from this file, remove the following two typedefs:
typedef std::unordered_set< Vertex_handle >                        Incident_vertices;
typedef std::unordered_map< Vertex_handle, Incident_vertices >        Adjacency_list;



// For sorting containers of Points, as with std::sort()
// Not used by anything right now, but useful enough to keep around
struct LTR { // comparator for Point type
    bool operator () ( const Point &lhs, const Point &rhs ) {
        return lhs.x() < rhs.x()
            || ( lhs.x() == rhs.x() && lhs.y() < rhs.y() );
    }
};

class TriangulationPrinter {
  public:

    double radiusOfPoints;
    DelaunayTriangulation T;

    TriangulationPrinter(DelaunayTriangulation T, double radiusOfPoints = 0.09) {
        this-> T = T;
        this->radiusOfPoints = radiusOfPoints;
    }

    void draw() {
        std::string fName = "temp";
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
    void drawSpanningGraph( Adjacency_list &graph ) {
        std::string fName = "temp";
        FILE *fp = fopen(  fName.c_str() ,"w");
        fprintf(fp,"\\documentclass{standalone} \n\\usepackage{tikz} \n \n\n\\begin{document}\n");
        fprintf(fp,"\n\n\n\\begin{tikzpicture}\n\n");

        // Draw triangulation
        for(Finite_edges_iterator eit = T.finite_edges_begin(); eit != T.finite_edges_end(); ++eit) {
            DelaunayTriangulation::Edge e = *eit;
            double x1 = e.first->vertex( (e.second+1)%3 )->point().x();
            double y1 = e.first->vertex( (e.second+1)%3 )->point().y();
            double x2 = e.first->vertex( (e.second+2)%3 )->point().x();
            double y2 = e.first->vertex( (e.second+2)%3 )->point().y();
            fprintf(fp, "\\draw [thin,stroke=gray] (%f,%f) -- (%f,%f);\n", x1,y1, x2,y2);
        }

        // Draw vertices
        for(Finite_vertices_iterator it = T.finite_vertices_begin(); it != T.finite_vertices_end(); ++it)
            fprintf(fp,"\\draw [fill=red,stroke=red] (%f,%f) circle [radius=%f];\n",it->point().x(),it->point().y(),radiusOfPoints);

        // draw spanning graph
        for( auto el : graph ) {
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

}; // class TriangulationPrinter

#endif // GSNUNF_CGALCOMPONENTS_H
