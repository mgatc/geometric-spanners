#ifndef GSNUNF_GEOMETRICSPANNERPRINTER_H
#define GSNUNF_GEOMETRICSPANNERPRINTER_H

#include "CGALComponents.h"
#include "Graph.h"

namespace gsnunf {

class GeometricSpannerPrinter {
  public:

    double radiusOfPoints;

    GeometricSpannerPrinter(double radiusOfPoints = 0.09);

    void drawTriangulation(const DelaunayTriangulation &T, std::string fName);
    void drawGraph(Graph &graph, std::string fName);

}; // class TriangulationPrinter

} // namespace gsnunf

#endif // GSNUNF_GEOMETRICSPANNERPRINTER_H
