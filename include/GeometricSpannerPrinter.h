#ifndef GSNUNF_GEOMETRICSPANNERPRINTER_H
#define GSNUNF_GEOMETRICSPANNERPRINTER_H

#include "CGALComponents.h"
#include "DelaunayGraph.h"

namespace gsnunf {

class GeometricSpannerPrinter {
  public:

    double radiusOfPoints;

    GeometricSpannerPrinter(double radiusOfPoints = 0.09);

    void draw( const DelaunayTriangulation& T, std::string fName );
    void draw( const DelaunayGraph& graph, std::string fName );

}; // class TriangulationPrinter

} // namespace gsnunf

#endif // GSNUNF_GEOMETRICSPANNERPRINTER_H
