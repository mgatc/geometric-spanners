#ifndef GSNUNF_PLANARSPANNER_H
#define GSNUNF_PLANARSPANNER_H

#include "CGALComponents.h"
#include "DelaunayGraph.h"
#include "SpanningGraph.h"



namespace gsnunf {

class PlanarSpanner : public DelaunayGraph {
  public:
    PlanarSpanner( shared_ptr<DelaunayTriangulation> DT, double epsilon );

  protected:

}; // class PlanarSpanner

} // namespace gsnunf

#endif // GSNUNF_PLANARSPANNER_H
