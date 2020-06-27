#ifndef GSNUNF_PLANARSPANNER_H
#define GSNUNF_PLANARSPANNER_H

#include "CGALComponents.h"
#include "Graph.h"
#include "SpanningGraph.h"



namespace gsnunf {

class PlanarSpanner : public Graph {
  public:
    PlanarSpanner( Graph &G, double epsilon );

  protected:

}; // class PlanarSpanner

} // namespace gsnunf

#endif // GSNUNF_PLANARSPANNER_H
