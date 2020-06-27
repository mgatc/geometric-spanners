#ifndef GSNUNF_SPANNINGGRAPH_H
#define GSNUNF_SPANNINGGRAPH_H

#include <list>

#include "CGALComponents.h"
#include "Graph.h"



namespace gsnunf {

class SpanningGraph : public Graph {
  public:

    SpanningGraph() {}
    ~SpanningGraph() {}

    SpanningGraph( const std::list<Point> &S );
    SpanningGraph( Graph& G );

  protected:

    void add_first_edge( Vertex_handle v, Vertex_circulator C );
    void add_second_edge( Vertex_handle v, Vertex_circulator C );
    void add_last_edge( Vertex_handle v, Vertex_circulator C );

    void compute_spanning_graph();

    void remove_first_edge( Vertex_circulator C );
    void remove_second_edge( Vertex_circulator C );
    void remove_last_edge( Vertex_circulator C );


  private:


}; // class SpanningGraph

} // namespace gsnunf

#endif // GSNUNF_SPANNINGGRAPH_H
