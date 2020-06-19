#ifndef SPANNINGGRAPH_H
#define SPANNINGGRAPH_H

#include "CGALComponents.h"
#include "Graph.h"

#include <list>



namespace gsnunf {

class SpanningGraph : public Graph {
public:

    SpanningGraph() {}
    ~SpanningGraph() {}

    SpanningGraph( std::list<Point> &S );

    void remove_first_edge( Vertex_circulator C );
    void remove_second_edge( Vertex_circulator C );
    void remove_last_edge( Vertex_circulator C );

    void add_first_edge( Vertex_handle v, Vertex_circulator C );
    void add_second_edge( Vertex_handle v, Vertex_circulator C );
    void add_last_edge( Vertex_handle v, Vertex_circulator C );

protected:

private:


};

}

#endif

