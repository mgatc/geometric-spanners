#ifndef ALGORITHMEVENT_H
#define ALGORITHMEVENT_H

//#include "DelaunayGraph.h"

template< class T >
struct GraphAlgorithmEvent {
    typename T::Vertex_handle v1;
    bool active;
};

template< class T >
struct VertexEvent : public GraphAlgorithmEvent {
    int focus;
};

template< class T >
struct EdgeEvent : public GraphAlgorithmEvent {
    typename T::Vertex_handle v2;
};

#endif // ALGORITHMEVENT_H
