#ifndef GRAPHALGORITHMEVENT_H
#define GRAPHALGORITHMEVENT_H

#include <vector>

enum EventType { Vertex, Edge, Focus };

struct GraphAlgorithmEvent {
    EventType type;
    GraphAlgorithmEvent( EventType type )
        : type( type ) {}
};

struct StatusEvent {
    bool active;
    StatusEvent( bool active ) : active(active) {}
};

struct FocusEvent {
    size_t level;
    FocusEvent( size_t level ) : level(level) {}
};

template< class T >
struct VertexEvent {
    typename T::Vertex_handle v;
    VertexEvent( typename T::Vertex_handle v ) : v(v) {}
};

template< class T >
struct EdgeEvent {
    std::pair<
        typename T::Vertex_handle,
        typename T::Vertex_handle
    > e;
    EdgeEvent( std::pair<typename T::Vertex_handle,typename T::Vertex_handle> e ) : e(e) {}
};

template< class T >
struct VertexStatusEvent : public GraphAlgorithmEvent, VertexEvent<T>, StatusEvent {
    VertexStatusEvent<T>( typename T::Vertex_handle vertex, bool status )
        : GraphAlgorithmEvent(EventType::Vertex), VertexEvent<T>(vertex), StatusEvent(status) {}
    VertexStatusEvent<T>( const VertexStatusEvent& other )
        : GraphAlgorithmEvent(EventType::Vertex), VertexEvent<T>( other.v ), StatusEvent( other.active ) {}
};

template<class T>
struct EdgeStatusEvent : public GraphAlgorithmEvent, EdgeEvent<T>, StatusEvent {
    EdgeStatusEvent<T>( std::pair<typename T::Vertex_handle,typename T::Vertex_handle> edge, bool status )
        : GraphAlgorithmEvent(EventType::Edge), EdgeEvent<T>(edge), StatusEvent(status) {}
    EdgeStatusEvent<T>( const EdgeStatusEvent& other )
        : GraphAlgorithmEvent( EventType::Edge ), EdgeEvent<T>( other.e ), StatusEvent( other.active ) {}
};

template< class T >
struct VertexFocusEvent : public GraphAlgorithmEvent, VertexEvent<T>, FocusEvent {
    VertexFocusEvent<T>( typename T::Vertex_handle vertex, size_t lvl )
        : GraphAlgorithmEvent(EventType::Focus), VertexEvent<T>(vertex), FocusEvent(lvl) {}
    VertexFocusEvent<T>( const VertexFocusEvent& other )
        : GraphAlgorithmEvent(EventType::Focus), VertexEvent<T>( other.v ), FocusEvent( other.level ) {}
};


#endif // GRAPHALGORITHMEVENT_H
