#ifndef GSNUNF_GRAPHALGOTV_H
#define GSNUNF_GRAPHALGOTV_H

#include <vector>

namespace gsnunf {

using namespace std;

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

template< class T >
class GraphAlgoTV {
public:
    using Vertex_handle = typename T::Vertex_handle;
    using EventQueue = std::queue< GraphAlgorithmEvent* >;

    EventQueue _eventQueue;

    std::list< VertexStatusEvent<T> > _vertexStatusEvents;
    std::list< EdgeStatusEvent<T> > _edgeStatusEvents;
    std::list< VertexFocusEvent<T> > _focusEvents;
    std::list< GraphAlgorithmEvent* > _eventList;

    // Change status of vertex v
    void addToEventQueue( Vertex_handle v, bool status ) {
        VertexStatusEvent<T> add( v, status );
        _vertexStatusEvents.push_back( add ); // put the event in a container that won't slice it
        _eventQueue.push( &_vertexStatusEvents.back() );           // put the event's address in the event queue
    }

    // Change focus[level] to highlight vertex v
    void addToEventQueue( Vertex_handle v, int level ) {
        VertexFocusEvent<T> add( v, level );
        _focusEvents.push_back( add ); // put the event in a container that won't slice it
        _eventQueue.push( &_focusEvents.back() );          // put the event's address in the event queue
    }

    // Change status of edge e
    void addToEventQueue( std::pair<Vertex_handle,Vertex_handle> e, bool status ) {
        EdgeStatusEvent<T> add( e, status );
        _edgeStatusEvents.push_back( add ); // put the event in a container that won't slice it
        _eventQueue.push( &_edgeStatusEvents.back() );           // put the event's address in the event queue
    }


    typename std::list<GraphAlgorithmEvent* >::iterator eventListEnd() {
        return _eventList.end();
    }

    typename std::list<GraphAlgorithmEvent* >::iterator processEventQueue() {
        GraphAlgorithmEvent* e;
        EventQueue events( _eventQueue );
        AdjacencyList E;
        std::vector<Vertex_handle> focus;
        VertexHash V;

        for( auto it=_DT.finite_vertices_begin(); it!=_DT.finite_vertices_end(); ++it )
            V.insert(it);

        while( !events.empty() ) {
            e = events.front();

            switch( e->type ) {
            // Only add an event that actually makes a change
            // I.e. current state and event state must be different, use xor
            case EventType::Vertex: {
                // cast the event to the proper type, e is a pointer
                auto event = static_cast< VertexStatusEvent<DelaunayGraph>* >(e);
                if( contains( V, event->v ) ^ event->active ) {
                    _eventList.push_back(e);
                    toggle( V, event->v );
                }
            }
            break;
            case EventType::Edge: {
                // cast the event to the proper type, e is a pointer
                auto event = static_cast< EdgeStatusEvent<DelaunayGraph>* >(e);
                Vertex_handle v1 = event->e.first,
                              v2 = event->e.second;

                if( !contains( E, v1 ) )
                    E.emplace( v1, VertexSet() );
                auto incident = &E.find(v1)->second;

                if( contains( *incident, v2 ) ^ event->active ) {
                    _eventList.push_back(e);

                    toggle( *incident, v2 );
                    if( !contains( E, v2 ) )
                        E.emplace( v2, VertexSet() );
                    incident = &E.find(v2)->second;
                    toggle( *incident, v1 );
                }
            }
            break;
            case EventType::Focus: {
                // cast the event to the proper type, e is a pointer
                auto event = static_cast< VertexFocusEvent<DelaunayGraph>* >(e);
                if( event->level <= focus.size() ) { // only react to the event if
                    _eventList.push_back(e);

                    // remove all from e.nextFocus to end
                    focus.erase( focus.begin()+event->level, focus.end() );
                    focus.push_back( event->v );
                }
            }
            break;
            default:
                std::cout<<"Invalid event type.\n";
            }
            events.pop(); // remove first in line
        }
        return _eventList.begin();
    }
};

} // namespace gsnunf

#endif // GSNUNF_GRAPHALGOTV_H

