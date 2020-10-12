#ifndef GSNUNF_GRAPHALGOTV_H
#define GSNUNF_GRAPHALGOTV_H

#include <vector>

#include <boost/functional/hash.hpp>
//#include <QApplication>

#include "../qt/graphscene.h"
#include "../qt/mainwindow.h"
#include "utilities.h"

namespace gsnunf {

using namespace std;

enum ItemType { Vertex=0, Edge, ITEM_TYPE_SIZE };

enum EventType { Symbol=0, Focus, Label, TextOut, EVENT_TYPE_SIZE };

enum FocusLevel { First=1, Second, Third, FOCUS_LEVEL_MAX=Third };

struct TV_Event {
    EventType eventType;
    TV_Event( EventType type )
      : eventType(type) {}
};

template< typename T >
struct TV_EventData {
    T value;
    TV_EventData( T value ) : value(value) {}
};

struct TV_ItemEvent : public TV_Event {
    ItemType itemType;
    size_t id;
    TV_ItemEvent( EventType event_t, ItemType item_t, size_t id )
      : TV_Event(event_t), itemType(item_t), id(id) {}
};

struct TV_SymbolEvent : public TV_ItemEvent, TV_EventData< pair<string,bool> > {
    TV_SymbolEvent( ItemType item_t, size_t id, string name, bool status )
      : TV_ItemEvent( EventType::Symbol, item_t, id ),
        TV_EventData( make_pair(name,status) ) {}
};

struct TV_FocusEvent : public TV_ItemEvent, TV_EventData<FocusLevel> {
    TV_FocusEvent( size_t id, FocusLevel level = First )
      : TV_ItemEvent( EventType::Focus, ItemType::Vertex, id ),
        TV_EventData( level ) {}
};

struct TV_LabelEvent : public TV_ItemEvent, TV_EventData<string> {
    TV_LabelEvent( ItemType item_t, size_t id, string label = "" )
      : TV_ItemEvent( EventType::Label, item_t, id ),
        TV_EventData( label ) {}
};

struct TV_ConsoleEvent : public TV_Event, TV_EventData<string> {
    TV_ConsoleEvent( string out )
      : TV_Event( EventType::TextOut ),
        TV_EventData( out ) {}
};

class GraphAlgoTV {
  public:
    using TV_Point = pair<double,double>;
    using TV_Edge = pair<size_t,size_t>;
    using EventQueue = std::queue< TV_Event* >;

    EventQueue _eventQueue;

    // Holds the raw list of events
    std::list< TV_SymbolEvent > _SymbolEvents;
    std::list< TV_FocusEvent > _FocusEvents;
    std::list< TV_LabelEvent > _LabelEvents;
    std::list< TV_ConsoleEvent > _ConsoleEvents;

    std::list< TV_Event* > _MainEventList;

    std::vector< vector< string > > symbolNames; // this one holds the proper ordering
    std::vector< unordered_map<string,string> > symbolColors; // symbols[ItemType][symbol name] = symbol color

    std::vector< TV_Point > V; // vertices
    std::unordered_map< TV_Point, size_t, boost::hash<TV_Point> > V_map;

    std::vector< TV_Edge > E; // edges
    std::unordered_map< TV_Edge, size_t, boost::hash<TV_Edge> > E_map;




    GraphAlgoTV() : symbolNames( ItemType::ITEM_TYPE_SIZE, {"status"} ), symbolColors( ItemType::ITEM_TYPE_SIZE, { make_pair( "status", "default" ) } ) {

    }

    bool registerSymbol( ItemType type, string name, string color = "default" ) {
        bool inserted = false;

        tie(ignore,inserted) = symbolColors.at(type).emplace( name, color );
        if(inserted) symbolNames.at(type).emplace_back(name);

        return inserted;
    }

    bool registerVertex( const TV_Point& point, string label = "" ) {
        bool inserted = false;

        tie(ignore,inserted) = V_map.emplace( point, V.size() );
        if(inserted) V.push_back(point);

        return inserted;
    }

    bool registerEdge( const TV_Point& p, const TV_Point& q, string label = "", bool force = false ) {
        bool inserted = false;

        if( force ) {
            ignore = registerVertex(p);
            ignore = registerVertex(q);
        }

        if( contains( V_map, p ) && contains( V_map, q ) ) {
            TV_Edge e = make_pair( V_map.at(p), V_map.at(q) );
            tie(ignore,inserted) = E_map.emplace( e, E.size() );
            if(inserted) E.push_back(e);
        }
        return inserted;
    }

    template< typename Tri >
    bool registerTriangulation( const Tri& Triangulation, vector<string> labels = {} ) {
        size_t i=0;
        // register vertices
        for( auto v=Triangulation.finite_vertices_begin(); v!=Triangulation.finite_vertices_end(); ++v ) {
            string label = i < labels.size() ? to_string(i) : "";
            registerVertex( make_pair( v->point().x(), v->point().y() ), label );
        }
        // register edges
        for( typename Tri::Finite_edges_iterator e=Triangulation.finite_edges_begin(); e!=Triangulation.finite_edges_end(); ++e ) {

            typename Tri::Face_handle fh = e->first;
            typename Tri::Vertex_handle v0 = fh->vertex(0);
            typename Tri::Vertex_handle v1 = fh->vertex(1);
            double x1 = v0->point().x();
            double y1 = v0->point().y();
            double x2 = v1->point().x();
            double y2 = v1->point().y();
            registerEdge( { x1, y1 }, { x2, y2 } );
        }
    }

    void activate( TV_Point p, string label = "status" ) {

    }

//    // Change status of vertex v
//    void addToEventQueue( Vertex_handle v, bool status ) {
//        VertexStatusEvent<T> add( v, status );
//        _vertexStatusEvents.push_back( add ); // put the event in a container that won't slice it
//        _eventQueue.push( &_vertexStatusEvents.back() );           // put the event's address in the event queue
//    }
//
//    // Change focus[level] to highlight vertex v
//    void addToEventQueue( Vertex_handle v, int level ) {
//        VertexFocusEvent<T> add( v, level );
//        _focusEvents.push_back( add ); // put the event in a container that won't slice it
//        _eventQueue.push( &_focusEvents.back() );          // put the event's address in the event queue
//    }
//
//    // Change status of edge e
//    void addToEventQueue( std::pair<Vertex_handle,Vertex_handle> e, bool status ) {
//        EdgeStatusEvent<T> add( e, status );
//        _edgeStatusEvents.push_back( add ); // put the event in a container that won't slice it
//        _eventQueue.push( &_edgeStatusEvents.back() );           // put the event's address in the event queue
//    }
//
//
//    typename std::list<GraphAlgorithmEvent* >::iterator eventListEnd() {
//        return _eventList.end();
//    }
//
    int play() {
        auto start = processEventQueue();
        QApplication a(argc, argv);
        MainWindow w;
        w.show();
        return a.exec();
    }

    std::list<TV_Event*>::iterator processEventQueue() {
//        TV_Event* e;
//        EventQueue events( _eventQueue );
//        AdjacencyList E;
//        std::vector<Vertex_handle> focus;
//        VertexHash V;
//
//        for( auto it=_DT.finite_vertices_begin(); it!=_DT.finite_vertices_end(); ++it )
//            V.insert(it);
//
//        while( !events.empty() ) {
//            e = events.front();
//
//            switch( e->type ) {
//            // Only add an event that actually makes a change
//            // I.e. current state and event state must be different, use xor
//            case EventType::Vertex: {
//                // cast the event to the proper type, e is a pointer
//                auto event = static_cast< VertexStatusEvent<DelaunayGraph>* >(e);
//                if( contains( V, event->v ) ^ event->active ) {
//                    _eventList.push_back(e);
//                    toggle( V, event->v );
//                }
//            }
//            break;
//            case EventType::Edge: {
//                // cast the event to the proper type, e is a pointer
//                auto event = static_cast< EdgeStatusEvent<DelaunayGraph>* >(e);
//                Vertex_handle v1 = event->e.first,
//                              v2 = event->e.second;
//
//                if( !contains( E, v1 ) )
//                    E.emplace( v1, VertexSet() );
//                auto incident = &E.find(v1)->second;
//
//                if( contains( *incident, v2 ) ^ event->active ) {
//                    _eventList.push_back(e);
//
//                    toggle( *incident, v2 );
//                    if( !contains( E, v2 ) )
//                        E.emplace( v2, VertexSet() );
//                    incident = &E.find(v2)->second;
//                    toggle( *incident, v1 );
//                }
//            }
//            break;
//            case EventType::Focus: {
//                // cast the event to the proper type, e is a pointer
//                auto event = static_cast< VertexFocusEvent<DelaunayGraph>* >(e);
//                if( event->level <= focus.size() ) { // only react to the event if
//                    _eventList.push_back(e);
//
//                    // remove all from e.nextFocus to end
//                    focus.erase( focus.begin()+event->level, focus.end() );
//                    focus.push_back( event->v );
//                }
//            }
//            break;
//            default:
//                std::cout<<"Invalid event type.\n";
//            }
//            events.pop(); // remove first in line
//        }
//        return _eventList.begin();
    }
};

} // namespace gsnunf

#endif // GSNUNF_GRAPHALGOTV_H

