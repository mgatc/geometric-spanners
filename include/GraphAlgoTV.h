#ifndef GSNUNF_GRAPHALGOTV_H
#define GSNUNF_GRAPHALGOTV_H
//
//#include <vector>
//
//#include <boost/functional/hash.hpp>
//
//#include "utilities.h"
//
//namespace unf_spanners {
//
//using namespace std;
//
//enum ItemType { Vertex=0, Edge, ITEM_TYPE_SIZE };
//
//enum EventType { Symbol=0, Focus, Label, TextOut, EVENT_TYPE_SIZE };
//
//enum FocusLevel { First=1, Second, Third, FOCUS_LEVEL_MAX=Third };
//
//template< typename T >
//struct TV_Event {
//    EventType eventType;
//    T value;
//    TV_Event( EventType type, T value )
//      : eventType(type), value(value) {}
//};
//
//template< typename T >
//struct TV_ItemEvent : public TV_Event<T> {
//    ItemType itemType;
//    size_t id;
//    TV_ItemEvent( EventType event_t, T value, ItemType item_t, size_t id )
//      : TV_Event<T>(event_t,value), itemType(item_t), id(id) {}
//};
//
//struct TV_SymbolEvent : public TV_ItemEvent< pair<string,bool> > {
//    TV_SymbolEvent( ItemType item_t, size_t id, string name, bool status )
//      : TV_ItemEvent( EventType::Symbol, make_pair(name,status), item_t, id ) {}
//};
//
//struct TV_FocusEvent : public TV_ItemEvent<FocusLevel> {
//    TV_FocusEvent( size_t id, FocusLevel level = First )
//      : TV_ItemEvent( EventType::Focus, level, ItemType::Vertex, id ) {}
//};
//
//struct TV_LabelEvent : public TV_ItemEvent<string> {
//    TV_LabelEvent( ItemType item_t, size_t id, string label = "" )
//      : TV_ItemEvent( EventType::Label, label, item_t, id ) {}
//};
//
//struct TV_ConsoleEvent : public TV_Event<string> {
//    TV_ConsoleEvent( string out )
//      : TV_Event( EventType::TextOut, out ) {}
//};
//
//
////template< class P >
////struct VertexEvent {
////    P index;
////    VertexEvent( P v ) : v(v) {}
////};
////
////template< class P >
////struct EdgeEvent {
////    std::pair<
////        typename T::VertexHandle,
////        typename T::VertexHandle
////    > e;
////    EdgeEvent( std::pair<typename T::VertexHandle,typename T::VertexHandle> e ) : e(e) {}
////};
////
////template< class T >
////struct VertexStatusEvent : public GraphAlgorithmEvent, VertexEvent<T>, StatusEvent {
////    VertexStatusEvent<T>( typename T::VertexHandle vertex, bool status )
////        : GraphAlgorithmEvent(EventType::Vertex), VertexEvent<T>(vertex), StatusEvent(status) {}
////    VertexStatusEvent<T>( const VertexStatusEvent& other )
////        : GraphAlgorithmEvent(EventType::Vertex), VertexEvent<T>( other.v ), StatusEvent( other.active ) {}
////};
////
////template<class T>
////struct EdgeStatusEvent : public GraphAlgorithmEvent, EdgeEvent<T>, StatusEvent {
////    EdgeStatusEvent<T>( std::pair<typename T::VertexHandle,typename T::VertexHandle> edge, bool status )
////        : GraphAlgorithmEvent(EventType::Edge), EdgeEvent<T>(edge), StatusEvent(status) {}
////    EdgeStatusEvent<T>( const EdgeStatusEvent& other )
////        : GraphAlgorithmEvent( EventType::Edge ), EdgeEvent<T>( other.e ), StatusEvent( other.active ) {}
////};
////
////template< class T >
////struct VertexFocusEvent : public GraphAlgorithmEvent, VertexEvent<T>, FocusEvent {
////    VertexFocusEvent<T>( typename T::VertexHandle vertex, size_t lvl )
////        : GraphAlgorithmEvent(EventType::Focus), VertexEvent<T>(vertex), FocusEvent(lvl) {}
////    VertexFocusEvent<T>( const VertexFocusEvent& other )
////        : GraphAlgorithmEvent(EventType::Focus), VertexEvent<T>( other.v ), FocusEvent( other.level ) {}
////};
//
//class GraphAlgoTV {
//  public:
//    using TV_Point = pair<double,double>;
//    using TV_Edge = pair<size_t,size_t>;
//    using EventQueue = std::queue< TV_Event<void>* >;
//
//    EventQueue _eventQueue;
//
//    // Holds the raw list of events
//    std::list< TV_SymbolEvent > _SymbolEvents;
//    std::list< TV_FocusEvent > _FocusEvents;
//    std::list< TV_LabelEvent > _LabelEvents;
//    std::list< TV_ConsoleEvent > _ConsoleEvents;
//
//    std::list< TV_Event<void>* > _MainEventList;
//
//    std::vector< vector< string > > symbolNames; // this one holds the proper ordering
//    std::vector< unordered_map<string,string> > symbolColors; // symbols[ItemType][symbol name] = symbol color
//
//    std::vector< TV_Point > V; // vertices
//    std::unordered_map< TV_Point, size_t, boost::hash<TV_Point> > V_map;
//
//    std::vector< TV_Edge > E; // edges
//    std::unordered_map< TV_Edge, size_t, boost::hash<TV_Edge> > E_map;
//
//
//    GraphAlgoTV() : symbolNames( ItemType::ITEM_TYPE_SIZE, {"status"} ), symbolColors( ItemType::ITEM_TYPE_SIZE, { make_pair( "status", "default" ) } ) {
//
//    }
//
//    bool registerSymbol( ItemType type, string name, string color = "default" ) {
//        bool inserted = false;
//
//        tie(ignore,inserted) = symbolColors.at(type).emplace( name, color );
//        if(inserted) symbolNames.at(type).emplace_back(name);
//
//        return inserted;
//    }
//
//    bool registerVertex( const TV_Point& point, string label = "" ) {
//        bool inserted = false;
//
//        tie(ignore,inserted) = V_map.emplace( point, V.size() );
//        if(inserted) V.push_back(point);
//
//        return inserted;
//    }
//
//    bool registerEdge( const TV_Point& p, const TV_Point& q, string label = "", bool force = false ) {
//        bool inserted = false;
//
//        if( force ) {
//            ignore = registerVertex(p);
//            ignore = registerVertex(q);
//        }
//
//        if( contains( V_map, p ) && contains( V_map, q ) ) {
//            TV_Edge e = make_pair( V_map.at(p), V_map.at(q) );
//            tie(ignore,inserted) = E_map.emplace( e, E.size() );
//            if(inserted) E.push_back(e);
//        }
//        return inserted;
//    }
//
//    template< typename Tri >
//    bool registerTriangulation( const Tri& Triangulation, vector<string> labels = {} ) {
//        size_t i=0;
//        // register vertices
//        for( auto v=Triangulation.finite_vertices_begin(); v!=Triangulation.finite_vertices_end(); ++v ) {
//            string label = i < labels.size() ? to_string(i) : "";
//            registerVertex( { v->point.x(), v->point.y() }, label );
//        }
//        // register edges
//        for( auto e=Triangulation.finite_edges_begin(); e!=Triangulation.finite_edges_end(); ++e ) {
//            double x1 = e.first->vertex( (e.second+1)%3 )->point().x();
//            double y1 = e.first->vertex( (e.second+1)%3 )->point().y();
//            double x2 = e.first->vertex( (e.second+2)%3 )->point().x();
//            double y2 = e.first->vertex( (e.second+2)%3 )->point().y();
//            registerEdge( { x1, y1 }, { x2, y2 } );
//        }
//    }
//
//    void activate( TV_Point p, string label = "status" ) {
//
//    }
//
////    // Change status of vertex v
////    void addToEventQueue( VertexHandle v, bool status ) {
////        VertexStatusEvent<T> add( v, status );
////        _vertexStatusEvents.push_back( add ); // put the event in a container that won'stretchFactor slice it
////        _eventQueue.push( &_vertexStatusEvents.back() );           // put the event's address in the event queue
////    }
////
////    // Change focus[level] to highlight vertex v
////    void addToEventQueue( VertexHandle v, int level ) {
////        VertexFocusEvent<T> add( v, level );
////        _focusEvents.push_back( add ); // put the event in a container that won'stretchFactor slice it
////        _eventQueue.push( &_focusEvents.back() );          // put the event's address in the event queue
////    }
////
////    // Change status of edge e
////    void addToEventQueue( std::pair<VertexHandle,VertexHandle> e, bool status ) {
////        EdgeStatusEvent<T> add( e, status );
////        _edgeStatusEvents.push_back( add ); // put the event in a container that won'stretchFactor slice it
////        _eventQueue.push( &_edgeStatusEvents.back() );           // put the event's address in the event queue
////    }
////
////
////    typename std::list<GraphAlgorithmEvent* >::iterator eventListEnd() {
////        return _eventList.end();
////    }
////
////    typename std::list<GraphAlgorithmEvent* >::iterator processEventQueue() {
////        GraphAlgorithmEvent* e;
////        EventQueue events( _eventQueue );
////        AdjacencyList E;
////        std::vector<VertexHandle> focus;
////        VertexHash V;
////
////        for( auto it=_DT.finite_vertices_begin(); it!=_DT.finite_vertices_end(); ++it )
////            V.insert(it);
////
////        while( !events.empty() ) {
////            e = events.front();
////
////            switch( e->type ) {
////            // Only add an event that actually makes a change
////            // I.e. current state and event state must be different, use xor
////            case EventType::Vertex: {
////                // cast the event to the proper type, e is a pointer
////                auto event = static_cast< VertexStatusEvent<DelaunayGraph>* >(e);
////                if( contains( V, event->v ) ^ event->active ) {
////                    _eventList.push_back(e);
////                    toggle( V, event->v );
////                }
////            }
////            break;
////            case EventType::Edge: {
////                // cast the event to the proper type, e is a pointer
////                auto event = static_cast< EdgeStatusEvent<DelaunayGraph>* >(e);
////                VertexHandle v1 = event->e.first,
////                              v2 = event->e.second;
////
////                if( !contains( E, v1 ) )
////                    E.emplace( v1, VertexSet() );
////                auto incident = &E.find(v1)->second;
////
////                if( contains( *incident, v2 ) ^ event->active ) {
////                    _eventList.push_back(e);
////
////                    toggle( *incident, v2 );
////                    if( !contains( E, v2 ) )
////                        E.emplace( v2, VertexSet() );
////                    incident = &E.find(v2)->second;
////                    toggle( *incident, v1 );
////                }
////            }
////            break;
////            case EventType::Focus: {
////                // cast the event to the proper type, e is a pointer
////                auto event = static_cast< VertexFocusEvent<DelaunayGraph>* >(e);
////                if( event->level <= focus.size() ) { // only react to the event if
////                    _eventList.push_back(e);
////
////                    // remove all from e.nextFocus to end
////                    focus.erase( focus.begin()+event->level, focus.end() );
////                    focus.push_back( event->v );
////                }
////            }
////            break;
////            default:
////                std::cout<<"Invalid event type.\n";
////            }
////            events.pop(); // remove first in line
////        }
////        return _eventList.begin();
////    }
//};
//
//} // namespace unf_spanners

#endif // GSNUNF_GRAPHALGOTV_H

