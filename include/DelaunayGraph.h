#ifndef GSNUNF_GRAPH_H
#define GSNUNF_GRAPH_H

#include <list>
#include <memory>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Vector_2.h>

#include "GraphAlgorithmEvent.h"



namespace gsnunf {

const double PI = 3.14159265359;
const double EPSILON = 0.000001;

template< class T >
bool contains( const T& V, const typename T::key_type& v ) {
    return V.find(v) != V.end();
}

/* If V contains v, remove v.
 * If V does not contain v, add it.
 * Return new value.
 */
template< class T >
bool toggle( T& V, const typename T::key_type& v ) {
    bool found = contains( V, v );
    if( found ) V.erase( v );
    else V.insert( v );
    return !found;
}

template< class T >
class DelaunayGraph {
  public:
    /* Types */
    typedef typename T::Geom_traits K;
    typedef typename T::Vertex_handle Vertex_handle;
    typedef typename T::Vertex_circulator Vertex_circulator;
    typedef typename T::Finite_vertices_iterator Finite_vertices_iterator;
    typedef typename T::Finite_edges_iterator Finite_edges_iterator;
    typedef typename T::Face_handle Face_handle;
    typedef typename K::FT FT;

    typedef CGAL::Vector_2<K> Vector_2;
    typedef CGAL::Container_from_circulator<Vertex_circulator> Vertex_container;

    typedef std::set<Vertex_handle> VertexSet;
    typedef std::unordered_set<Vertex_handle> VertexHash;
    template< typename V >
    using VertexMap = std::unordered_map< Vertex_handle, V >;
    using AdjacencyList = VertexMap< VertexSet >;
    typedef std::queue<GraphAlgorithmEvent* > EventQueue;
    template< typename N >
    using EdgeInfoMap = VertexMap< VertexMap< std::optional<N> > >;

    /* Data */
    const T& _DT;
    AdjacencyList _E;
    EventQueue _eventQueue;

    std::list< VertexStatusEvent<DelaunayGraph> > _vertexStatusEvents;
    std::list< EdgeStatusEvent<DelaunayGraph> > _edgeStatusEvents;
    std::list< VertexFocusEvent<DelaunayGraph> > _focusEvents;
    std::list< GraphAlgorithmEvent* > _eventList;

    /* Functions */
    DelaunayGraph( const T& DT ) : _DT(DT) {}
    DelaunayGraph( const DelaunayGraph<T>& G ) : _DT(G._DT), _E(G._E) {}

    // Change status of vertex v
    void addToEventQueue( Vertex_handle v, bool status ) {
        VertexStatusEvent<DelaunayGraph> add( v, status );
        _vertexStatusEvents.push_back( add ); // put the event in a container that won't slice it
        _eventQueue.push( &_vertexStatusEvents.back() );           // put the event's address in the event queue
    }

    // Change focus[level] to highlight vertex v
    void addToEventQueue( Vertex_handle v, int level ) {
        VertexFocusEvent<DelaunayGraph> add( v, level );
        _focusEvents.push_back( add ); // put the event in a container that won't slice it
        _eventQueue.push( &_focusEvents.back() );          // put the event's address in the event queue
    }

    // Change status of edge e
    void addToEventQueue( std::pair<Vertex_handle,Vertex_handle> e, bool status ) {
        EdgeStatusEvent<DelaunayGraph> add( e, status );
        _edgeStatusEvents.push_back( add ); // put the event in a container that won't slice it
        _eventQueue.push( &_edgeStatusEvents.back() );           // put the event's address in the event queue
    }

    void add_edge( const Vertex_handle v1, const Vertex_handle v2 ) {
        add_half_edge( v1, v2 );
        add_half_edge( v2, v1 );
        addToEventQueue( std::make_pair(v1,v2), true );
    }

    void remove_edge( const Vertex_handle v1, const Vertex_handle v2 ) {
        remove_half_edge( v1, v2 );
        remove_half_edge( v2, v1 );
        addToEventQueue( std::make_pair(v1,v2), false );
    }

    void add_all_edges() {
        for( auto eit = _DT.finite_edges_begin(); eit != _DT.finite_edges_end(); ++eit ) {
            typename T::Edge e = *eit;
            Vertex_handle p = e.first->vertex( (e.second+1)%3 ),
                          q = e.first->vertex( (e.second+2)%3 );
            add_edge( p, q );
        }
    }

    double get_angle( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
        // TODO: use math/angle from CGAL, not stl
        assert( !_DT.is_infinite(p) && !_DT.is_infinite(q) && !_DT.is_infinite(r) );

        Vector_2 pq( p->point(), q->point() );
        Vector_2 rq( r->point(), q->point() );

        double result = atan2( rq.y(), rq.x()) - atan2(pq.y(), pq.x() );

        // atan() returns a value between -PI and PI. From zero ("up"), CCW rotation is negative and CW is positive.
        // Our zero is also "up," but we only want positive values between 0 and 2*PI:

        result *= -1; // First, invert the result. This will associate CW rotation with positive values.
        if( result < EPSILON ) // Then, if the result is less than 0 (or epsilon for floats) add 2*PI.
            result += 2*PI;

        return result;
    }

    /**
     *  Given a Delaunay Triangulation DT and an output list out, compute the canonical ordering of
     *  the underlying point set.
     */
    void canonical_order( std::list<Vertex_handle> &out ) {
        out.clear();

        if( _DT.number_of_vertices() <= 3 ) {
            for( auto v=_DT.finite_vertices_begin(); v!=_DT.finite_vertices_end(); ++v ) {
                out.push_back( v->handle() );
                addToEventQueue( v, 0 );// focus0 on v->handle
                addToEventQueue( v, false ); // inactivate v->handle
            }
            return;
        }

        VertexHash eligible_vertices, empty_set, has_incident_chords, is_removed, on_outer_face, reserved_vertices;
        Vertex_handle v_1,v_2,v_k;

        Vertex_circulator v_convex_hull = _DT.incident_vertices( _DT.infinite_vertex() ), // create a circulator of the convex hull
                          done( v_convex_hull );

        do { // cycle through convex hull to set vertex info
            on_outer_face.insert( v_convex_hull );  // incident_chords (called in next loop) relies on accurate on_outer_face values
            addToEventQueue( v_convex_hull, 0 );
        } while( ++v_convex_hull != done );

        do { // cycle through convex hull again to add eligible vertices to list
            update_incident_chords( v_convex_hull, has_incident_chords, eligible_vertices, is_removed, on_outer_face, reserved_vertices, _DT.number_of_vertices(), false );
            update_eligibility( v_convex_hull, eligible_vertices, reserved_vertices, has_incident_chords );
            addToEventQueue( v_convex_hull, 0 );
        } while( ++v_convex_hull != done );

        // find two consecutive valid vertices to be v_1 and v_2
        while( !contains( eligible_vertices, v_convex_hull ) && !contains( eligible_vertices, ++v_convex_hull ) ); // loop until reaching valid vertex

        v_1 = v_convex_hull;
        addToEventQueue( v_1, 0 );

        v_2 = --v_convex_hull;
        addToEventQueue( v_2, 0 );

        // v_1 and v_2 must be reserved as the last two vertices to be added
        reserved_vertices.insert(v_1);
        reserved_vertices.insert(v_2);
        eligible_vertices.erase(v_1);
        eligible_vertices.erase(v_2);

        while( !eligible_vertices.empty() ) {

            v_k = *eligible_vertices.begin();
            addToEventQueue( v_k, 0 );
            addToEventQueue( v_k, false );

            eligible_vertices.erase(v_k);
            is_removed.insert(v_k);
            on_outer_face.erase(v_k);
            out.push_front(v_k);

            // Update neighbors
            Vertex_circulator v_n = _DT.incident_vertices( v_k ),
                              done( v_n );
            do { // update outer face
                addToEventQueue( v_n, 1 );
                if( !_DT.is_infinite(v_n) && !contains( is_removed, v_n ) )
                    on_outer_face.insert( v_n );
            } while( ++v_n != done );

            do { // test incidence of new outer face vertices
                addToEventQueue( v_n, 1 );
                if( !_DT.is_infinite(v_n) && !contains( is_removed, v_n ) && !contains( reserved_vertices, v_n ) ) {
                    update_incident_chords( v_n, has_incident_chords, eligible_vertices, is_removed, on_outer_face, reserved_vertices, _DT.number_of_vertices()-out.size(), true );
                    update_eligibility( v_n, eligible_vertices, reserved_vertices, has_incident_chords );
                }
            } while( ++v_n != done );
        }
        out.push_front(v_2);
        addToEventQueue( v_2, 0 );
        addToEventQueue( v_2, false );
        out.push_front(v_1);
        addToEventQueue( v_1, 0 );
        addToEventQueue( v_1, false );

        assert( _DT.number_of_vertices() == out.size() );
    }

    /**
     *  Given a Delaunay Triangulation DT, a vertex handle v_update, an unordered set of eligible
     *  vertices, an unordered set of reserved vertices, the number of vertices left in DT,
     *  and a boolean value for controlling recursive behavior, update v_update's incident chord
     *  count. Additionally, if update_neighbors is true, recursive call update_incident_chords on
     *  each neighbor.
     */
    void update_incident_chords( Vertex_handle v_update,
                                 VertexHash& has_incident_chords,
                                 VertexHash& eligible_vertices,
                                 const VertexHash& is_removed,
                                 const VertexHash& on_outer_face,
                                 const VertexHash& reserved_vertices,
                                 int vertices_remaining,
                                 bool update_neighbors ) {

        // Conditions that will result in 0 incident chords to v_update
        if(   _DT.is_infinite( v_update )            // Vertex is the infinite vertex
           || contains( is_removed, v_update )       // Vertex is removed already
           || !contains( on_outer_face, v_update )   // If a vertex is not on the outer face, it cannot be part of a chord
           || vertices_remaining <= 3 ) {            // Graphs must contain at least 4 vertices to contain a chord

            has_incident_chords.erase( v_update );
            return;
        }

        // Otherwise, v_update is guaranteed to be incident to 2 other vertices with on_outer_face==true.
        // Any additional vertex incident to v_update with on_outer_face==true forms a chord with v_update.
        // These vertices must be updated too.

        int incident_chords = -2;

        Vertex_circulator v_neighbor = _DT.incident_vertices( v_update ),
                                       done( v_neighbor );
        do {
            // focus2 on v_neighbor
            if( !_DT.is_infinite( v_neighbor ) && !contains( is_removed, v_neighbor ) && contains( on_outer_face, v_neighbor ) ) {
                if( update_neighbors )
                    update_incident_chords( v_neighbor, has_incident_chords, eligible_vertices, is_removed, on_outer_face, reserved_vertices, vertices_remaining, false );

                incident_chords++;

                update_eligibility( v_neighbor, eligible_vertices, reserved_vertices, has_incident_chords );
            }
        } while( ++v_neighbor != done );

        assert( incident_chords >= 0 ); // Make sure our assumptions are correct

        if( incident_chords == 0 )
            has_incident_chords.erase( v_update );
        else
            has_incident_chords.insert( v_update );
    }

    DelaunayGraph<T>& operator=( DelaunayGraph<T>&& other ) { // move assignment
        if(&other != *this) {
            _DT = other._DT;
            _E = other._E;
        }
        return *this;
    }

    int count_valid_neighbors( Vertex_circulator C, const VertexHash& invalid ) {
        Vertex_circulator done(C);
        int k = 0; // count N_i path length k

        do {
            // ++focus on C
            if( !contains( invalid, C ) && !_DT.is_infinite(C) )
                k++;
        } while( ++C != done );

        return k;
    }

    void normalize_circulator( Vertex_circulator &C, const VertexHash& invalid ) {
        Vertex_circulator done = C;
        // Position circulator so that we are guaranteed to be on the first vertex on the path N_i
        // First, loop until the circulator reaches an invalid vertex or completes a full rotation
        while( !contains( invalid, C ) && !_DT.is_infinite(C) && ++C != done );// cout<<v_n->point()<<"\n";
        done = C;
        // Loop until the circulator reaches a valid vertex
        while( ( contains( invalid, C ) || _DT.is_infinite(C) ) && ++C != done );// cout<<v_n->point()<<"\n";
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



  private:

    void add_half_edge( typename T::Vertex_handle v1, typename T::Vertex_handle v2 ) {
        typename AdjacencyList::iterator v1_incident_vertices = _E.find(v1);

        if( v1_incident_vertices == _E.end() ) // v1 not found in g
            std::tie( v1_incident_vertices, std::ignore ) = _E.insert( make_pair( v1, VertexSet() ) );

        v1_incident_vertices->second.insert(v2);
    }

    void remove_half_edge( Vertex_handle v1, Vertex_handle v2 ) {

        typename AdjacencyList::iterator v1_incident_vertices = _E.find(v1);

        // if v1 in present in _E, and v2 is present in v1, remove v2 from v1
        if( contains( _E, v1 ) && contains( v1_incident_vertices->second, v2 ) )
            v1_incident_vertices->second.erase(v2);

    }

    /**
     *  Given an unordered set of eligible vertices, a vertex handle v, and an unordered set of
     *  reserved vertices, perform any necessary insertion or removal of v from eligible
     *  based on its current presence in eligible, its incident chord status, and its reserved status.
     */
    void update_eligibility( const Vertex_handle v, VertexHash& eligible, const VertexHash& reserved, const VertexHash& has_incident_chords ) {
        bool is_reserved = contains( reserved, v );
        bool is_in_eligible = contains( eligible, v );
        bool does_have_incident_chords = contains( has_incident_chords, v );

        if( !does_have_incident_chords && !is_reserved && !is_in_eligible )
            eligible.insert(v);

        if( ( does_have_incident_chords || is_reserved ) && is_in_eligible )
            eligible.erase(v);
    }

}; // class DelaunayGraph

} // namespace gsnunf

#endif // GSNUNF_GRAPH_H
