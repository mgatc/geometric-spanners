#ifndef GSNUNF_GRAPH_H
#define GSNUNF_GRAPH_H

#include <list>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <CGAL/circulator.h>

#include "CGALComponents.h"



namespace gsnunf {

template<typename T>
bool contains( const T& V, const typename T::key_type& v ) {
    return V.find(v) != V.end();
}

template< class T >
class DelaunayGraph {
  public:

    typedef typename T::Vertex_handle Vertex_handle;
    typedef typename T::Vertex_circulator Vertex_circulator;
    typedef typename T::Finite_vertices_iterator Finite_vertices_iterator;
    typedef typename T::Finite_edges_iterator Finite_edges_iterator;
    typedef typename T::Face_handle Face_handle;

    typedef CGAL::Vector_2<typename T::Geom_traits> Vector_2;
    typedef CGAL::Container_from_circulator<Vertex_circulator> Vertex_container;

    typedef std::set<Vertex_handle> VertexSet;
    typedef std::unordered_set<Vertex_handle> VertexHash;

    template <typename V>
    using VertexMap = std::unordered_map<Vertex_handle, V>;

    using AdjacencyList = VertexMap<VertexSet>;



    const T& _DT;
    AdjacencyList _E;



    DelaunayGraph( const T& DT ) : _DT(DT) {}
    DelaunayGraph( const DelaunayGraph<T>& G ) : _DT(G._DT), _E(G._E) {}



    void add_edge( const Vertex_handle v1, const Vertex_handle v2 ) {
        add_half_edge( v1, v2 );
        add_half_edge( v2, v1 );
        // activate edge( v1, v2 )
    }

    void remove_edge( const Vertex_handle v1, const Vertex_handle v2 ) {
        remove_half_edge( v1, v2 );
        remove_half_edge( v2, v1 );
        // inactivate edge( v1, v2 )
    }

    double get_angle( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
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
                // focus0 on v->handle
                // inactivate v->handle
            }
            return;
        }

        VertexHash eligible_vertices, empty_set, has_incident_chords, is_removed, on_outer_face, reserved_vertices;
        Vertex_handle v_1,v_2,v_k;

        Vertex_circulator v_convex_hull = _DT.incident_vertices( _DT.infinite_vertex() ), // create a circulator of the convex hull
                          done( v_convex_hull );

        do { // cycle through convex hull to set vertex info
            on_outer_face.insert( v_convex_hull );  // incident_chords (called in next loop) relies on accurate on_outer_face values
            // focus0 on v_convex_hull
        } while( ++v_convex_hull != done );

        do { // cycle through convex hull again to add eligible vertices to list
            update_incident_chords( v_convex_hull, has_incident_chords, eligible_vertices, is_removed, on_outer_face, reserved_vertices, _DT.number_of_vertices(), false );
            update_eligibility( v_convex_hull, eligible_vertices, reserved_vertices, has_incident_chords );
        } while( ++v_convex_hull != done );

        // find two consecutive valid vertices to be v_1 and v_2
        while( !contains( eligible_vertices, v_convex_hull ) && !contains( eligible_vertices, ++v_convex_hull ) ); // loop until reaching valid vertex
        v_1 = v_convex_hull;
        // focus0 on v_1
        v_2 = --v_convex_hull;
        // focus0 on v_2

        // v_1 and v_2 must be reserved as the last two vertices to be added
        reserved_vertices.insert(v_1);
        reserved_vertices.insert(v_2);
        eligible_vertices.erase(v_1);
        eligible_vertices.erase(v_2);

        while( !eligible_vertices.empty() ) {

            v_k = *eligible_vertices.begin();
            // focus0 on v_k

            // inactivate v_k

            eligible_vertices.erase(v_k);
            is_removed.insert(v_k);
            on_outer_face.erase(v_k);
            out.push_front(v_k);

            // Update neighbors
            Vertex_circulator v_n = _DT.incident_vertices( v_k ),
                              done( v_n );
            do { // update outer face
                // focus1 on v_n
                if( !_DT.is_infinite(v_n) && !contains( is_removed, v_n ) )
                    on_outer_face.insert( v_n );
            } while( ++v_n != done );

            do { // test incidence of new outer face vertices
                if( !_DT.is_infinite(v_n) && !contains( is_removed, v_n ) && !contains( reserved_vertices, v_n ) ) {
                    update_incident_chords( v_n, has_incident_chords, eligible_vertices, is_removed, on_outer_face, reserved_vertices, _DT.number_of_vertices()-out.size(), true );
                    update_eligibility( v_n, eligible_vertices, reserved_vertices, has_incident_chords );
                }
            } while( ++v_n != done );
        }
        out.push_front(v_2);
        // focus0 on v_2
        // inactivate v_2
        out.push_front(v_1);
        // focus0 on v_1
        // inactivate v_1

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
