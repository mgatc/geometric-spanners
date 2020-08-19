#ifndef GSNUNF_DELAUNAYGRAPH_H
#define GSNUNF_DELAUNAYGRAPH_H

#include <iostream>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Vector_2.h>

#include "GraphAlgorithmEvent.h"
#include "Timer.h"

namespace gsnunf {

using namespace std;

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

class DelaunayGraph {
  public:
    /* Types */
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation_2;
    typedef Delaunay_triangulation_2::Vertex_handle Vertex_handle;
    typedef Delaunay_triangulation_2::Vertex_circulator Vertex_circulator;
    typedef Delaunay_triangulation_2::Finite_vertices_iterator Finite_vertices_iterator;
    typedef Delaunay_triangulation_2::Finite_edges_iterator Finite_edges_iterator;
    typedef Delaunay_triangulation_2::Face_handle Face_handle;
    typedef K::FT FT;
    typedef CGAL::Vector_2<K> Vector_2;
    typedef CGAL::Container_from_circulator<Vertex_circulator> Vertex_container;

    typedef set<Vertex_handle> VertexSet;
    typedef unordered_set<Vertex_handle> VertexHash;
    template< typename V >
    using VertexMap = unordered_map< Vertex_handle, V >;
    using AdjacencyList = VertexMap< VertexSet >;
    template< typename N >
    using EdgeInfoMap = VertexMap< VertexMap< optional<N> > >;



    /* Data */
    Delaunay_triangulation_2 _DT;
    AdjacencyList _E;

    //GraphAlgoTV _algoTV;



    /* Functions */
    template< typename RandomAccessIterator >
    DelaunayGraph( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd ) : _DT( pointsBegin, pointsEnd ) {}

    void add_edge( const Vertex_handle& v1, const Vertex_handle& v2 ) {
        add_half_edge( v1, v2 );
        add_half_edge( v2, v1 );

        //_algoTV.addToEventQueue( std::make_pair(v1,v2), true );
    }

    void remove_edge( const Vertex_handle& v1, const Vertex_handle& v2 ) {
        remove_half_edge( v1, v2 );
        remove_half_edge( v2, v1 );
        //_algoTV.addToEventQueue( std::make_pair(v1,v2), false );
    }

    void add_all_edges() {
        for( auto eit = _DT.finite_edges_begin(); eit != _DT.finite_edges_end(); ++eit ) {
            Delaunay_triangulation_2::Edge e = *eit;
            Vertex_handle p = e.first->vertex( (e.second+1)%3 ),
                          q = e.first->vertex( (e.second+2)%3 );
            add_edge( p, q );
        }
    }

    double get_angle( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) const {
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

    Vertex_circulator orient_circulator( const Vertex_circulator& C, const Vertex_handle& v ) const {
        Vertex_circulator out(C),
                          done(C);
        do {
            if( out->handle() == v )
                return out;
        } while( --out != done );
        return Vertex_circulator();
    }

    /**
     *  Given a Delaunay Triangulation DT and an output list out, compute the canonical ordering of
     *  the underlying point set.
     */
    template< typename OutputIterator >
    void canonical_order( OutputIterator ordering ) const {
        Timer t(",");
        VertexHash on_outer_face, complete;
        queue<Vertex_handle> ready;
        size_t i = n();

        vector<Vertex_handle> order(i);

        Vertex_circulator v_convex_hull = _DT.incident_vertices( _DT.infinite_vertex() ), // create a circulator of the convex hull
                          done( v_convex_hull );

        do { // cycle through convex hull to set vertex info
            on_outer_face.insert( v_convex_hull->handle() );
            ready.push( v_convex_hull->handle() );
            //_algoTV.addToEventQueue( v_convex_hull, 0 );
        } while( ++v_convex_hull != done );

        // Reserve v_1 and v_2 so we can guarantee they are on the convex hull
        for( size_t i=0; i<2; ++i ) {
            order[i] = ready.front();
            ready.pop();
            //_algoTV.addToEventQueue(order[i], 0 );
        }

        Vertex_handle v_k = ready.front();

        while( !ready.empty() ) {
            v_k = ready.front();
            //std::cout<<v_k->point()<<" ";
            ready.pop();
            //_algoTV.addToEventQueue( v_k, 0 );
            //_algoTV.addToEventQueue( v_k, false );

            if( incident_chords( on_outer_face, v_k ) > 0 ) {
                ready.push(v_k);
                //std::cout<<"requeued";
            } else {
                //std::cout<<"processed";
                order[--i] = v_k;
                on_outer_face.erase(v_k);
                // add all neighbors not complete or on outer face to ready list
                // add all neighbors not complete to on outer face
                Vertex_circulator N = _DT.incident_vertices(v_k),
                                  done(N);
                do {
                    if( !contains( complete, N->handle() ) && !_DT.is_infinite( N->handle() ) ) {
                        if( !contains( on_outer_face, N->handle() ) )
                            ready.push( N->handle() );
                        on_outer_face.insert( N->handle() );
                    }
                } while( ++N != done );
                complete.insert(v_k);
            }
            //std::cout<<"\n";
        }
        std::copy( order.begin(), order.end(), ordering );
    }

    size_t incident_chords( const VertexHash& on_outer_face, const Vertex_handle& v_k ) const {
        Vertex_circulator N = _DT.incident_vertices(v_k),
                          done(N);
        size_t c = 0;
        do {
            if( contains( on_outer_face, N->handle() ) )
                ++c;
        } while( ++N != done );
        // v_k is guaranteed to be incident to two vertices in on_outer_face, its neighbors.
        // Anything >2 is a chord
        assert( c >= 2 );
        return (c - 2);
    }

    int count_valid_neighbors( Vertex_circulator C, const VertexHash& invalid ) const {
        Vertex_circulator done(C);
        int k = 0; // count N_i path length k

        do {
            if( !contains( invalid, C ) && !_DT.is_infinite(C) )
                k++;
        } while( ++C != done );

        return k;
    }

    void normalize_circulator( Vertex_circulator &C, const VertexHash& invalid ) const {
        Vertex_circulator done = C;
        // Position circulator so that we are guaranteed to be on the first vertex on the path N_i
        // First, loop until the circulator reaches an invalid vertex or completes a full rotation
        while( !contains( invalid, C ) && !_DT.is_infinite(C) && ++C != done );// cout<<v_n->point()<<"\n";
        done = C;
        // Loop until the circulator reaches a valid vertex
        while( ( contains( invalid, C ) || _DT.is_infinite(C) ) && ++C != done );// cout<<v_n->point()<<"\n";
    }

    size_t n() const {
        return _DT.number_of_vertices();
    }

  private:

    void add_half_edge( const Vertex_handle& v1, const Vertex_handle& v2 ) {
        typename AdjacencyList::iterator v1_incident_vertices = _E.find(v1);

        if( v1_incident_vertices == _E.end() ) // v1 not found in g
            tie( v1_incident_vertices, ignore ) = _E.insert( make_pair( v1, VertexSet() ) );

        v1_incident_vertices->second.insert(v2);
    }

    void remove_half_edge( const Vertex_handle& v1, const Vertex_handle& v2 ) {
        typename AdjacencyList::iterator v1_incident_vertices = _E.find(v1);

        // if v1 in present in _E, and v2 is present in v1, remove v2 from v1
        if( contains( _E, v1 ) && contains( v1_incident_vertices->second, v2 ) )
            v1_incident_vertices->second.erase(v2);
    }

}; // class DelaunayGraph


// Make some key types of DelaunayGraph public in the namespace

template<class V >
using VertexMap = DelaunayGraph::template VertexMap<V>;
using VertexSet = DelaunayGraph::VertexSet;
using VertexHash = DelaunayGraph::VertexHash;
using Vertex_handle = DelaunayGraph::Vertex_handle;
using Vertex_circulator = DelaunayGraph::Vertex_circulator;
using Point = DelaunayGraph::K::Point_2;


} // namespace gsnunf

#endif // GSNUNF_GRAPH_H
