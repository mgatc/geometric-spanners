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
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Vector_2.h>
#include <CGAL/utils.h> // min, max

#include "utilities.h"

namespace gsnunf {

using namespace std;




typedef CGAL::Triangulation_vertex_base_with_info_2<index_t, K> Vb;
typedef CGAL::Triangulation_face_base_2<K>                      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>            Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                  Delaunay_triangulation;

typedef Delaunay_triangulation::Vertex_handle                   Vertex_handle;
typedef Delaunay_triangulation::Vertex_circulator               Vertex_circulator;
typedef Delaunay_triangulation::Finite_vertices_iterator        Finite_vertices_iterator;
typedef Delaunay_triangulation::Finite_edges_iterator           Finite_edges_iterator;
typedef Delaunay_triangulation::Face_handle                     Face_handle;

typedef set<Vertex_handle> VertexSet;
typedef unordered_set<Vertex_handle> VertexHash;
template< typename V >
using VertexMap = unordered_map< Vertex_handle, V >;
using AdjacencyList = VertexMap< VertexSet >;
//template< typename N >
//using EdgeInfoMap = VertexMap< VertexMap< optional<N> > >;

//typedef CGAL::Container_from_circulator<Vertex_circulator> Vertex_container;

class DelaunayGraph {
  public:

    /* Data */
    Delaunay_triangulation m_DT;
    AdjacencyList m_E;


    /* Functions */
    DelaunayGraph() = default;
    template< typename RandomAccessIterator >
    DelaunayGraph( RandomAccessIterator pointsBegin,
                   RandomAccessIterator pointsEnd )
    {
        vector<Point> P(pointsBegin, pointsEnd);
        vector<size_t> index;
        spatialSort<K>(P, index);

        /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
          (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
        Face_handle hint;
        for(size_t entry : index) {
            auto vh = m_DT.insert(P[entry], hint);
            hint = vh->face();
            vh->info() = entry;
        }
    }

    template< typename RandomAccessIterator >
    // kitty: bnnnnn23333333333333333333333,m000000000000000000000000000000000000000000000000000000000re	IIIKKKKKKKKKKKHUG0
    void buildFromEdgeList( RandomAccessIterator edgesBegin, RandomAccessIterator edgesEnd ) {
        vector<pair<Point,Point>> edges( edgesBegin, edgesEnd );
        std::sort( edges.begin(), edges.end() );

        for( auto e=edgesBegin; e!=edgesEnd; ++e ) {
            Vertex_handle u = m_DT.insert(e->first);
            Vertex_handle v = m_DT.insert(e->second);
            add_edge(u,v);
        }
    }

    inline void add_edge( const Vertex_handle& v1, const Vertex_handle& v2 ) {
        add_half_edge( v1, v2 );
        add_half_edge( v2, v1 );
        //cout<<"    add edge ("<<v1->point()<<", "<<v2->point()<<")\n";
        //_algoTV.addToEventQueue( std::make_pair(v1,v2), true );
    }

    inline void remove_edge( const Vertex_handle& v1, const Vertex_handle& v2 ) {
        remove_half_edge( v1, v2 );
        remove_half_edge( v2, v1 );
        //_algoTV.addToEventQueue( std::make_pair(v1,v2), false );
    }

    void add_all_edges() {
        for( auto eit = m_DT.finite_edges_begin(); eit != m_DT.finite_edges_end(); ++eit ) {
            Delaunay_triangulation::Edge e = *eit;
            Vertex_handle p = e.first->vertex( (e.second+1)%3 ),
                          q = e.first->vertex( (e.second+2)%3 );
            add_edge( p, q );
        }
    }
//
//    number_t get_angle( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) const {
//        assert( !m_DT.is_infinite(p) && !m_DT.is_infinite(q) && !m_DT.is_infinite(r) );
//
//        Vector_2 pq( p->point(), q->point() );
//        Vector_2 rq( r->point(), q->point() );
//
//        number_t result = atan2( rq.y(), rq.x()) - atan2(pq.y(), pq.x() );
//
//        // atan() returns a value between -PI and PI. From zero ("up"), CCW rotation is negative and CW is positive.
//        // Our zero is also "up," but we only want positive values between 0 and 2*PI:
//
//        result *= -1; // First, invert the result. This will associate CW rotation with positive values.
//        if( result < EPSILON ) // Then, if the result is less than 0 (or epsilon for floats) add 2*PI.
//            result += 2*PI;
//
//        return result;
//    }

    inline Vertex_circulator orient_circulator( const Vertex_circulator& C, const Vertex_handle& v ) const {
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
        //Timer t(",");
        VertexHash on_outer_face, complete;
        queue<Vertex_handle> ready;
        size_t i = n();

        vector<Vertex_handle> order(i);

        Vertex_circulator v_convex_hull = m_DT.incident_vertices( m_DT.infinite_vertex() ), // create a circulator of the convex hull
                          done( v_convex_hull );

        do { // cycle through convex hull to set vertex info
            on_outer_face.insert( v_convex_hull->handle() );
            ready.push( v_convex_hull->handle() );
            //_algoTV.addToEventQueue( v_convex_hull, 0 );
        } while( ++v_convex_hull != done );

        // Reserve v_1 and v_2 so we can guarantee they are on the convex hull
        for(size_t j=0; j < 2; ++j ) {
            order[j] = ready.front();
            ready.pop();
            //_algoTV.addToEventQueue(order[j], 0 );
        }

        Vertex_handle v_k = ready.front();

        while( !ready.empty() ) {
            v_k = ready.front();
            //std::cout<<v_k->point()<<" ";
            ready.pop();
//            _algoTV.addToEventQueue( v_k, 0 );
//            _algoTV.addToEventQueue( v_k, false );

            if( incident_chords( on_outer_face, v_k ) > 0 ) {
                ready.push(v_k);
                //std::cout<<"requeued";
            } else {
                //std::cout<<"processed";
                order[--i] = v_k;
                on_outer_face.erase(v_k);
                // add all neighbors not complete or on outer face to ready list
                // add all neighbors not complete to on outer face
                Vertex_circulator N = m_DT.incident_vertices(v_k);
                done = N;
                do {
                    if( !contains( complete, N->handle() ) && !m_DT.is_infinite( N->handle() ) ) {
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

    inline size_t incident_chords( const VertexHash& on_outer_face, const Vertex_handle& v_k ) const {
        Vertex_circulator N = m_DT.incident_vertices(v_k),
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

    inline int count_valid_neighbors( Vertex_circulator C, const VertexHash& invalid ) const {
        Vertex_circulator done(C);
        int k = 0; // count N_i path length k

        do {
            if( !contains( invalid, C ) && !m_DT.is_infinite(C) )
                k++;
        } while( ++C != done );

        return k;
    }

    inline void normalize_circulator( Vertex_circulator &C, const VertexHash& invalid ) const {
        Vertex_circulator done = C;
        // Position circulator so that we are guaranteed to be on the first vertex on the path N_i
        // First, loop until the circulator reaches an invalid vertex or completes a full rotation
        while( !contains( invalid, C ) && !m_DT.is_infinite(C) && ++C != done );// cout<<v_n->point()<<"\n";
        done = C;
        // Loop until the circulator reaches a valid vertex
        while( ( contains( invalid, C ) || m_DT.is_infinite(C) ) && ++C != done );// cout<<v_n->point()<<"\n";
    }

    inline size_t n() const {
        return m_DT.number_of_vertices();
    }

    inline size_t degree() const {
        size_t max_d = 0;
        for( auto& v : m_E )
            max_d = std::max( v.second.size(), max_d );
        return max_d;
    }

  private:

    inline void add_half_edge( const Vertex_handle& v1, const Vertex_handle& v2 ) {
        auto v1_incident_vertices = m_E.find(v1);

        if( v1_incident_vertices == m_E.end() ) // v1 not found in g
            tie( v1_incident_vertices, ignore ) = m_E.insert( make_pair( v1, VertexSet() ) );

        v1_incident_vertices->second.insert(v2);
    }

    inline void remove_half_edge( const Vertex_handle& v1, const Vertex_handle& v2 ) {
        auto v1_incident_vertices = m_E.find(v1);

        // if v1 in present in _E, and v2 is present in v1, remove v2 from v1
        if( contains( m_E, v1 ) && contains( v1_incident_vertices->second, v2 ) )
            v1_incident_vertices->second.erase(v2);
    }

}; // class DelaunayGraph

} // namespace gsnunf

#endif // GSNUNF_DELAUNAYGRAPH_H
