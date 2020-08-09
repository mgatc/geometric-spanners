#ifndef GSNUNF_TRANSFORMPOLYGON_H
#define GSNUNF_TRANSFORMPOLYGON_H

#include <algorithm>
#include <experimental/memory_resource>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DelaunayGraph.h"

namespace gsnunf {

/* Represents a single split of a split vertex where
 * first = the vertex to be split
 * second = the split vertex's s_1
 */
template< class T >
struct SplitVertex {
    using v_type = Vertex_handle<T>;
    using s_1_type = SplitVertex<T>*;
    v_type v;
    s_1_type s_1;
    SplitVertex() {}
    SplitVertex( const v_type v, const s_1_type s_1 ) : v(v), s_1(s_1) {}
    SplitVertex( const SplitVertex<T>& other ) : v(other.v), s_1(other.s_1) {}

    SplitVertex& operator=( const SplitVertex& other ) { //copy assignment
        if( this != &other ) {
            v = other.v;
            s_1 = other.s_1;
        }
        return *this;
    }
};

struct SplitVertexHasher {
    template< class T >
    std::size_t operator()( const SplitVertex<T>& k ) const {
        // Compute individual hash values for first, second and third
            // http://stackoverflow.com/a/1646913/126995
        size_t res = 17;
        res = res * 31 + hash< Vertex_handle<T> >()( k.v );
        res = res * 31 + hash< Vertex_handle<T> >()( k.s_1->v );
        return res;
    }
    template< class T >
    std::size_t operator()( const SplitVertex<T>* k ) const {
        // Compute individual hash values for first, second and third
            // http://stackoverflow.com/a/1646913/126995
        size_t res = 17;
        res = res * 31 + hash< Vertex_handle<T> >()( k->v );
        res = res * 31 + hash< Vertex_handle<T> >()( k->s_1->v );
        return res;
    }
};

template< class T >
struct CompareSplitVertex {
    bool operator()( const SplitVertex<T>& a, const SplitVertex<T>& b ) {
        return a.v < b.v || ( a.v == b.v && a->s_1 < b->s_1 );
    }
};
//using SplitVertex = std::pair< Vertex_handle<T>, SplitVertex<T>& >;
// Compare two split vertices for equality
template< class T >
inline bool operator==( const SplitVertex<T>& lhs,
                        const SplitVertex<T>& rhs ) {
    return lhs.v == rhs.v && lhs.s_1->v == rhs.s_1->v;
}

// Compare two split vertices for inequality
template< class T >
inline bool operator!=( const SplitVertex<T>& lhs,
                        const SplitVertex<T>& rhs ) {
    return !(lhs == rhs);
}

template< class T >
using SplitVertexMap = VertexMap<T, VertexMap<T, SplitVertex<T> > >;

template< class T >
using IncidentSplitVertexContainer = std::unordered_set< SplitVertex<T>, SplitVertexHasher >;

template< class T >
using UnsplitVertex = VertexMap<T, IncidentSplitVertexContainer<T> >;

template< class T >
using SplitVertexEdgeMap = std::unordered_map< SplitVertex<T>*, IncidentSplitVertexContainer<T>, SplitVertexHasher >;
//    std::experimental::fundamentals_v2::pmr::polymorphic_allocator< SplitVertex<T> > > > >;

template< class T >
void //SplitVertex<T>
add_half_edge( SplitVertexEdgeMap<T>& E, SplitVertex<T>* a, SplitVertex<T>* b ) {
    E[a].emplace(*b);
}

template< class T >
inline SplitVertex<T>*
add_vertex( SplitVertexMap<T>& V, const Vertex_handle<T>& v, SplitVertex<T>* s_1 ) {
    //cout<<"add "<<a.v->point()<<" ("<< a.s_1->v->point()<<") - "<<b.v->point()<<" ("<<b.s_1->v->point()<<")"<<endl;
    V[v][s_1->v] = SplitVertex<T>( v, s_1 );

    return &V.at(v).at(s_1->v);
}

template< class T >
void//inline pair< SplitVertex<T>, SplitVertex<T> >
add_edge( SplitVertexEdgeMap<T>& E, SplitVertex<T>* a, SplitVertex<T>* b ) {
    cout<<"edge add "<<a->v->point()<<" ("<< a->s_1->v->point()<<") - "<<b->v->point()<<" ("<<b->s_1->v->point()<<")"<<endl;

   // return {
        add_half_edge( E, b, a );
        add_half_edge( E, a, b );
   // };
}

//template< class T >
//inline pair< SplitVertex<T>, SplitVertex<T> >
//add_edge( const SplitVertexEdgeMap<T>& V, SplitVertexEdgeMap<T>& E,
//                const Vertex_handle<T>& a, const Vertex_handle<T>& s_1_a,
//                const Vertex_handle<T>& b, const Vertex_handle<T>& s_1_b ) {
//    //cout<<"add "<<a.v->point()<<" ("<< a.s_1->v->point()<<") - "<<b.v->point()<<" ("<<b.s_1->v->point()<<")"<<endl;
//    // TODO///////////////////////////////////////////////////////////////
//    auto& list_a = V[a][s_1_a];
//    *list_a.emplace( b, nullptr ).first;
//
//    auto& list_b = E[b][s_1_b];
//    list_b.emplace( a, &list_a.back() );
//    list_a.back().s_1 = &list_b.back();
//
//    return {
//        &list_a.back(),
//        &list_b.back()
//    };
//}

// Print the contents of a split vertex edge list
template< class T >
void print_edges( const SplitVertexEdgeMap<T>& E ) {
    for( auto v1 : E ) {
        cout<< "u: "<<v1.first->v->point()<<" s1: "<<v1.first->s_1->v->point()<<"\n";
        for( auto v2 : v1.second )
            cout<<"  v: "<<v2->v->point()<<" s1: "<<v2->s_1->v->point()<<"\n";
        cout<<"\n";

    }
}

// Print the contents of a split vertex list
template< class T >
void print_vertices( SplitVertexMap<T>& V ) {
    for( auto v : V ) {
        cout<< "u: "<<v.first->point()<<"\n";
        for( auto s1 : v.second )
            cout<<"  s1: "<<s1.first->point()<<"\n";
        cout<<"\n";

    }
}

template< class T >
void TransformPolygon( const DelaunayGraph<T>& SG, SplitVertexMap<T>& V, SplitVertexEdgeMap<T>& P ) {
    P.clear();
    V.clear();
    Vertex_handle<T> v_inf = SG._DT.infinite_vertex();

    // Convex hull circulator
    Vertex_circulator<T> s_1_finder = SG._DT.incident_vertices( v_inf ),
                         v_1_finder = SG._DT.incident_vertices( s_1_finder );

    v_1_finder = SG.orient_circulator( v_1_finder, v_inf );
    ++v_1_finder; // rotate once CCW

    SplitVertex<T> s_first( s_1_finder, nullptr ), // v_1 will need an s_1 to point to
                       v_1( v_1_finder, &s_first ),
                    v_next( v_1 ); // used in the traversal to hold vertex before insertion

    SplitVertex<T>* s_1 = &s_first;
    SplitVertex<T>* v_i = add_vertex<T>( V, v_1_finder, s_1 );

    Vertex_circulator<T> N;

    do {
        N = SG._DT.incident_vertices( v_i->v );
        while( (--N)->handle() != v_i->s_1->v ); // rotate N until it points to s_1

        const VertexSet<T>& N_SG = SG._E.find( v_i->v )->second; // get neighbors of v_i in SG._E
        while( !contains( N_SG, --N ) ); // rotate CW until reaching a neighbor in SG
        s_1 = v_i;
        v_next = SplitVertex<T>( N, s_1 );
        if( v_next == v_1 ) V.at( v_1.v ).erase( v_1.s_1->v ); // v_1's key points to a temporary
        v_i = add_vertex( V, N, s_1 );//[v_next.v][v_next.s_1->v] = v_next;

        assert( !SG._DT.is_infinite(N) && !SG._DT.is_infinite(s_1->v) );
        add_edge<T>( P, v_i, s_1 );

    } while( *v_i != v_1 );


//    cout<<endl;
//    print_edges<T>(P);
//    cout<<endl;
}


} // namespace gsnunf

#endif // GSNUNF_TRANSFORMPOLYGON_H

