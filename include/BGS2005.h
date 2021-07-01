#ifndef GSNUNF_BGS2005_H
#define GSNUNF_BGS2005_H

#include <algorithm> // min, swap
#include <functional> // hash
#include <iostream>
#include <list>
#include <memory>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "DelaunayGraph.h"
#include "GeometricSpannerPrinter.h"
#include "metrics.h"

namespace gsnunf {

namespace bgs2005 {

struct SplitVertexSet;
struct SplitVertexHasher;

/* Represents a single split of a split vertex where
 * first = the vertex to be split
 * second = the split vertex's s_1
 */
struct SplitVertex {
    using v_type = Vertex_handle;
    using key_type = size_t;

    SplitVertexSet* V;
    v_type v;
    key_type s_1;
    v_type s_1_handle;
    key_type key;

    SplitVertex() {}
    SplitVertex( const v_type& v )
        : v(v) {}
    SplitVertex( const v_type& v, const SplitVertex& s_1 )
        : v(v), s_1( s_1.key ), s_1_handle(s_1.v) {}
    SplitVertex( const v_type& v, const key_type& s_1 )
        : v(v), s_1( s_1 ) {}
    SplitVertex( const SplitVertex& other )
        : V(other.V), v(other.v), s_1(other.s_1), s_1_handle( other.s_1_handle ), key(other.key) {}

    SplitVertex& operator=( const SplitVertex& other ) { //copy assignment
        if( this != &other ) {
            V = other.V;
            v = other.v;
            key = other.key;
            s_1 = other.s_1;
            s_1_handle = other.s_1_handle;
        }
        return *this;
    }
};

ostream& operator<<( ostream& os, const SplitVertex& v ){
    os << v.v->point() << " (" << v.s_1_handle->point() << ")";
    return os;
}

struct CompareSplitVertex {
    bool operator()( const SplitVertex& a, const SplitVertex& b ) {
        return a.v < b.v || ( a.v == b.v && a.s_1 < b.s_1 );
    }
};

// Compare two split vertices for equality
inline bool operator==( const SplitVertex& lhs,
                        const SplitVertex& rhs ) {
    return lhs.v == rhs.v && lhs.s_1_handle == rhs.s_1_handle;
}

// Compare two split vertices for inequality
inline bool operator!=( const SplitVertex& lhs,
                        const SplitVertex& rhs ) {
    return !(lhs == rhs);
}

using key_type = typename SplitVertex::key_type;
using IncidentSplitVertexContainer = std::set< key_type >;
using SplitVertexEdgeMap = std::unordered_map< key_type, IncidentSplitVertexContainer >;

struct SplitVertexSet {
    VertexMap< VertexMap< size_t > > index;

    vector< SplitVertex > V;

    size_t insert( const SplitVertex& v ) {
        size_t vector_key;
        SplitVertex v_insert = v;
        if( contains( index, v_insert.v ) && contains( index.at(v_insert.v), v_insert.s_1_handle ) ) {
            vector_key = index[v_insert.v][v_insert.s_1_handle];
            v_insert.key = vector_key;
            V[vector_key] = v_insert;
        } else {
            vector_key = insert_in_container(v_insert);
            insert_in_map(v_insert);
        }
        return vector_key;
    }
    size_t insert_in_container( SplitVertex& v ) {
        size_t vector_key = V.size();
        v.key = vector_key;
        v.V = this;
        V.push_back(v);
        return vector_key;
    }
    void insert_in_map( const SplitVertex& v ) {
        index[v.v].insert_or_assign( v.s_1_handle, v.key );
    }
    SplitVertex at( const size_t i ) const {
        return V.at(i);
    }
    size_t at( const SplitVertex& v ) const {
        return index.at(v.v).at( V.at(v.s_1).v );
    }
    SplitVertex& at( const Vertex_handle v, const Vertex_handle s_1 ) {
        return V.at( index.at(v).at(s_1) );
    }
};

inline void add_half_edge( SplitVertexEdgeMap& E, const SplitVertex& a, const SplitVertex& b ) {
    E[a.key].emplace(b.key);
}

inline size_t add_vertex( SplitVertexSet& V, const Vertex_handle& v, const SplitVertex s_1 ) {
    return V.insert( SplitVertex( v, s_1 ) );
}

inline void add_edge( const DelaunayGraph& SG, SplitVertexEdgeMap& E, const SplitVertex& a, const SplitVertex& b ) {
    assert( SG._DT.is_edge( a.v, b.v ) );

    add_half_edge( E, b, a );
    add_half_edge( E, a, b );
}

// Print the contents of a split vertex edge list
void print( const SplitVertexSet& V, const SplitVertexEdgeMap& E ) {
    for( auto v1 : E ) {
        key_type k = v1.first;
        SplitVertex u = V.at(k);
        cout<< "u: "<<u.v->point()<<" s1: "<<V.at(u.s_1).v->point()<<"\n";
        for( auto v2 : v1.second )
            cout<<"  v: "<<V.at(v2).v->point()<<" s1: "<<V.at(V.at(v2).s_1).v->point()<<"\n";
        cout<<"\n";

    }
}

// Print the contents of a split vertex list
void print( const SplitVertexSet& V ) {
    for( auto v : V.V ) {
        cout<<"v"<<v.key<<": "<<v;
        cout<<" s1_key: "<<v.s_1<<"\n";
    }
}

struct SplitVertexHasher {
    std::size_t operator()( const SplitVertex& k ) const {
        // Compute individual hash values for first, second and third
            // http://stackoverflow.com/a/1646913/126995
        size_t res = 17;
         res = res * 31 + hash< Vertex_handle >()( k.v );
         res = res * 31 + hash< Vertex_handle >()( k.V->at(k.s_1).v );
        return res;
    }
};



namespace spanning_graph {

inline void add_first_edge( DelaunayGraph& G, Vertex_handle v, Vertex_circulator C ) {
    Vertex_handle v2 = C->handle();
    G.add_edge( v, v2 );

    //G.addToEventQueue( v, 1 );
    //G.addToEventQueue( v2, 2 );
    //G.addToEventQueue( { v, v2 }, true );
}

inline void add_second_edge( DelaunayGraph& G, Vertex_handle v, Vertex_circulator C ) {
    while( G._DT.is_infinite(++C) );
    Vertex_handle v2 = C->handle();
    G.add_edge( v, v2 );

//        G.addToEventQueue( v, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v, v2 }, true );
}

inline void add_last_edge( DelaunayGraph& G, const Vertex_handle& v, Vertex_circulator C, const VertexHash& is_removed ) {
    --C;
    Vertex_circulator done(C);

    while( ( contains( is_removed, C ) || G._DT.is_infinite(C) ) && --C != done );

    Vertex_handle v2 = C->handle();

    G.add_edge( v, v2 );

//        G.addToEventQueue( v, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v, v2 }, true );
}

inline void remove_first_edge( DelaunayGraph& G, Vertex_circulator C ) {
    Vertex_handle v1 = C->handle(),
                  v2 = (++C)->handle();
    G.remove_edge( v1, v2 );

//        G.addToEventQueue( v1, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v1, v2 }, false );
}

inline void remove_second_edge( DelaunayGraph& G, Vertex_circulator C ) {
    Vertex_handle v1 = (++C)->handle(),
                  v2 = (++C)->handle();
    G.remove_edge( v1, v2 );

//        G.addToEventQueue( v1, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v1, v2 }, false );
}

inline void remove_last_edge( DelaunayGraph& G, Vertex_circulator C, const VertexHash& is_removed ) {
    --C;

    Vertex_circulator done(C);

    while( ( contains( is_removed, C ) || G._DT.is_infinite(C) ) && --C != done );

    Vertex_handle v1 = C->handle(),
                     v2 = (--C)->handle();
    G.remove_edge( v1, v2 );

//        G.addToEventQueue( v1, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v1, v2 }, false );
}

}; // namespace spanning_graph

void SpanningGraph( DelaunayGraph& G ) {
    using namespace spanning_graph;

    vector< Vertex_handle > canonical;

    G.canonical_order( inserter( canonical, canonical.end() ) );
    //Timer timer(",");

    Vertex_circulator v_n, done;
    size_t i;

    VertexHash is_removed( canonical.begin(), canonical.end() );

    // Add first three vertices from canonical
    for( i=0; i<3; ++i ) { // Add edges of triangle
        is_removed.erase( canonical.at(i) );
        G.add_edge( canonical.at(i), canonical.at((i+1)%3) );
        //G.addToEventQueue( { canonical.at(i), canonical.at((i+1)%3) }, true );
        //G.addToEventQueue( { canonical.at(i), canonical.at((i+1)%3) }, true );
    }
    // Add the rest of the vertices from canonical
    for( i=i; i<canonical.size(); ++i ) {
        is_removed.erase( canonical.at(i) );
//        G.addToEventQueue( *c_iter, 0 );        // activate c_iter
//        G.addToEventQueue( *c_iter, true );

        v_n = G._DT.incident_vertices( canonical.at(i) );
        done = v_n;

        G.normalize_circulator( v_n, is_removed );
        done = v_n;

        int k = G.count_valid_neighbors( v_n, is_removed );

        if( k == 2 ) {
            // remove edge between first two vertices
            remove_first_edge( G, v_n );
            // add edge between canonical iterator and first vertex
            add_first_edge( G, canonical.at(i), v_n );
            // add edge between canonical iterator and second vertex
            add_second_edge( G, canonical.at(i), v_n );

        } else if( k > 2 ) {
            // remove edge between first two vertices
            remove_first_edge( G, v_n );
            // remove edge between last two vertices
            remove_last_edge( G, v_n, is_removed );
            // add edge between canonical iterator and first vertex
            add_first_edge( G, canonical.at(i), v_n );
            // add edge between canonical iterator and second vertex
            add_second_edge( G, canonical.at(i), v_n );
            // add edge between canonical iterator and last vertex
            add_last_edge( G, canonical.at(i), v_n, is_removed );
        }
    }

    // Test assumption
    for( auto it=G._E.begin(); it!=G._E.end(); ++it ) { // for all v_i, 1<=i<=n
        assert( it->second.size() <= 3 );                 // |v_i| <= 3
    }
}



namespace transform_polygon {

inline Vertex_handle find_s_1_handle( const DelaunayGraph& SG, const pair< const Vertex_handle, VertexMap< size_t > >& unsplit, const Vertex_handle& v_n ) {
    Vertex_handle v_i = unsplit.first;
    Vertex_circulator N = SG._DT.incident_vertices(v_i); // get circulator around unsplit.first
    while( (++N)->handle() != v_n );
    while( !contains( unsplit.second, N->handle() ) ) ++N; // orient to a neighbor in unsplit.second

    return N;
}

} // namespace transform_polygon

void TransformPolygon( const DelaunayGraph& SG, SplitVertexSet& V, SplitVertexEdgeMap& P ) {
    using namespace transform_polygon;
    P.clear();
    Vertex_handle v_inf = SG._DT.infinite_vertex();

    // Convex hull circulator
    Vertex_circulator s_1_finder = SG._DT.incident_vertices( v_inf ),
                      v_1_finder = SG._DT.incident_vertices( s_1_finder );

    v_1_finder = SG.orient_circulator( v_1_finder, v_inf );
    ++v_1_finder; // rotate once CCW

    SplitVertex s_1( s_1_finder, 0 ); // v_1 will need an s_1 to point to

    SplitVertex v_1 = V.at( add_vertex( V, v_1_finder, s_1 ) ),
                v_i = v_1,
                v_next = v_1;
    //cout<<v_i<<"\n";

    Vertex_circulator N;

    do {
        N = SG._DT.incident_vertices( v_i.v );
        while( (--N)->handle() != s_1.v ); // rotate N until it points to s_1

        const VertexSet& N_SG = SG._E.find( v_i.v )->second; // get neighbors of v_i in SG._E
        while( !contains( N_SG, --N ) ); // rotate CW until reaching a neighbor in SG
        s_1 = v_i;
        v_next = SplitVertex( N, s_1 );

        if( v_next == v_1 ) {       // if we've looped all the way around
            SplitVertex& v_ref = V.V.at(v_1.key);
            v_ref.s_1 = s_1.key; // update v_1's s_1
            v_i = v_ref; // update v_i to explicitly point to v_1
        } else { // add new vertex
            v_i = V.at( add_vertex( V, N, s_1 ) );
        }

        assert( !SG._DT.is_infinite(N) );
        add_edge( SG, P, v_i, s_1 );

    } while( v_i != v_1 );

    for( auto e=SG._DT.finite_edges_begin(); e!=SG._DT.finite_edges_end(); ++e ) {
        Vertex_handle u = e->first->vertex( (e->second+1)%3 ),
                         v = e->first->vertex( (e->second+2)%3 );
        if( !contains( SG._E.at(u), v ) ) { // only add edges that are not in SG
            pair< const Vertex_handle, VertexMap< key_type > >&
                unsplit_u = *V.index.find( u ),
                unsplit_v = *V.index.find( v );
            Vertex_handle s_1_u = find_s_1_handle( SG, unsplit_u, v ),
                             s_1_v = find_s_1_handle( SG, unsplit_v, u );
            add_edge( SG, P, V.at( u, s_1_u ), V.at( v, s_1_v ) );
        }
    }
}



namespace polygon_spanner {

enum VertexStatus { unknown, known, complete };
using VertexStatusMap = unordered_map< key_type, VertexStatus >;

/*
 * Iterates from s_1 to s_m (inclusive) and performs "foreach" on each vertex.
 * Provides a SplitVertex of the currently pointed-to vertex.
 */
template< typename F >
void foreach_neighbor( const DelaunayGraph& SG, const SplitVertexSet& V, const SplitVertexEdgeMap& E, const SplitVertex& v_i, F foreach ) {
    assert( !SG._DT.is_infinite( v_i.v ) );
    /*
     *    The circulator provided by SG._DT will provide the correct ordering.
     * However, the vertices provided by the circulator are not split. Aside
     * from s_1 and s_m, which can have equal primary vertices (with differing
     * s_1s of course), the primary vertices of these neighbors will be unique.
     * Therefore, we will create an unordered set of incident vertices to v_i
     * and provide a custom hash function that will simply pass along the
     * CGAL Vertex_handle hash. This way, we can match a split vertex from the
     * vertex handle provided by the circulator in constant time.
     */

    struct HashOnVertexHandle {
        size_t operator()( const SplitVertex& k ) const {
            return hash< Vertex_handle >()( k.v );
        }
    };
    struct CompareOnVertexHandle {
        bool operator()( const SplitVertex& a, const SplitVertex& b ) const {
            return a.v == b.v;
        }
    };
    unordered_set< SplitVertex, HashOnVertexHandle, CompareOnVertexHandle > N_E;

    for( auto k: E.at( v_i.key ) ) {
        if( k != v_i.s_1 ) { // don't enter s_1
            N_E.insert( V.at(k) );
        }
    }
    //cout<<"N_E size:"<<N_E.size()<<"\n";

    Vertex_circulator N = SG._DT.incident_vertices( v_i.v );
    N = SG.orient_circulator( N, v_i.s_1_handle );


    SplitVertex v_neighbor = V.at( v_i.s_1 );
    bool include_s_1 = contains( E.at(v_i.key), v_neighbor.key );
    //cout<< "    start foreach: "<<v_i<<" n: "<<(N_E.size()+size_t(include_s_1))<<"\n";

    if( include_s_1 ) {
    //    cout<<"    v_neighbor: "<<v_neighbor<<"\n";
        foreach( v_neighbor );
    }

    do {
        --N;
        v_neighbor = SplitVertex( N->handle() );
        if( contains( N_E, v_neighbor ) ) {
            v_neighbor = *N_E.find( v_neighbor );
            //cout<<"    v_neighbor: "<<v_neighbor<<"\n";
            foreach( v_neighbor );
            N_E.erase( v_neighbor );
        }
    } while( !N_E.empty() );
}

inline SplitVertex get_s_m( const DelaunayGraph& SG, const SplitVertexSet& V, const SplitVertexEdgeMap& E, const SplitVertex& v_split ) {
    key_type s_i;
    foreach_neighbor( SG, V, E, v_split, [&]( SplitVertex& v_n ) {
        s_i = v_n.key;
    });
    return V.at(s_i);
}

/*
 * Adds a split vertex edge, asserting that neither vertices a nor b are marked as complete
 */
inline void add_polygon_spanner_edge( const DelaunayGraph& SG, SplitVertexEdgeMap& E,
                              SplitVertex& a, SplitVertex& b,
                              const VertexStatusMap& status ) {
    assert( !contains( status, a.key ) || status.at(a.key) != complete );
    assert( !contains( status, b.key ) || status.at(b.key) != complete );
    return add_edge( SG, E, a, b );
}

void add_cross_edges( const DelaunayGraph& SG, SplitVertexSet& V, SplitVertexEdgeMap& E, SplitVertexEdgeMap& E_P,
                      SplitVertex& p, SplitVertex& q, SplitVertex& r,
                      const VertexStatusMap& status ) {
    optional< SplitVertex > v_last = nullopt;
    bool isInZone = false;
    foreach_neighbor( SG, V, E, q, [&]( SplitVertex& v_n ) {
        if( isInZone && v_last ) {
            add_polygon_spanner_edge( SG, E_P, *v_last, v_n, status );
            //SG.addToEventQueue( {s_last.first, N}, true );
        }
        if( v_n.key == p.key ) isInZone = true;
        if( v_n.key == r.key ) isInZone = false;

        v_last = { v_n };
    });
}

void add_forward_edges( const DelaunayGraph& SG, SplitVertexSet& V, SplitVertexEdgeMap& E, SplitVertexEdgeMap& E_P,
                        SplitVertex& p, SplitVertex& q, SplitVertex& r,
                        const VertexStatusMap& status ) {

    //double deg = 180/PI; // used for displaying angles in degree

    double alpha = SG.get_angle( p.v, q.v, r.v );
    short subangles = rint( ceil( alpha / (PI/2) ) );
    double beta = alpha / subangles;

    vector< SplitVertex > add( subangles, SplitVertex( SG._DT.infinite_vertex() ) ); // initialize add to infinite vertex

    double theta;
    short i;
    bool isInZone = false;

    foreach_neighbor( SG, V, E, q, [&]( SplitVertex& v_n ) {
        // r is guaranteed to be already added or will be added immediately after this step
        // so set isInZone to false before processing it
        if( v_n.key == r.key ) isInZone = false;
        //SG.addToEventQueue( N, 1 );// focus1 on N
        if( isInZone ) {
            theta = SG.get_angle( p.v, q.v, v_n.v );
            if( theta > 2*PI-EPSILON )
                theta = 0;
            i = std::min( int(theta/beta), subangles-1 );

            if( SG._DT.is_infinite( add.at(i).v )
              || Vector_2( v_n.v->point(), q.v->point() ).squared_length() < Vector_2( add.at(i).v->point(), q.v->point() ).squared_length() )
                add.at(i) = v_n;   // if the saved vertex is infinite or longer than the current one, update
        }
        // p is guaranteed to be already added or will be added immediately after this step
        // so don't set isInZone to true until after it's passed
        if( v_n.key == p.key ) isInZone = true;
    } );

    for( SplitVertex v : add )
        if( !SG._DT.is_infinite( v.v ) ) {
            add_polygon_spanner_edge( SG, E_P, q, v, status );
            //SG.addToEventQueue( {q.first, v.first}, true ); // add edge
        }
}

inline void add_polygon_edges( const DelaunayGraph& SG, SplitVertexEdgeMap& E_P,
                        SplitVertex& p, SplitVertex& q, SplitVertex& r,
                        const VertexStatusMap& status ) {
    /*
     * Only add edges from P if they have not already been added.
     * Failure to check for this will result in breaking the
     * assumption that edges are only added to non-complete
     * vertices.
     */
    if( !contains( E_P, p.key) || !contains( E_P.at(p.key), q.key ) )
        add_polygon_spanner_edge( SG, E_P, p, q, status );
    if( !contains( E_P, q.key) || !contains( E_P.at(q.key), r.key ) )
        add_polygon_spanner_edge( SG, E_P, q, r, status );
}

inline void add_polygon_spanner_edges( const DelaunayGraph& SG, SplitVertexSet& V, SplitVertexEdgeMap& E, SplitVertexEdgeMap& E_P,
                                SplitVertex& p, SplitVertex& q, SplitVertex& r,
                                const VertexStatusMap& status ) {
    assert( !SG._DT.is_infinite(p.v) && !SG._DT.is_infinite(q.v) && !SG._DT.is_infinite(r.v) );
    if( p == r ) return;

    //cout<< " adding forward edges...\n";
    add_forward_edges( SG, V, E, E_P, p, q, r, status );
    //cout<< " adding cross edges...\n";
    add_cross_edges( SG, V, E, E_P, p, q, r, status );
}

void process_vertex( const DelaunayGraph& SG, SplitVertexSet& V, SplitVertexEdgeMap& E,
                     SplitVertexEdgeMap& E_P, SplitVertex& v_i, const VertexStatusMap& status ) {
    SplitVertex s_1 = V.at(v_i.s_1),
                   s_m = get_s_m( SG, V, E, v_i ),
                   s_j,
                   s_k;

//        cout <<  " v_i:" << v_i
//             << " s_1:" << s_1
//             <<"\n";
    if( E_P.empty() ) {
//             cout<< " adding edges between s_1 and s_m...\n";
        add_polygon_spanner_edges( SG, V, E, E_P, s_1, v_i, s_m, status );
    } else {
        optional< key_type > s_j_key;
        key_type s_k_key = s_m.key;
        foreach_neighbor( SG, V, E_P, v_i, [&]( SplitVertex& v_n ) {
            if( !s_j_key ) s_j_key = { v_n.key };
            s_k_key = v_n.key;
        } );
        assert(s_j_key);
        s_j =  V.at(*s_j_key );
        s_k =  V.at( s_k_key );
//            cout<< " s_j:" << s_j << " s_k:" << s_k <<"\n";
//            cout<< " adding edges between s_1 and s_j...\n";
        add_polygon_spanner_edges( SG, V, E, E_P, s_1, v_i, s_j, status );
//            cout<< " adding edges between s_k and s_m...\n";
        add_polygon_spanner_edges( SG, V, E, E_P, s_k, v_i, s_m, status );
    }

//        cout<< " s_m:" << s_m << "\n";
//        cout<< " adding edges from P...\n";
    add_polygon_edges( SG, E_P, s_1, v_i, s_m, status );
}

} // namespace polygon_spanner

void PolygonSpanner( DelaunayGraph& SG, SplitVertexSet& V, SplitVertexEdgeMap& E ) {
    using namespace polygon_spanner;

    // Create a vertex status map
    VertexStatusMap status;

    queue< key_type > level; // BFS queue
    SplitVertexEdgeMap E_P;

    SplitVertex v_i = V.V.front();

    level.push( v_i.key );

    //SG.addToEventQueue( v_i, 0 );

    do { // loop through level queue
        v_i = V.at( level.front() );
//        //SG.addToEventQueue( v_i, 0 ); // focus0 on v_i
//        cout<<"\nprocessing "<<v_i<<"\n";

        if( !E_P.empty() ) assert( E_P.at( v_i.key ).size() <= 5 );

        process_vertex( SG, V, E, E_P, v_i, status );

        // BFS housekeeping

        foreach_neighbor( SG, V, E, v_i, [&]( SplitVertex v_n ) {
            if( !contains( status, v_n.key ) ) { // If N[v_i] is NOT known, queue it and add to known
                level.push(v_n.key);
                status[v_n.key] = known;
//                SG.addToEventQueue( N, 1 ); // focus1 on C
            }
        });
        level.pop();
        status[v_i.key] = complete;
    } while( !level.empty() ); // level is not empty

    //SG.addToEventQueue( SG._DT.infinite_vertex(), 0 ); // focus0 on infinite
    Vertex_handle v1;

    // Add all edges from E_P to SG
    for( auto& edge : E_P ) {
        // Each edge is a pair of size_t, unordered_set<size_t>
        v1 = V.at( edge.first ).v;
        for( auto v2 : edge.second ) {
            SG.add_edge( v1, V.at(v2).v );
        }
    }
    // Lemma 3.4   ((PI+1)*(2*PI/(3*cos(PI/6)))) = 10.01602416
    // CAUTION: StretchFactor(SG) takes O(n^3)!
    //assert( StretchFactor(SG) <= ((PI+1)*(2*PI/(3*cos(PI/6)))) );

    // Test degree assumption given after lemma 3.4
    for( auto it=SG._E.begin(); it!=SG._E.end(); ++it ) {
        assert( it->second.size() <= 27 );
    }

    swap( E, E_P );

} // PolygonSpanner( SpanningGraph &P )

} // namespace bgs2005


template< typename RandomAccessIterator, typename OutputIterator >
void BGS2005( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false ) {
    using namespace bgs2005;

    DelaunayGraph G( pointsBegin, pointsEnd ); // Step 1

    SpanningGraph( G ); // Step 2

    SplitVertexSet V;
    SplitVertexEdgeMap P;
    TransformPolygon( G, V, P ); // Step 3

    PolygonSpanner( G, V, P ); // Step 4

    // send resulting edge list to output iterator
    for( auto const& adj : G._E ) {
        Vertex_handle v_1 = adj.first;
        for( auto const& v_2 : adj.second ) {
            *result = make_pair( v_1->info(), v_2->info() );
            ++result;
        }
    }
}

} // namespace gsnunf

#endif // GSNUNF_BGS2005_H
