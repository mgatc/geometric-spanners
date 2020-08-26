#ifndef GSNUNF_POLYGONSPANNER_H
#define GSNUNF_POLYGONSPANNER_H

#include <algorithm> // min, swap
#include <functional> // hash
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DelaunayGraph.h"
#include "GeometricSpannerPrinter.h"
#include "SpanningGraph.h"
#include "SplitVertex.h"
#include "StretchFactor.h"
#include "TransformPolygon.h"

namespace gsnunf {

using namespace std;

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

SplitVertex get_s_m( const DelaunayGraph& SG, const SplitVertexSet& V, const SplitVertexEdgeMap& E, const SplitVertex& v_split ) {
    key_type s_i;
    foreach_neighbor( SG, V, E, v_split, [&]( SplitVertex& v_n ) {
        s_i = v_n.key;
    });
    return V.at(s_i);
}

/*
 * Adds a split vertex edge, asserting that neither vertices a nor b are marked as complete
 */
void add_polygon_spanner_edge( const DelaunayGraph& SG, SplitVertexEdgeMap& E,
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
    using Vector_2 = typename DelaunayGraph::Vector_2;

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

void add_polygon_edges( const DelaunayGraph& SG, SplitVertexEdgeMap& E_P,
                        SplitVertex& p, SplitVertex& q, SplitVertex& r,
                        const VertexStatusMap& status ) {
    /*
     * Only add edges from P if they have not already been added.
     * Failure to check for this will result in breaking the
     * assumption that edges are only added to non-complete
     * vertices.
     */
    if( p.key >= E_P.size() || !contains( E_P.at(p.key), q.key ) )
        add_polygon_spanner_edge( SG, E_P, p, q, status );
    if( q.key >= E_P.size() || !contains( E_P.at(q.key), r.key ) )
        add_polygon_spanner_edge( SG, E_P, q, r, status );
}

void add_polygon_spanner_edges( const DelaunayGraph& SG, SplitVertexSet& V, SplitVertexEdgeMap& E, SplitVertexEdgeMap& E_P,
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
    for( size_t i=0; i<E_P.size(); ++i ) { //auto& edge : E_P ) {
        // Each edge is a pair of size_t, unordered_set<size_t>
        v1 = V.at(i).v;
        for( auto v2 : E_P.at(i) ) {
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

} // namespace gsnunf

#endif // GSNUNF_POLYGONSPANNER_H
