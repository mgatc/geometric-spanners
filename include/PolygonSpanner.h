#ifndef GSNUNF_POLYGONSPANNER_H
#define GSNUNF_POLYGONSPANNER_H

#include <algorithm>
#include <queue>

#include "DelaunayGraph.h"
#include "GeometricSpannerPrinter.h"
#include "SpanningGraph.h"
#include "SplitVertex.h"
#include "StretchFactor.h"
#include "TransformPolygon.h"



namespace gsnunf {

namespace polygon_spanner {

    enum VertexStatus { unknown, known, complete };
    template< class T >
    using VertexStatusMap = std::unordered_map< key_type<T>, VertexStatus >;

    /*
     * Iterates from s_1 to s_m (inclusive) and performs "foreach" on each vertex.
     * Provides a SplitVertex of the currently pointed-to vertex.
     */
    template< typename T, typename F >
    void foreach_neighbor( DelaunayGraph<T>& SG, SplitVertexSet<T>& V, const SplitVertexEdgeMap<T>& E, const SplitVertex<T>& v_i, F foreach ) {
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
            std::size_t operator()( const SplitVertex<T>& k ) const {
                return hash< Vertex_handle<T> >()( k.v );
            }
        };
        struct CompareOnVertexHandle {
            bool operator()( const SplitVertex<T>& a, const SplitVertex<T>& b ) const {
                return a.v == b.v;
            }
        };
        unordered_set< SplitVertex<T>, HashOnVertexHandle, CompareOnVertexHandle > N_E;

        for( auto k: E.at( v_i.key ) ) {
            if( k != v_i.s_1 ) { // don't enter s_1
                N_E.insert( V.at(k) );
            }
        }
        //cout<<"N_E size:"<<N_E.size()<<"\n";

        Vertex_circulator<T> N = SG._DT.incident_vertices( v_i.v );
        N = SG.orient_circulator( N, v_i.s_1_handle );


        SplitVertex<T> v_neighbor = V.at( v_i.s_1 );
        bool include_s_1 = contains( E.at(v_i.key), v_neighbor.key );
        //cout<< "    start foreach: "<<v_i<<" n: "<<(N_E.size()+size_t(include_s_1))<<"\n";

        if( include_s_1 ) {
        //    cout<<"    v_neighbor: "<<v_neighbor<<"\n";
            foreach( v_neighbor );
        }

        do {
            --N;
            v_neighbor = SplitVertex<T>( N->handle() );
            if( contains( N_E, v_neighbor ) ) {
                v_neighbor = *N_E.find( v_neighbor );
                //cout<<"    v_neighbor: "<<v_neighbor<<"\n";
                foreach( v_neighbor );
                N_E.erase( v_neighbor );
            }
        } while( !N_E.empty() );
    }

    template< typename T >
    SplitVertex<T> get_s_m( DelaunayGraph<T>& SG, SplitVertexSet<T>& V, const SplitVertexEdgeMap<T>& E,
                            SplitVertex<T>& v_split ) {
        key_type<T> s_i;
        foreach_neighbor( SG, V, E, v_split, [&]( SplitVertex<T>& v_n ) {
            s_i = v_n.key;
        });
        return V.at(s_i);
    }

    /*
     * Adds a split vertex edge, asserting that neither vertices a nor b are marked as complete
     */
    template< typename T >
    void add_polygon_spanner_edge( const DelaunayGraph<T>& SG, SplitVertexEdgeMap<T>& E,
                                  SplitVertex<T>& a, SplitVertex<T>& b,
                                  const VertexStatusMap<T>& status ) {
        assert( !contains( status, a.key ) || status.at(a.key) != complete );
        assert( !contains( status, b.key ) || status.at(b.key) != complete );
        return add_edge<T>( SG, E, a, b );
    }

    template< typename T >
    void add_cross_edges( DelaunayGraph<T>& SG, SplitVertexSet<T>& V, SplitVertexEdgeMap<T>& E, SplitVertexEdgeMap<T>& E_P,
                          SplitVertex<T>& p, SplitVertex<T>& q, SplitVertex<T>& r,
                          const VertexStatusMap<T>& status ) {
        optional< SplitVertex<T> > v_last = nullopt;
        bool isInZone = false;
        foreach_neighbor( SG, V, E, q, [&]( SplitVertex<T>& v_n ) {
            if( isInZone && v_last ) {
                add_polygon_spanner_edge<T>( SG, E_P, *v_last, v_n, status );
                //SG.addToEventQueue( {s_last.first, N}, true );
            }
            if( v_n.key == p.key ) isInZone = true;
            if( v_n.key == r.key ) isInZone = false;
            //cout<<(isInZone ? "in zone" : "not in zone" )<<"\n";
            //cout<<"  v_n "<<v_n<<" v_last "<<v_last<<"\n";
            v_last = { v_n };
        });
    }

    template< typename T >
    void add_forward_edges( DelaunayGraph<T>& SG, SplitVertexSet<T>& V, SplitVertexEdgeMap<T>& E, SplitVertexEdgeMap<T>& E_P,
                            SplitVertex<T>& p, SplitVertex<T>& q, SplitVertex<T>& r,
                            const VertexStatusMap<T>& status ) {
        using Vector_2 = typename DelaunayGraph<T>::Vector_2;

        double deg = 180/PI;

        double alpha = SG.get_angle( p.v, q.v, r.v );
        short subangles = rint( ceil( alpha / (PI/2) ) );
        double beta = alpha / subangles;

        std::vector< SplitVertex<T> > add( subangles, SplitVertex<T>( SG._DT.infinite_vertex() ) ); // initialize add to infinite vertex

        double theta;
        short i;
        bool isInZone = false;

        foreach_neighbor( SG, V, E, q, [&]( SplitVertex<T>& v_n ) {
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

        for( SplitVertex<T> v : add )
            if( !SG._DT.is_infinite( v.v ) ) {
                add_polygon_spanner_edge<T>( SG, E_P, q, v, status );
                //SG.addToEventQueue( {q.first, v.first}, true ); // add edge
            }
    }
    template< typename T >
    void add_polygon_edges( DelaunayGraph<T>& SG, SplitVertexEdgeMap<T>& E_P,
                            SplitVertex<T>& p, SplitVertex<T>& q, SplitVertex<T>& r,
                            const VertexStatusMap<T>& status ) {
        /*
         * Only add edges from P if they have not already been added.
         * Failure to check for this will result in breaking the
         * assumption that edges are only added to non-complete
         * vertices.
         */
        if( !contains( E_P, p.key) || !contains( E_P.at(p.key), q.key ) )
            add_polygon_spanner_edge<T>( SG, E_P, p, q, status );
        if( !contains( E_P, q.key) || !contains( E_P.at(q.key), r.key ) )
            add_polygon_spanner_edge<T>( SG, E_P, q, r, status );
    }

    template< typename T >
    void add_polygon_spanner_edges( DelaunayGraph<T>& SG, SplitVertexSet<T>& V, SplitVertexEdgeMap<T>& E, SplitVertexEdgeMap<T>& E_P,
                                    SplitVertex<T>& p, SplitVertex<T>& q, SplitVertex<T>& r,
                                    const VertexStatusMap<T>& status ) {
        assert( !SG._DT.is_infinite(p.v) && !SG._DT.is_infinite(q.v) && !SG._DT.is_infinite(r.v) );
        if( p == r ) return;

            //cout<< " adding forward edges...\n";

        add_forward_edges( SG, V, E, E_P, p, q, r, status );
                    //cout<< " adding cross edges...\n";

        add_cross_edges( SG, V, E, E_P, p, q, r, status );
    }

    template< typename T >
    void process_vertex( DelaunayGraph<T>& SG, SplitVertexSet<T>& V, SplitVertexEdgeMap<T>& E,
                         SplitVertexEdgeMap<T>& E_P, SplitVertex<T>& v_i, const VertexStatusMap<T>& status ) {
        SplitVertex<T> s_1 = V.at(v_i.s_1),
                       s_m = get_s_m( SG, V, E, v_i ),
                       s_j,
                       s_k;

//        cout <<  " v_i:" << v_i
//             << " s_1:" << s_1
//             <<"\n";
        //cout << "E_P:\n";
        //print( V, E_P );
        if( E_P.empty() ) {
//             cout<< " adding edges between s_1 and s_m...\n";
            add_polygon_spanner_edges( SG, V, E, E_P, s_1, v_i, s_m, status );
        } else {
            optional< key_type<T> > s_j_key;
            key_type<T> s_k_key = s_m.key;
//            cout<< " E_P.at( v_i.key ).size: "<<E_P.at( v_i.key ).size()<<"\n";
            foreach_neighbor( SG, V, E_P, v_i, [&]( SplitVertex<T>& v_n ) {
                if( !s_j_key ) s_j_key = { v_n.key };
                s_k_key = v_n.key;
//                cout<<" s_j:" << V.at(*s_j_key) << " s_k:" << V.at(s_k_key)<<"\n";

            } );
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
        add_polygon_edges<T>( SG, E_P, s_1, v_i, s_m, status );
    }

} // namespace polygon_spanner

template< typename T >
void PolygonSpanner( DelaunayGraph<T>& SG, SplitVertexSet<T>& V, SplitVertexEdgeMap<T>& E ) {

    using namespace polygon_spanner;

    // Create a vertex status map
    VertexStatusMap<T> status;

    std::queue< key_type<T> > level; // BFS queue
    SplitVertexEdgeMap<T> E_P;

    SplitVertex<T> v_i = V.V.front();

    level.push( v_i.key );

    //SG.addToEventQueue( v_i, 0 );

    do { // loop through level queue
        v_i = V.at( level.front() );
//        //SG.addToEventQueue( v_i, 0 ); // focus0 on v_i
//        cout<<"\nprocessing "<<v_i<<"\n";

        if( !E_P.empty() ) assert( E_P.at( v_i.key ).size() <= 5 );

        process_vertex( SG, V, E, E_P, v_i, status );

        // BFS housekeeping

        foreach_neighbor( SG, V, E, v_i, [&]( SplitVertex<T> v_n ) {
            if( !contains( status, v_n.key ) ) { // If N[v_i] is NOT known, queue it and add to known
                level.push(v_n.key);
                status[v_n.key] = known;
//                SG.addToEventQueue( N, 1 ); // focus1 on C
            }
        } );

        level.pop();
        status[v_i.key] = complete;
    } while( !level.empty() ); // level is not empty

    //SG.addToEventQueue( SG._DT.infinite_vertex(), 0 ); // focus0 on infinite
    Vertex_handle<T> v1;

    // Add all edges from E_P to SG
    for( auto& edge : E_P ) {
        // Each edge is a pair of size_t, unordered_set<size_t>
        v1 = V.at( edge.first ).v;
        for( auto v2 : edge.second ) {
            SG.add_edge( v1, V.at(v2).v );
        }
    }
    // Lemma 3.4   ((PI+1)*(2*PI/(3*cos(PI/6)))) = 10.01602416
    assert( StretchFactor(SG) <= ((PI+1)*(2*PI/(3*cos(PI/6)))) );

    // Test degree assumption given after lemma 3.4
    for( auto it=SG._E.begin(); it!=SG._E.end(); ++it ) {
        assert( it->second.size() <= 12 );
    }

    swap( E, E_P );

} // PolygonSpanner( SpanningGraph &P )

} // namespace gsnunf

#endif // GSNUNF_POLYGONSPANNER_H
