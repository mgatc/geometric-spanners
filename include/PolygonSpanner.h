#ifndef GSNUNF_POLYGONSPANNER_H
#define GSNUNF_POLYGONSPANNER_H

#include <algorithm>
#include <queue>

#include "DelaunayGraph.h"
#include "GeometricSpannerPrinter.h"
#include "SpanningGraph.h"
#include "StretchFactor.h"



namespace gsnunf {

namespace polygon_spanner {

    enum VertexStatus { unknown, known, complete };
    template< class T >
    using VertexStatusMap = std::unordered_map< SplitVertex<T>, VertexStatus, SplitVertexHasher >;
    //VertexMap< T, VertexMap< T, VertexStatus > >;


    /* Given a list of incident split vertices and a Vertex handle,
     * return the first split vertex with the given handle as its v.
     */
    template <typename T>
    SplitVertex<T> find_split( const IncidentSplitVertexContainer<T>& I, const Vertex_handle<T> v_n ) {
        for( SplitVertex<T> v : I )
            if( v.v == v_n )
                return v;
    }

    /*
     * Iterates from s_1 to s_m (inclusive) and performs "foreach" on each vertex.
     * Provides a vertex handle of the currently pointed-to vertex.
     */
    template< typename T, typename F >
    void foreach_neighbor( DelaunayGraph<T>& SG, SplitVertexMap<T>& V, const SplitVertexEdgeMap<T>& E, const SplitVertex<T>& v_i, F foreach ) {
        assert( !SG._DT.is_infinite( v_i.v ) );
        SplitVertex<T>* v_ptr;
        v_ptr = &V.at(v_i.v).at(v_i.s_1->v);

        IncidentSplitVertexContainer<T> incident_to_v_i = E.at( v_ptr ); // copy so we can make local edits

        Vertex_circulator<T> N = SG._DT.incident_vertices( v_i.v );
        N = SG.orient_circulator( N, v_i.s_1->v );
        SplitVertex<T> v_s = *v_i.s_1;
        incident_to_v_i.erase(v_s);

        size_t encounteredEdgesFromP = 0;

        do {
                //cout<< " N: "<<N->point()<<endl;
            if( encounteredEdgesFromP > 0 ) // if this is not s_1
                v_s = find_split( incident_to_v_i, N );

            assert( !SG._DT.is_infinite( N->handle() ) );
            if( gsnunf::contains( SG._E, v_i.v ) && gsnunf::contains( SG._E.at(v_i.v), N->handle() ) )
                ++encounteredEdgesFromP;
            //cout<<"v_s:"<<v_s.v->point()<<"("<<v_s.s_1->v->point()<<")\n";
            foreach( v_s );
            --N;
        } while( encounteredEdgesFromP < 2 );
    }

//    template< typename T >
//    Vertex_handle<T> find_s_1( const DelaunayGraph<T>& SG, const SplitVertexEdgeMap<T>& E,
//                               const Vertex_handle<T>& v, const Vertex_handle<T>& s ) {
//        //cout<<"find_s_1:"<<v->point()<<" ("<<s->point()<<")"<<endl;
//        Vertex_circulator<T> N = SG._DT.incident_vertices(v);
//        N = SG.orient_circulator( N, s ); //while( (--N_N)->handle() != v );// cout<<"until "<<v_i->point()<<": "<<N_N->point()<<"\n"; // rotate until N_N points to v_i
//        cout<<" "<<N->point()<<endl;
//        while( !contains( E.at(v), N->handle() ) ) cout<<" "<<(++N)->point()<<endl; // rotate until N_N points to a vertex in E[N]
//        //cout<<"  done"<<endl;
//        return N->handle();
//    }

//    template< typename T >
//    SplitVertex<T> split( const DelaunayGraph<T>& SG, const SplitVertexEdgeMap<T>& E,
//                          const Vertex_handle<T>& v, const Vertex_handle<T>& s ) {
//        //cout<<"split:"<<v->point()<<" ("<<s->point()<<")"<<endl;
//
//        return { v, find_s_1( SG, E, v, s ) };
//    }

//    template< typename T >
//    Vertex_handle<T> find_s_m( const DelaunayGraph<T>& SG, const SplitVertex<T>& v_split ) {
//        Vertex_handle<T> last;
//        foreach_neighbor( SG, v_split, [&]( Vertex_handle<T> N ) {
//            last = N;
//        } );
//        return last;
//    }

    template< typename T >
    SplitVertex<T> get_s_m( DelaunayGraph<T>& SG, SplitVertexMap<T>& V, const SplitVertexEdgeMap<T>& E,
                            SplitVertex<T>& v_split ) {
        SplitVertex<T> s_i;
        foreach_neighbor( SG, V, E, v_split, [&]( SplitVertex<T> v_n ) {
            s_i = v_n;
            //cout<< "s_i = "<<v_n.v->point()<<"("<<v_n.s_1->v->point()<<")\n";
        } );
        return s_i;
    }

    template< typename T >
    void add_cross_edges( DelaunayGraph<T>& SG, SplitVertexMap<T>& V, SplitVertexEdgeMap<T>& DT_split, SplitVertexEdgeMap<T>& E_P,
                          SplitVertex<T> p, SplitVertex<T> q, SplitVertex<T> r ) {
        SplitVertex<T> v_last;
        foreach_neighbor( SG, V, DT_split, q, [&]( SplitVertex<T> v_n ) {
            if( v_n != p ) {
                add_edge<T>( E_P, &v_last, &v_n );
                //SG.addToEventQueue( {s_last.first, N}, true );
            }
            v_last = v_n;
        } );
    }

    template< typename T >
    void add_forward_edges( DelaunayGraph<T>& SG, SplitVertexMap<T>& V, SplitVertexEdgeMap<T>& DT_split, SplitVertexEdgeMap<T>& E_P,
                            SplitVertex<T> p, SplitVertex<T> q, SplitVertex<T> r ) {
        using Vector_2 = typename DelaunayGraph<T>::Vector_2;

        double deg = 180/PI;

        double alpha = SG.get_angle( p.v, q.v, r.v );
        //cout<< "alpha " << alpha*deg;
        short subangles = rint( ceil( alpha / (PI/2) ) );
        //cout<< " subangles "<<subangles;
        double beta = alpha / subangles;
        //cout<< " beta "<<beta*deg<<endl;

        std::vector< SplitVertex<T> > add( subangles, SplitVertex<T>( SG._DT.infinite_vertex(), nullptr ) ); // initialize add to infinite vertex

        double theta;
        short i;

        foreach_neighbor( SG, V, DT_split, q, [&]( SplitVertex<T> v_n ) {
            //SG.addToEventQueue( N, 1 );// focus1 on N
            theta = SG.get_angle( p.v, q.v, v_n.v );
            if( theta > 2*PI-EPSILON )
                theta = 0;
            i = std::min( int(theta/beta), subangles-1 );
            //assert( i < subangles );
            //cout<<"  theta:"<<theta*deg<<" i:"<<i<<"\n";

            if( SG._DT.is_infinite( add.at(i).v )
              || Vector_2( v_n.v->point(), q.v->point() ).squared_length() < Vector_2( add.at(i).v->point(), q.v->point() ).squared_length() )
                add.at(i) = v_n;   // if the saved vertex is infinite or longer than the current one, update
        } );

        for( SplitVertex<T> v : add )
            if( !SG._DT.is_infinite( v.v ) && !SG._DT.is_infinite( v.s_1->v ) ) {
                add_edge<T>( E_P, &q, &v );
                //SG.addToEventQueue( {q.first, v.first}, true ); // add edge
            }
    }
    template< typename T >
    void add_polygon_edges( SplitVertexEdgeMap<T>& E_P,
                            SplitVertex<T>& p, SplitVertex<T>& q, SplitVertex<T>& r ) {
        add_edge<T>( E_P, &p, &q );
        add_edge<T>( E_P, &q, &r );
    }

    template< typename T >
    void add_polygon_spanner_edges( DelaunayGraph<T>& SG, SplitVertexMap<T>& V, SplitVertexEdgeMap<T>& DT_split, SplitVertexEdgeMap<T>& E_P,
                                    SplitVertex<T> p, SplitVertex<T> q, SplitVertex<T> r ) {
        assert( !SG._DT.is_infinite(p.v) && !SG._DT.is_infinite(q.v) && !SG._DT.is_infinite(r.v) );

        //cout<<"add polygon spanner edges\n";
        //cout<<r->point()<<endl;
        if( p == r ) return;

        add_forward_edges( SG, V, DT_split, E_P, p, q, r );
        add_cross_edges( SG, V, DT_split, E_P, p, q, r );
    }

    template< typename T >
    void process_vertex( DelaunayGraph<T>& SG, SplitVertexMap<T>& V, SplitVertexEdgeMap<T>& DT_split, SplitVertexEdgeMap<T>& E,
                         SplitVertexEdgeMap<T>& E_P, SplitVertex<T>& v_i ) {
        SplitVertex<T> s_1 = *v_i.s_1,
                       s_m = get_s_m( SG, V, DT_split, v_i ),
                       s_j = s_m,
                       s_k = s_m;
        if( E_P.empty() ) {
            add_polygon_spanner_edges( SG, V, DT_split, E_P, s_1, v_i, s_m );
        } else {
            SplitVertex<T> s_j = s_m,
                           s_k = s_m;
            bool first = true;
            foreach_neighbor( SG, V, E_P, v_i, [&]( SplitVertex<T> v_n ) {
                if( first ) {
                    s_j = v_n;
                    first = false;
                }
                s_k = v_n;
            } );
            add_polygon_spanner_edges( SG, V, DT_split, E_P, s_1, v_i, s_j );
            add_polygon_spanner_edges( SG, V, DT_split, E_P, s_k, v_i, s_m );
        }
        add_polygon_edges<T>( E_P, s_1, v_i, s_m );

        //cout<<"    sm: "<< s_m.v->point() << " ("<<s_m.s_1->v->point()<<")\n";
//        Vertex_circulator<T> unsplit_s_j = SG._DT.incident_vertices( v_i.first ),
//                             unsplit_s_k = SG._DT.incident_vertices( v_i.first );
//        cout<<"processing "<<v_i.first->point()<<endl;
//
//        if( E_P.size() == 0 ) { // processing v_1
//            add_polygon_spanner_edges( E, E_P, SG, s_1, v_i, s_m );
//        } else { // processing v_i, where i>1
//            auto& N_PS = E_P.at(v_i.first).at(s_1.first);  // neighbors of v_i in E_P edges
//            cout<<N_PS.size()<<"\n";
//            assert( N_PS.size() <= 5 ); // Lemma 3.3
//            // find s_j
//            unsplit_s_j = SG.orient_circulator( unsplit_s_j, s_1.first );
//            while( !contains( N_PS, unsplit_s_j->handle() ) ) --unsplit_s_j;
//            SplitVertex<T> s_j = split( SG, E, unsplit_s_j, v_i.first );
//            add_polygon_spanner_edges( E, E_P, SG, s_1, v_i, s_j );
//            // find s_k
//            unsplit_s_k = SG.orient_circulator( unsplit_s_k, s_m.first );
//            while( !contains( N_PS, unsplit_s_k->handle() ) ) ++unsplit_s_k;
//            SplitVertex<T> s_k = split( SG, E, unsplit_s_k, v_i.first );
//            add_polygon_spanner_edges( E, E_P, SG, s_k, v_i, s_m );
//            cout<<"s_1:"<<s_1.first->point()<<" ("<<s_1.second->point()<<") "
//                <<"s_j:"<<s_j.first->point()<<" ("<<s_j.second->point()<<") "
//                <<"s_k:"<<s_k.first->point()<<" ("<<s_k.second->point()<<") "
//                <<"s_m:"<<s_m.first->point()<<" ("<<s_m.second->point()<<") "<<endl;
//        }
    }

    template< class T >
    Vertex_handle<T> find_s_1_handle( DelaunayGraph<T>& SG, const pair< const Vertex_handle<T>, VertexMap<T, SplitVertex<T> > >& unsplit, const Vertex_handle<T>& v_n ) {
        Vertex_handle<T> v_i = unsplit.first;
        Vertex_circulator<T> N = SG._DT.incident_vertices(v_i); // get circulator around unsplit.first
        while( (++N)->handle() != v_n );
        while( !contains( unsplit.second, N->handle() ) ) ++N; // orient to a neighbor in unsplit.second

        return N;
    }

    template< class T >
    void TransformDelaunayTriangulation( DelaunayGraph<T>& SG, SplitVertexMap<T>& V, SplitVertexEdgeMap<T>& P ) {
        TransformPolygon( SG, V, P );

        for( auto e=SG._DT.finite_edges_begin(); e!=SG._DT.finite_edges_end(); ++e ) {
            Vertex_handle<T> u = e->first->vertex( (e->second+1)%3 ),
                             v = e->first->vertex( (e->second+2)%3 );
            pair< const Vertex_handle<T>, VertexMap<T, SplitVertex<T> > >& unsplit_u = *V.find( u ),
                                                              unsplit_v = *V.find( v );
            Vertex_handle<T> s_1_u = find_s_1_handle( SG, unsplit_u, v ),
                             s_1_v = find_s_1_handle( SG, unsplit_v, u );
            add_edge<T>( P, &V[u][s_1_u], &V[v][s_1_v] );
        }
//        cout<<endl;
//        print_edges<T>(P);
//        cout<<endl;
    }
} // namespace polygon_spanner

template< typename T >
void PolygonSpanner( DelaunayGraph<T>& SG, SplitVertexMap<T>& V, SplitVertexEdgeMap<T>& E ) {

    using namespace polygon_spanner;
    SplitVertexEdgeMap<T> DT_split;
    TransformDelaunayTriangulation( SG, V, DT_split );

    // Create a vertex status map
    VertexStatusMap<T> status;

    std::queue< SplitVertex<T> > level; // BFS queue
    //DelaunayGraph PS( SG._DT );         // resultant graph object
    SplitVertexEdgeMap<T> E_P;

    Vertex_handle<T> v_1_handle = SG._DT.finite_vertices_begin();
    SplitVertex<T> v_i = V[v_1_handle].begin()->second;

    level.push( v_i );

    //SG.addToEventQueue( v_i, 0 );

    do { // loop through level queue
        v_i = level.front();
        //SG.addToEventQueue( v_i, 0 ); // focus0 on v_i
        cout<<"processing "<<v_i.v->point()<<" s1: "<<v_i.s_1->v->point()<<"\n";

        process_vertex( SG, V, DT_split, E, E_P, v_i );

        // BFS housekeeping

        foreach_neighbor( SG, V, DT_split, v_i, [&]( SplitVertex<T> v_n ) {
            if( !contains( status, v_n ) ) { // If N[v_i] is NOT known, queue it and add to known
                level.push(v_n);
                status[v_n] = known;
//                SG.addToEventQueue( N, 1 ); // focus1 on C
            }
        } );

        level.pop();
        status[v_i] = complete;
    } while( !level.empty() ); // level is not empty

    //SG.addToEventQueue( SG._DT.infinite_vertex(), 0 ); // focus0 on infinite


    // Lemma 3.4
    // ((PI+1)*(2*PI/(3*cos(PI/6)))) = 10.01602416
    //assert( StretchFactor(SG) <= (PI+1) ); // fails
    //assert( StretchFactor(SG) <= ((PI+1)*(2*PI/(3*cos(PI/6)))) ); // fails

    // Test degree assumption given after lemma 3.4
//    for( auto it=SG._E.begin(); it!=SG._E.end(); ++it ) {
//        assert( it->second.size() <= 12 );
//    }

} // PolygonSpanner( SpanningGraph &P )

} // namespace gsnunf

#endif // GSNUNF_POLYGONSPANNER_H
