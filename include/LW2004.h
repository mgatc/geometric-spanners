#ifndef GSNUNF_LW2004_H
#define GSNUNF_LW2004_H

#include <iostream>
#include <optional>
#include <vector>

#include "Timer.h"
#include "DelaunayGraph.h"
#include "SpanningGraph.h"
#include "TransformPolygon.h"
#include "PolygonSpanner.h"
#include "GeometricSpannerPrinter.h"

namespace gsnunf {

namespace lw2004 {

void add_cross_edges( DelaunayGraph& Del, const Vertex_handle& v_1, const Vertex_handle& u, const Vertex_handle& v_k ) {

    Vertex_circulator N = Del._DT.incident_vertices(u);
    N = Del.orient_circulator( N, v_1 );
    Vertex_handle last = N,
        first = N;

    do {
        ++N;
        if( (last != first && N != v_k) || v_1 == v_k )
            Del.add_edge( N, last );
        last = N;
    } while( N != v_k );

//    while( ++N != v_k ) {
//        if( last )
//            Del.add_edge( N, *last );
//        last = { N->handle() };
//    }
//    // don't add v_1 or v_k unless v_1 == v_k
//    if( v_1 == v_k ) {
//        Del.add_edge( N, v_1 );
//        Del.add_edge( N, v_k );
//    }
}


void add_forward_edges( DelaunayGraph& Del, const Vertex_handle& v_1, const Vertex_handle& u, const Vertex_handle& v_k, double angle = PI/2 ) {

    using Vector_2 = typename DelaunayGraph::Vector_2;

    //double deg = 180/PI; // used for displaying angles in degree

    double alpha = Del.get_angle( v_1, u, v_k );
    short subangles = rint( ceil( alpha / angle ) );
    double beta = alpha / subangles;

    vector< Vertex_handle > add( subangles, Del._DT.infinite_vertex() ); // initialize add to infinite vertex

    double theta;
    short i;
    bool isInZone = false;

    Vertex_circulator N = Del._DT.incident_vertices(u);
    N = Del.orient_circulator( N, v_1 );
    optional<Vertex_handle> last = nullopt;

    while( ++N != v_k ) {
        // r is guaranteed to be already added or will be added immediately after this step
        // so set isInZone to false before processing it
        if( N->handle() == v_k ) isInZone = false;
        //SG.addToEventQueue( N, 1 );// focus1 on N
        if( isInZone ) {
            theta = Del.get_angle( v_1, u, N );
            if( theta > 2*PI-EPSILON )
                theta = 0;
            i = std::min( int(theta/beta), subangles-1 );

            if( Del._DT.is_infinite( add.at(i) )
              || Vector_2( N->point(), u->point() ).squared_length() < Vector_2( add.at(i)->point(), u->point() ).squared_length() )
                add.at(i) = N;   // if the saved vertex is infinite or longer than the current one, update
        }
        // p is guaranteed to be already added or will be added immediately after this step
        // so don't set isInZone to true until after it's passed
        if( N->handle() == v_1 ) isInZone = true;

    }

    for( Vertex_handle v : add )
        if( !Del._DT.is_infinite( v ) ) {
            Del.add_edge( u, v );
            //SG.addToEventQueue( {q.first, v.first}, true ); // add edge
        }
}

void add_edges( DelaunayGraph& Del, const Vertex_handle& v_1, const Vertex_handle& u, const Vertex_handle& v_k, double alpha = PI/2 ) {
    cout<<"  add edges between "<< v_1->point()<<" and "<<v_k->point()<<"\n";

    cout<<"   add forward edges\n";
    add_forward_edges( Del, v_1, u, v_k, alpha );
    cout<<"   add cross edges\n";
    add_cross_edges( Del, v_1, u, v_k );
}

}

template< typename RandomAccessIterator, typename OutputIterator >
void LW2004( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result ) {
    using namespace lw2004;
    //GeometricSpannerPrinter printer( .25f );
    //Timer t(",");

    // Step 1

    DelaunayGraph Del( pointsBegin, pointsEnd );
    Del.add_all_edges(); // add all edges from the underlying DT to the edge list of the graph

    // Step 2

    vector<Vertex_handle> pi( Del.n(), Del._DT.infinite_vertex() );
    auto order = pi.begin();

    using SingleVertexIncidenceListIterator = DelaunayGraph::AdjacencyList::iterator;
    struct SingleVertexIncidenceListCompare {
        bool operator()( const SingleVertexIncidenceListIterator& lhs, const SingleVertexIncidenceListIterator& rhs ) const {
            //bool less_than = true;
            return (
                lhs->second.size() > rhs->second.size()
              ||
                ( lhs->second.size() == rhs->second.size() && lhs->first->point() > rhs->first->point() )
            );
        }
    };
    std::set< SingleVertexIncidenceListIterator, SingleVertexIncidenceListCompare > G;
    for( auto it=Del._E.begin(); it != Del._E.end(); ++it ) {
        G.insert(it);
    }
    SingleVertexIncidenceListIterator u;
    while( !G.empty() ) {
        u = *G.begin();
        G.erase(u);
        *order = u->first;
        ++order;
        Vertex_circulator N = Del._DT.incident_vertices( u->first ),
            done(N);

        do {
            //cout<<"v:"<<N->point()<<"\n";
            SingleVertexIncidenceListIterator temp = Del._E.find( N->handle() );
            if( temp != Del._E.end() && !temp->second.empty() ) {
                //const SingleVertexIncidenceListIterator& temp_G = *G.find(temp_E);
                G.erase( temp );
                Del.remove_edge( N->handle(), u->first );
                G.insert( temp );
            }
        } while( --N != done );
        assert( Del._E[u->first].size() == 0 );
    }

    Del._E.clear();

    // Step 3

    VertexSet processed;

    for( Vertex_handle u : pi ) {
        assert( Del._E[u].size() <= 5 );
        cout<<u->point()<<"\n";

        Vertex_circulator N = Del._DT.incident_vertices(u);
        Vertex_handle done = N->handle(), last;
        // orient to a processed vertex
        while( !contains( processed, (--N)->handle() ) && N->handle() != done );
        done = N;
        last = N;

        do {
            --N;
            if( contains( processed, N->handle() ) || N->handle() == done ) {
                add_edges( Del, last, u, N->handle() );
                last = N;
            } else if( N->handle() == done ) {

            }
        } while( N->handle() != done );

        processed.insert(u);
    }

    //printer.print( "BGS2002" );

    // send resulting edge list to output iterator
//    for( auto const& adj : G._E ) {
//        Vertex_handle v_1 = adj.first;
//        for( auto const& v_2 : adj.second ) {
//            *result = make_pair( v_1->point(), v_2->point() );
//            ++result;
//        }
//    }
}

} // namespace gsnunf

#endif // GSNUNF_LW2004_H
