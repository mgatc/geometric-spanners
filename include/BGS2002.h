#ifndef GSNUNF_BGS2002_H
#define GSNUNF_BGS2002_H

#include <algorithm> //sort
#include <iostream>
#include <list>
#include <vector>

#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/algorithm.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
//#include <boost/property_map.hpp>
#include <boost/ref.hpp>
#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>


#include "Timer.h"
//#include "DelaunayGraph.h"
//#include "SpanningGraph.h"
#include "TransformPolygon.h"
#include "PolygonSpanner.h"
#include "GeometricSpannerPrinter.h"

namespace gsnunf {

namespace bgs2002 {

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Vb;
typedef CGAL::Triangulation_face_base_2<K>                          Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef Delaunay::Point                                             Point;
typedef Delaunay::Finite_vertices_iterator Finite_vertices_iterator;
typedef Delaunay::Finite_edges_iterator    Finite_edges_iterator;

typedef  CGAL::Aff_transformation_2<K> Transformation;
typedef  CGAL::Vector_2<K> Vector;


typedef std::pair<unsigned,unsigned> unsignedPair;

typedef boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    boost::property< boost::vertex_index_t, size_t >,
    boost::property< boost::edge_index_t, size_t >
> AdjacencyList;

template<class ArgumentType, class ResultType>
struct unary_funct {
    typedef ArgumentType argument_type;
    typedef ResultType result_type;
};

struct AutoCount : public unary_funct<const Point&,std::pair<Point,unsigned> > {
    mutable unsigned i;
    AutoCount() : i(0) {}
    pair<Point,unsigned> operator()(const Point& p) const {
        return make_pair(p,i++);
    }
};

namespace spanning_graph {

void add_first_edge( DelaunayGraph& G, Vertex_handle v, Vertex_circulator C ) {
    Vertex_handle v2 = C->handle();
    G.add_edge( v, v2 );

    //G.addToEventQueue( v, 1 );
    //G.addToEventQueue( v2, 2 );
    //G.addToEventQueue( { v, v2 }, true );
}

void add_second_edge( DelaunayGraph& G, Vertex_handle v, Vertex_circulator C ) {
    while( G._DT.is_infinite(++C) );
    Vertex_handle v2 = C->handle();
    G.add_edge( v, v2 );

//        G.addToEventQueue( v, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v, v2 }, true );
}

void add_last_edge( DelaunayGraph& G, Vertex_handle v, Vertex_circulator C, const VertexHash& is_removed ) {
    --C;
    Vertex_circulator done(C);

    while( ( contains( is_removed, C ) || G._DT.is_infinite(C) ) && --C != done );

    Vertex_handle v2 = C->handle();

    G.add_edge( v, v2 );

//        G.addToEventQueue( v, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v, v2 }, true );
}

void remove_first_edge( DelaunayGraph& G, Vertex_circulator C ) {
    Vertex_handle v1 = C->handle(),
                  v2 = (++C)->handle();
    G.remove_edge( v1, v2 );

//        G.addToEventQueue( v1, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v1, v2 }, false );
}

void remove_second_edge( DelaunayGraph& G, Vertex_circulator C ) {
    Vertex_handle v1 = (++C)->handle(),
                  v2 = (++C)->handle();
    G.remove_edge( v1, v2 );

//        G.addToEventQueue( v1, 1 );
//        G.addToEventQueue( v2, 2 );
//        G.addToEventQueue( { v1, v2 }, false );
}

void remove_last_edge( DelaunayGraph& G, Vertex_circulator C, const VertexHash& is_removed ) {
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

void SpanningGraph( const Delaunay& DT, AdjacencyList& G ) {
    using namespace spanning_graph;
    using namespace boost;

    size_t N = DT.number_of_vertices();

    AdjacencyList G_temp(N);
    // Add all edges to graph
    for( auto e = DT.finite_edges_begin(); e!=DT.finite_edges_end(); ++e ) {
        auto v_1 = e->first->vertex( (e->second+1)%3 );
        auto v_2 = e->first->vertex( (e->second+2)%3 );
        cout<< v_1->info()<< " " << v_2->info() << "\n";
        add_edge( v_1->info(), v_2->info(), G_temp );

    }
    // Initialize the interior edge index
    property_map<AdjacencyList, edge_index_t>::type e_index = get( edge_index, G_temp );
    graph_traits<AdjacencyList>::edges_size_type edge_count = 0;
    graph_traits<AdjacencyList>::edge_iterator ei, ei_end;

    for( tie(ei, ei_end) = edges(G_temp); ei != ei_end; ++ei )
        put(e_index, *ei, edge_count++);

    // Test for planarity - we know it is planar, we just want to
    // compute the planar embedding as a side-effect
    typedef vector< graph_traits<AdjacencyList>::edge_descriptor > vec_t;
    vector<vec_t> embedding(N);
    if( !boyer_myrvold_planarity_test(
            boyer_myrvold_params::graph = G_temp,
            boyer_myrvold_params::embedding = &embedding[0]
        )
    ) cout << "Error: Input graph is not planar\n";

    vector<graph_traits<AdjacencyList>::vertex_descriptor> canonical;
    planar_canonical_ordering( G_temp, &embedding[0], std::back_inserter(canonical) );

//    Vertex_circulator v_n, done;
    size_t i;

    vector<bool> is_removed( N, true );

    cout<<"canonical size:"<<canonical.size()<<"\n";
    for( size_t i : canonical ) {
        cout<< i<<"\n";
    }

    // Add first three vertices from canonical
    for( i=0; i<3; ++i ) { // Add edges of triangle
        is_removed[canonical.at(i)] = false;
        add_edge( canonical.at(i), canonical.at((i+1)%3), G );
        //G.addToEventQueue( { canonical.at(i), canonical.at((i+1)%3) }, true );
        //G.addToEventQueue( { canonical.at(i), canonical.at((i+1)%3) }, true );
    }
    // Add the rest of the vertices from canonical
    for( i=i; i<N; ++i ) {
        assert( canonical.at(i) < is_removed.size() );
        is_removed[canonical.at(i)] = false;
//        G.addToEventQueue( *c_iter, 0 );        // activate c_iter
//        G.addToEventQueue( *c_iter, true );

//        v_n = G._DT.incident_vertices( canonical.at(i) );
//        done = v_n;
//
//        G.normalize_circulator( v_n, is_removed );
//        done = v_n;
//
//        int k = G.count_valid_neighbors( v_n, is_removed );
//
//        if( k == 2 ) {
//            // remove edge between first two vertices
//            remove_first_edge( G, v_n );
//            // add edge between canonical iterator and first vertex
//            add_first_edge( G, canonical.at(i), v_n );
//            // add edge between canonical iterator and second vertex
//            add_second_edge( G, canonical.at(i), v_n );
//
//        } else if( k > 2 ) {
//            // remove edge between first two vertices
//            remove_first_edge( G, v_n );
//            // remove edge between last two vertices
//            remove_last_edge( G, v_n, is_removed );
//            // add edge between canonical iterator and first vertex
//            add_first_edge( G, canonical.at(i), v_n );
//            // add edge between canonical iterator and second vertex
//            add_second_edge( G, canonical.at(i), v_n );
//            // add edge between canonical iterator and last vertex
//            add_last_edge( G, canonical.at(i), v_n, is_removed );
//        }
    }

//    // Test assumption
//    for( auto it=G._E.begin(); it!=G._E.end(); ++it ) { // for all v_i, 1<=i<=n
//        assert( it->second.size() <= 3 );                 // |v_i| <= 3
//    }
}

} // namespace bgs2002

template< typename RandomAccessIterator, typename OutputIterator >
void BGS2002( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result ) {
    using namespace bgs2002;
    //GeometricSpannerPrinter printer;

    //Timer t(",");
    cout<<"Delaunay\n";

    // Make Delaunay triangulation and assign IDs
    vector<Point> points( pointsBegin, pointsEnd );
    Delaunay DT;
    DT.insert(
        boost::make_transform_iterator( points.begin(), AutoCount() ),
        boost::make_transform_iterator(   points.end(), AutoCount() )
    );


//    printer.drawEdges( G._DT, {
//        { "line width", "2pt" },
//        { "color", "gray" }
//    });
    cout<<"SpanningGraph\n";

    AdjacencyList G( points.size() );

    SpanningGraph( DT, G );
//    printer.drawEdges( G, {
//        { "line width", "13pt" },
//        { "color", "cyan" }
//    });
//
    //size_t split_size_estimate = G._DT.number_of_vertices();
    //SplitVertexSet V;
    //SplitVertexEdgeMap P;
//    {
//        //Timer timer(",");
           // cout<<"TransformPolygon\n";
//
      //  TransformPolygon( G, V, P );
//    }
    //cout<<sizeof(*P.begin())<<"\n";
////    {
////        //Timer timer(",");
            //cout<<"PolygonSpanner\n";

     //   PolygonSpanner( G, V, P );

//        cout<<"V:"<< sizeof(*V.index.begin())*V.V.size()<<" "
//            <<"E:"<< sizeof(*P.begin())*P.size()<<" \n";
//    }
//    cout<<"\n";
//    printer.drawEdges( G, {
//        { "line width", "5pt" },
//        { "color", "cyan" }
//    });
//    printer.drawVertices( G._DT, 5, {
//        { "fill", "blue" },
//        { "draw", "white" }
//    });
//    printer.print( "PolygonSpanner" );
//    cout<<StretchFactor(G)<<",";
//    cout<<G.degree()<<",";

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

#endif // GSNUNF_BGS2002_H
