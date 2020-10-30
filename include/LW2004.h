#ifndef GSNUNF_LW2004_H
#define GSNUNF_LW2004_H

#include <algorithm> // min, max
#include <cmath> // ceil
#include <unordered_set> // hashed adjacency list
#include <vector> // vertex containers

#include <boost/functional/hash.hpp> // hashing pairs
#include <boost/heap/fibonacci_heap.hpp> // ordering

//#include <CGAL/algorithm.h>
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/utils.h> // min, max
#include <CGAL/Vector_2.h>

#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "utilities.h"

namespace gsnunf {

using namespace std;

namespace lw2004 {

using namespace CGAL;

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K>    Vb;
typedef CGAL::Triangulation_face_base_2<K>                          Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef CGAL::Aff_transformation_2<K>                               Transformation;
typedef Delaunay::Vertex_handle                                     Vertex_handle;
typedef Delaunay::Vertex_circulator                                 Vertex_circulator;
typedef CGAL::Vector_2<K>                                           Vector_2;
typedef Delaunay::Point                                             Point;
typedef Delaunay::Finite_vertices_iterator                          Finite_vertices_iterator;
typedef Delaunay::Finite_edges_iterator                             Finite_edges_iterator;

typedef pair<size_t,size_t>                                         size_tPair;
typedef boost::hash<size_tPair>                                     size_tPairHash;
typedef unordered_set<size_tPair,size_tPairHash>                    size_tPairSet;

struct comparatorForMinHeap {
    bool operator()(const size_tPair &n1, const size_tPair &n2) const {
        return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
    }
};

typedef boost::heap::fibonacci_heap<size_tPair,boost::heap::compare<comparatorForMinHeap>> Heap;
typedef Heap::handle_type handle;

inline size_tPair createEdge( const size_t i, const size_t j ) {
    return make_pair( std::min(i,j), std::max(i,j) );
}

inline void createNewEdge( const Delaunay& T, const vector<Delaunay::Vertex_handle>& handles, size_tPairSet &E, const size_t i, const size_t j, const size_t n, bool printLog = false ) {
    assert( std::max(i,j) < n );
    assert( T.is_edge( handles.at(i), handles.at(j) ) );
    //if( printLog ) cout<<"add:("<<i<<","<<j<<") ";
    E.insert( createEdge( i, j ) );
}

} // namespace lw2004

// alpha is set to pi/2
template< typename RandomAccessIterator, typename OutputIterator >
void LW2004( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, double alpha = PI/2 ) {
    using namespace lw2004;

    // ensure valid alpha
    alpha = CGAL::max( EPSILON, CGAL::min( alpha, PI/2 ) );

    // add points to vector and remove any duplicates or AutoCount() will fail
//    vector<Point> points( pointsBegin, pointsEnd );
//    std::sort( points.begin(), points.end() );
//    auto last = std::unique( points.begin(), points.end() );
//    points.erase( last, points.end() );

    //cout << "Step 1 starts...\n";
    Delaunay T;

    T.insert( pointsBegin, pointsEnd );

    // Add IDs
    size_t i=0;
    for( auto v=T.finite_vertices_begin(); v!=T.finite_vertices_end(); ++v )
        v->info() = i++;

    assert( i == T.number_of_vertices() );

    size_t n = i;

    Vertex_handle v_inf = T.infinite_vertex();

    //cout << "Step 1 is over...\n";
    // TriangulationPrinter tp(T);
    // tp.draw("del");
    //************* Step 2 ****************//
    vector<Vertex_handle> pointID2VertexHandle(n, v_inf);
    for( auto vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit ) {
        //assert( !T.is_infinite(vit) );
        pointID2VertexHandle[ vit->info() ] = vit;
    }

    Heap H;
    vector<handle> handleToHeap(n);
    vector<size_t> piIndexedByV(n), piIndexedByPiU(n);
    vector<unordered_set<size_t>> currentNeighbours(n);

   // size_t maxDegree = 0;
    // Initialize the vector currentNeighbours with appropriate neighbours for every vertex
    for( size_t it = 0; it < n; it++ ) {
        Vertex_circulator N = T.incident_vertices( pointID2VertexHandle.at(it) ),
            done(N);
        //if (vc != 0) {
        do {
            if( !T.is_infinite(N) ) {
                //degree++;
                currentNeighbours.at(it).insert( N->info() );
            }
        } while( ++N != done );
        //}
       // if(degree > maxDegree)
        //    maxDegree = degree;

        size_t degree = currentNeighbours.at(it).size();
        handleToHeap[it] = H.emplace( degree,it );
    }

    //cout << "Maximum degree in the Delaunay Triangulation: " << maxDegree << endl;
    // Use a heap to walk through G_0 to G_{n-1} and set up the Pi for every vertex
    i = n-1; // start at the last valid index
    while(!H.empty()) {
        size_tPair p = H.top();
        H.pop();
        // make sure our math is correct, e.g., degree from heap key matches neighbor container size
        assert( p.first == currentNeighbours.at( p.second ).size() );
        // make sure our assumptions about the graph G_i are correct (see p. 5 in LW2004)
        assert( 0 <= p.first && p.first <= 5 );

        for( size_t neighbour : currentNeighbours.at( p.second ) ) {
            currentNeighbours.at(neighbour).erase(p.second);
            handle h = handleToHeap.at(neighbour);
            size_tPair q = make_pair( currentNeighbours.at( neighbour ).size(), neighbour );
            H.update(h,q);
            H.update(h);
        }
        currentNeighbours.at(p.second).clear();
        piIndexedByV[p.second] = i;
        piIndexedByPiU[i] = p.second;
        --i;
    }

    handleToHeap.clear();
    H.clear();
    currentNeighbours.clear();
    //cout << "Step 2 is over...\n";



    //************* Step 3 ****************//
    // In this step we assume alpha = pi/2 in order to minimize the degree
    size_tPairSet ePrime; // without set duplicate edges could be inserted (use the example down below)
    vector<bool> isProcessed(n, false);
    Delaunay::Vertex_handle u_handle = v_inf;

    // Iterate through vertices by pi ordering
    for( size_t u : piIndexedByPiU ) {
        u_handle = pointID2VertexHandle.at(u);
        assert( !T.is_infinite(u_handle) );
        isProcessed[u] = true;
        //if( printLog ) cout<<"\nu:"<<u<<" ";

        // Get neighbors of u
        Vertex_circulator N = T.incident_vertices( u_handle );
        //if( printLog ) cout<<"N_init:"<<(T.is_infinite(N)?size_t(numeric_limits<size_t>::max):size_t(N->info()))<<" ";
        // find a processed neighbor if it exists or we reach the start again
        while( T.is_infinite(--N) );
        Vertex_circulator done(N); // set done to a vertex that is not infinite
        while( ( T.is_infinite(--N) || !isProcessed.at(N->info()) ) && N!=done );
        done = N; // update N
//        // Rotate N until reaching an infinite vertex (check first) or the original vertex
//        while( !T.is_infinite(--N) && N != done ); //
//        if( T.is_infinite(N) ) --N; // if we stopped on an infinite vertex, move to its successor
        //if( printLog ) cout<<"N_ready:"<<N->info()<<" ";

        // Find and store sector boundaries, start with N
        size_t processedNeighbors = isProcessed.at( N->info() ) ? 1 : 0;
        vector<Vertex_handle> sectorBoundaries{ N->handle() };
        while( --N != sectorBoundaries.front() ) {
            if( ( !T.is_infinite(N) && isProcessed.at( N->info() ) ) ) { // check for v_inf first or isProcessed will be out of range
                sectorBoundaries.push_back( N->handle() );
                ++processedNeighbors;
            }
//            else if( T.is_infinite(N) ) {
//                ++N; // move to predecessor
//                // Add predecessor if not already added.
//                if( N->handle() != sectorBoundaries.back() )
//                    sectorBoundaries.push_back( N->handle() );
//                // Note, if we've reached this branch, we know that
//                // v_inf's successor is sectorBoundaries.front()
//                --N; // move back to N (infinite vertex) to naturally exit the loop
//            }
        }
        assert( processedNeighbors <= 5 );

        // for debugging, print sector boundaries
//        if( printLog ) for ( size_t i=0; i<sectorBoundaries.size(); ++i ) {
//            cout<<"v"<<i<<":"<<sectorBoundaries.at(i)->info()<<" ";
//        }

        // Now, compute the angles of the sectors, the number of cones in each sector,
        // and the actual angles
        vector<double> alphaReal( sectorBoundaries.size() );
        vector< vector<Vertex_handle> > closest( sectorBoundaries.size() );

        for( size_t i=0; i<sectorBoundaries.size(); ++i ) {
            double sectorAngle = get_angle<K>(
                sectorBoundaries.at(i)->point(),
                u_handle->point(),
                sectorBoundaries.at( (i+1)%sectorBoundaries.size() )->point()
            );
            if( sectorAngle < EPSILON ) sectorAngle = 360.0;
            size_t numCones = rint( ceil( sectorAngle / alpha ) );
            assert( numCones > 0 ); // guard against /0
            alphaReal[i] = sectorAngle / numCones;
            closest.at(i).resize( numCones, v_inf );
        }

        Vertex_handle lastN = v_inf;
        //if( isProcessed.at( N->info() ) ) --N; // if N is processed, step
        size_t sector = -1; // The first sector boundary will increment sector, which should be 0 for the first sector

        do { // Loop through neighbors and add appropriate edges
            if( !T.is_infinite(N) ) {
                // If N is the next sectorBoundary, increment sector
                if( N->handle() == sectorBoundaries.at((sector+1)%sectorBoundaries.size()) )
                    ++sector;
                // It is possible for a sectorBoundary to be not processed,
                // in the case of no processed neighbors.
                if( !isProcessed.at( N->info() ) ) {
                    assert( sector < sectorBoundaries.size() );
                    // evaluate possible forward edges
                    double theta = get_angle<K>(
                        sectorBoundaries.at(sector)->point(),
                        u_handle->point(),
                        N->point()
                    );
                    // get angle will return 360 for any angle(vuv) (we want it to be 0 here)
                    size_t cone = size_t( (theta-EPSILON) / alphaReal.at(sector) );
                    if( cone >= closest.at(sector).size() )
                        cone = 0;
                    // Store value until after all neighbors are processed, then add
                    if( T.is_infinite( closest.at(sector).at(cone) )
                      || CGAL::squared_distance( u_handle->point(), N->point() )
                       < CGAL::squared_distance( u_handle->point(), closest.at(sector).at(cone)->point() ) ) {
                            closest.at(sector).at(cone) = N->handle();   // if the saved vertex is infinite or longer than the current one, update
//                            if( printLog ) cout<<"s_closest["<<sector<<"]["<<cone<<"]:"<<N->info()<<" ";
                    }
                    // cross edges
                    if( !T.is_infinite( lastN ) && !isProcessed.at( lastN->info() ) ) {
//                        if( printLog ) cout<<"cross_";
                        createNewEdge( T, pointID2VertexHandle, ePrime, lastN->info(), N->info(), n, false );
                    }
                }
            }
            lastN = N->handle();
        } while( --N != sectorBoundaries.front() );

        // If N and lastN are not processed, add final cross edge
        if( !T.is_infinite(     N ) && !isProcessed.at(     N->info() )
         && !T.is_infinite( lastN ) && !isProcessed.at( lastN->info() ) )
        {
//            if( printLog ) cout<<"cross_";
            createNewEdge( T, pointID2VertexHandle, ePrime, lastN->info(), N->info(), n, false );
        }

        // Add edges in closest
        for( auto segment : closest )
            for( auto v : segment )
                if( !T.is_infinite(v) ) {
//                    if( printLog ) cout<<"forward_";
                    createNewEdge( T, pointID2VertexHandle, ePrime, u, v->info(), n, false );
                }
    }

    // Edge list is only needed for printing. Remove for production.
//    vector< pair<Point,Point> > edgeList;
//    edgeList.reserve( ePrime.size() );

    // Send resultant graph to output iterator
    for( size_tPair e : ePrime ) {
        //edgeList.emplace_back( pointID2VertexHandle.at(e.first)->point(), pointID2VertexHandle.at(e.second)->point() );

        *result = make_pair( pointID2VertexHandle.at(e.first)->point(), pointID2VertexHandle.at(e.second)->point() );
        ++result;
        *result = make_pair( pointID2VertexHandle.at(e.second)->point(), pointID2VertexHandle.at(e.first)->point() );
        ++result;
    }


    //
    //
    // START PRINTER NONSENSE
    //
    //


//    if( printLog ) {
//        GraphPrinter printer(0.007);
//        GraphPrinter::OptionsList options;
//
//        options = {
//            { "color", printer.inactiveEdgeColor },
//            { "line width", to_string(printer.inactiveEdgeWidth) }
//        };
//        printer.drawEdges( T, options );
//
//        options = { // active edge options
//            { "color", printer.activeEdgeColor },
//            { "line width", to_string(printer.activeEdgeWidth) }
//        };
//        printer.drawEdges( edgeList.begin(), edgeList.end(), options );
//
//
//        options = {
//            { "vertex", make_optional( to_string(printer.vertexRadius) ) }, // vertex width
//            { "color", make_optional( printer.backgroundColor ) }, // text color
//            { "fill", make_optional( printer.activeVertexColor ) }, // vertex color
//            { "line width", make_optional( to_string(0) ) } // vertex border (same color as text)
//        };
//        GraphPrinter::OptionsList borderOptions = {
//            { "border", make_optional( to_string(printer.vertexRadius) ) }, // choose shape of vertex
//            { "color", printer.activeEdgeColor }, // additional border color
//            { "line width", to_string(printer.inactiveEdgeWidth) }, // additional border width
//        };
//        printer.drawVerticesWithInfo( T, options, borderOptions );
//
//        printer.print( "lw2004" );
//        cout<<"\n";
//    }




    //
    //
    // END PRINTER NONSENSE
    //
    //

} // function LW2004

} // namespace gsnunf

#endif // GSNUNF_LW2004_H
