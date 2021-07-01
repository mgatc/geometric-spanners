#ifndef GSNUNF_LW2004_H
#define GSNUNF_LW2004_H

#include <algorithm> // min, max
#include <cmath> // ceil
#include <unordered_set> // hashed adjacency list
#include <vector> // vertex containers

#include "DelaunayGraph.h"
#include "GeometricSpannerPrinter.h"
#include "utilities.h"
#include "metrics.h"


namespace gsnunf {

using namespace std;

namespace lw2004 {

using namespace CGAL;

inline size_tPair createEdge( const size_t i, const size_t j )
{
    return make_pair( std::min(i,j), std::max(i,j) );
}

inline void createNewEdge( const Delaunay_triangulation& T,
                           const vector<Vertex_handle>& handles,
                           size_tPairSet &E,
                           const size_t i,
                           const size_t j,
                           const size_t n,
                           bool printLog = false )
{
    assert( std::max(i,j) < n );
    assert( T.is_edge( handles.at(i), handles.at(j) ) );
    //if( printLog ) cout<<"add:("<<i<<","<<j<<") ";
    E.insert( createEdge( i, j ) );
}

} // namespace lw2004

// alpha is set to pi/2
template< typename RandomAccessIterator, typename OutputIterator >
void LW2004( RandomAccessIterator pointsBegin,
             RandomAccessIterator pointsEnd,
             OutputIterator result,
             double alpha = PI/2 )
{
    using namespace lw2004;

    // ensure valid alpha
    alpha = CGAL::max( EPSILON, CGAL::min( alpha, PI/2 ) );

    vector<Point> P(pointsBegin, pointsEnd);
    vector<size_t> index;
    spatialSort<K>(P, index);

    //Step 1: Construct Delaunay triangulation
    Delaunay_triangulation T;

    //N is the number of vertices in the delaunay triangulation.
    size_t n = P.size();
    if(n > SIZE_T_MAX - 1) return;

    //Stores all the vertex handles (CGAL's representation of a vertex, its properties, and data).
    vector<Vertex_handle> handles(n);

    /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
      (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
    Face_handle hint;
    for(size_t entry : index) {
        auto vh = T.insert(P[entry], hint);
        hint = vh->face();
        vh->info() = entry;
        handles[entry] = vh;
    }


    Vertex_handle v_inf = T.infinite_vertex();

    //cout << "Step 1 is over...\n";
    // TriangulationPrinter tp(T);
    // tp.draw("del");
    //************* Step 2 ****************//

    vector<size_t> piIndexedByPiU;
    piIndexedByPiU.reserve(n);
    reverseLowDegreeOrdering(T,back_inserter(piIndexedByPiU));




    //************* Step 3 ****************//
    // In this step we assume alpha = pi/2 in order to minimize the degree
    size_tPairSet ePrime; // without set duplicate edges could be inserted (use the example down below)
    vector<bool> isProcessed(n, false);
    Vertex_handle u_handle = v_inf;

    // Iterate through vertices by pi ordering
    for( size_t u : piIndexedByPiU ) {
        u_handle = handles.at(u);
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
                        createNewEdge( T, handles, ePrime, lastN->info(), N->info(), n, false );
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
            createNewEdge( T, handles, ePrime, lastN->info(), N->info(), n, false );
        }

        // Add edges in closest
        for( auto segment : closest )
            for( auto v : segment )
                if( !T.is_infinite(v) ) {
//                    if( printLog ) cout<<"forward_";
                    createNewEdge( T, handles, ePrime, u, v->info(), n, false );
                }
    }

    // Edge list is only needed for printing. Remove for production.
//    vector< pair<Point,Point> > edgeList;
//    edgeList.reserve( ePrime.size() );

    // Send resultant graph to output iterator
    for( size_tPair e : ePrime ) {
        //edgeList.emplace_back( handles.at(e.first)->point(), handles.at(e.second)->point() );

        *result = e;
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
