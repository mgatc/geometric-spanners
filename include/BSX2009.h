#ifndef GSNUNF_BSX2009_H
#define GSNUNF_BSX2009_H

#include <algorithm> // min, max
#include <cmath> // ceil
#include <unordered_set> // hashed adjacency list
#include <vector> // vertex containers

#include <boost/functional/hash.hpp> // hashing pairs
#include <boost/heap/fibonacci_heap.hpp> // ordering

#include <CGAL/circulator.h> // vertex circulators
#include <CGAL/Delaunay_triangulation_2.h> // DT
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // K
#include <CGAL/Triangulation_vertex_base_with_info_2.h> // DT
#include <CGAL/utils.h> // min, max

#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "metrics.h"
#include "utilities.h"


namespace gsnunf {

using namespace std;

namespace bsx2009 {

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K>      Vb;
typedef CGAL::Triangulation_face_base_2<K>                          Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef CGAL::Aff_transformation_2<K>                               Transformation;
typedef Delaunay::Vertex_handle                                     Vertex_handle;
typedef Delaunay::Vertex_circulator                                 Vertex_circulator;
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

inline bool createNewEdge( const Delaunay& T, const vector<Delaunay::Vertex_handle>& handles, size_tPairSet &E, const size_t i, const size_t j, const size_t n, bool printLog = false ) {
    assert( std::max(i,j) < n );
    assert( T.is_edge( handles.at(i), handles.at(j) ) );
    //if( printLog ) cout<<"add:("<<i<<","<<j<<") ";

    bool inserted = false;
    tie(ignore,inserted) = E.insert( makeNormalizedPair(i,j) );
    return inserted;
}

} // namespace bsx2009

template< typename RandomAccessIterator, typename OutputIterator >
void BSX2009( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, double alpha = 2*PI/3, bool printLog = false ) {
    using namespace bsx2009;

    // ensure valid alpha
    alpha = CGAL::max( EPSILON, CGAL::min( alpha, 2*PI/3 ) );
    size_t numCones = rint( ceil( 2*PI / alpha ) );
    assert( numCones > 0 ); // guard against /0
    double alphaReal = 2*PI / numCones;
    size_t FINAL_DEGREE_BOUND = 14 + numCones;


    vector<Point> P(pointsBegin, pointsEnd);
    vector<size_t> index;
    spatialSort<K>(P, index);

    //Step 1: Construct Delaunay triangulation
    bsx2009::Delaunay T;

    //N is the number of vertices in the delaunay triangulation.
    size_t n = P.size();
    if(n > SIZE_T_MAX - 1) return;

    //Stores all the vertex handles (CGAL's representation of a vertex, its properties, and data).
    vector<bsx2009::Vertex_handle> handles(n);

    /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
      (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
    Delaunay::Face_handle hint;
    for(size_t entry : index) {
        auto vh = T.insert(P[entry], hint);
        hint = vh->face();
        vh->info() = entry;
        handles[entry] = vh;
    }

    Vertex_handle v_inf = T.infinite_vertex();



    //************* Step 2 ****************//
    Heap H;
    vector<handle> handleToHeap(n);
    vector<size_t> piIndexedByV(n), piIndexedByPiU(n);
    vector<unordered_set<size_t>> currentNeighbors(n);

    // Initialize the vector currentNeighbors with appropriate neighbors for every vertex
    for( size_t it=0; it<n; ++it ) {
        Vertex_circulator N = T.incident_vertices( handles.at(it) ),
            done(N);
        do {
            if( !T.is_infinite(N) )
                currentNeighbors.at(it).insert( N->info() );
        } while( ++N != done );

        size_t degree = currentNeighbors.at(it).size();
        handleToHeap[it] = H.emplace( degree, it );
    }

    // Use a heap to walk through G_0 to G_{n-1} and set up the Pi for every vertex
    size_t i = n-1; // start at the last valid index

    while( !H.empty() ) {
        size_tPair p = H.top();
        H.pop();
        // make sure our math is correct, e.g., degree from heap key matches neighbor container size
        assert( p.first == currentNeighbors.at( p.second ).size() );
        assert( 0 <= p.first && p.first <= 5 ); // Lemma 1

        // Erase this vertex from incidence list of neighbors and update the neighbors' key in the heap
        for( size_t neighbor : currentNeighbors.at( p.second ) ) {
            currentNeighbors.at(neighbor).erase(p.second);
            handle h = handleToHeap.at(neighbor);
            size_tPair q = make_pair( currentNeighbors.at( neighbor ).size(), neighbor );
            H.update(h,q);
            H.update(h);
        }
        currentNeighbors.at(p.second).clear();
        piIndexedByV[p.second] = i;
        piIndexedByPiU[i] = p.second;
        --i;
    }

    handleToHeap.clear();
    H.clear();
    currentNeighbors.clear();
    //cout << "Step 2 is over...\n";



    //************* Step 3 ****************//
    size_tPairSet ePrime;
    vector<bool> isProcessed(n, false);
    Vertex_handle u_handle = v_inf;

    // Iterate through vertices by piU ordering
    for( size_t u : piIndexedByPiU ) {
        u_handle = handles.at(u);
        assert( !T.is_infinite(u_handle) );
        Point u_point = u_handle->point();
        isProcessed[u] = true;
        //if( printLog ) cout<<"\nu:"<<u<<" ";

        // Get neighbors of u
        Vertex_circulator N = T.incident_vertices( u_handle );
        //if( printLog ) cout<<"N_init:"<<(T.is_infinite(N)?size_t(numeric_limits<size_t>::max):size_t(N->info()))<<" ";
        if( T.is_infinite(N) ) --N;

        Vertex_circulator
            closest(N),
            done(N);

        // Track degree and number of processed neighbors to ensure correctness
        size_t processedNeighbors = 0;
        size_t degree = 0;
        bool inserted = false;

        do { // Find closest unprocessed neighbor, also count processed neighbors and neighbors in ePrime
            if( !T.is_infinite(N) ) {
                if( gsnunf::contains( ePrime, makeNormalizedPair( u, N->info() ) ) )
                    ++degree;

                if( isProcessed.at( N->info() ) )
                    ++processedNeighbors;
                else if( CGAL::squared_distance(u_point,N->point()) < CGAL::squared_distance(u_point,closest->point()) )
                    closest = N;
            }
        } while( --N != done );

        //if( printLog ) cout<<"degree:"<<degree<<",";
        assert( processedNeighbors <= 5 ); // Lemma 1, proof for Lemma 3
        assert( degree <= 15 ); // Proof for Lemma 3

        // We will add a max of numCones-1 since we are guaranteed to add the closest
        // but cannot add to the two cones touching closest.
        vector<Vertex_handle> closestInCones( numCones-1, v_inf );
        closestInCones.front() = closest; // add closest to "add" list
        while( --N != closest ); // start from the closest vertex

        // Loop through neighbors and consider forward edges
        while( --N != closest ) {
            if( !T.is_infinite(N) && !isProcessed.at(N->info()) ) {
                // evaluate possible forward edges
                double theta = get_angle<K>(
                    closest->point(),
                    u_point,
                    N->point()
                );
                size_t cone = size_t( (theta-EPSILON) / alphaReal );
                // trap neighbors in forbidden cones by putting them in 0 (which is already guaranteed to be closest)
                cone = ( cone < closestInCones.size() ? cone : 0 );

                if( cone > 0 // banish the forbidden cones
                    && ( T.is_infinite( closestInCones.at(cone) )
                        || CGAL::squared_distance(u_point,N->point()) < CGAL::squared_distance(u_point,closestInCones.at(cone)->point()) ) )
                {   // If we made it through all that, it's the current closestInCone!
                    closestInCones[cone] = N->handle();
                    //if( printLog ) cout<<"s_closest["<<cone<<"]:"<<N->info()<<" ";
                }
            }
        }
        // We've found all the closest neighbors in each now,
        // now add edges from each to the current vertex (u)
        for( auto v : closestInCones ) {
            if( !T.is_infinite(v) ) {
                //if( printLog ) cout<<"forward_";
                inserted = createNewEdge( T, handles, ePrime, u, v->info(), n, printLog );
                //degree += size_t(inserted);
                //if( printLog ) cout<<"degree:"<<degree<<",";
            }
        }

        // Loop through neighbors again and add cross edges between
        // consecutive neighbors that are NOT processed (or infinite).
        Vertex_handle lastN = N;
        do {
            --N; // Increment first, then check validity
            if( !( T.is_infinite(N) || isProcessed.at(N->info()) ) ) {
                if( !( T.is_infinite(lastN) || isProcessed.at(lastN->info()) ) ) {
                    // don't add to degree for cross edges, they are not incident on u!
                    //if( printLog ) cout<<"cross_";
                    inserted = createNewEdge( T, handles, ePrime, lastN->info(), N->info(), n, printLog );
                }
            }
            lastN = N->handle();
        } while( N != closest );



//        if( printLog && degree>FINAL_DEGREE_BOUND ) {
//            break; // exit to printer
//        } else if( !printLog && degree>FINAL_DEGREE_BOUND ) {
//            // recurse to print log
//            cout<<"BOUND:"<<FINAL_DEGREE_BOUND<<",deg:"<<degree<<"\n";
//            list<pair<Point,Point> > result;
//            BSX2009( pointsBegin, pointsEnd, back_inserter(result), alpha, true );
//        }
        assert( degree <= FINAL_DEGREE_BOUND ); // Lemma 3

    } // END OF STEP 3 LOOP



    // Edge list is only needed for printing. Remove for production.
//    vector< pair<Point,Point> > edgeList;
//    edgeList.reserve( ePrime.size() );

    // Send resultant graph to output iterator
    for( size_tPair e : ePrime ) {
        // Edge list is only needed for printing. Remove for production.
        //edgeList.emplace_back( handles.at(e.first)->point(), handles.at(e.second)->point() );

        *result = e;
        ++result;
//        *result = make_pair( handles.at(e.second)->point(), handles.at(e.first)->point() );
//        ++result;
    }


    //
    //
    // START PRINTER NONSENSE
    //
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
//        printer.print( "bsx2009" );
//        cout<<"\n";
//    }




    //
    //
    // END PRINTER NONSENSE
    //
    //

} // function BSX2009

} // namespace gsnunf

#endif // GSNUNF_BSX2009_H
