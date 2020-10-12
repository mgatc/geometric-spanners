
#ifndef GSNUNF_KPX2010_H
#define GSNUNF_KPX2010_H

#include <cmath>
#include <float.h>
#include <forward_list>
#include <fstream>
#include <limits>
#include <list>
#include <queue>
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include <CGAL/Aff_transformation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>


#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "metrics.h"
#include "utilities.h"


namespace gsnunf {

using namespace std;

namespace kpx2010 {

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

inline bool createNewEdge( const Delaunay& T, size_tPairSet &E, const Vertex_handle i, const Vertex_handle j, const size_t n, bool printLog = false ) {
    assert( T.is_edge( i, j ) );
    if( printLog ) cout<<"add:("<<i->info()<<","<<j->info()<<") ";

    bool inserted = false;
    tie(ignore,inserted) = E.insert( makeNormalizedPair( i->info(), j->info() ) );
    return inserted;
}

} // namespace kpx2010

template< typename RandomAccessIterator, typename OutputIterator >
void KPX2010( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, size_t k, bool printLog = false ) {
    using namespace kpx2010;


    // ensure k >= 14
    k = std::max( k, size_t(14) );
    double alpha = 2*PI / k;

    // Construct Delaunay triangulation
    Delaunay T( pointsBegin, pointsEnd );
    size_t n = T.number_of_vertices();

    vector<kpx2010::Vertex_handle> handles(n);

    // Add IDs
    size_t i=0;
    for( auto v=T.finite_vertices_begin(); v!=T.finite_vertices_end(); ++v ) {
        v->info() = i;
        handles[i] = v;
        ++i;
    }

    Delaunay::Vertex_handle v_inf = T.infinite_vertex();


    //GraphAlgoTV tv;
    vector<string> vertexInfo;
    vertexInfo.reserve(n);
//    getVertexInfo( T, back_inserter(vertexInfo) );
//    tv.registerTriangulation( T, vertexInfo );
//    tv.play();

    //cout<<n<<",";


    //
    //
    // START PRINTER NONSENSE
    //
    //

                GraphPrinter printer(1);
                GraphPrinter::OptionsList options;

                options = {
                    { "color", printer.inactiveEdgeColor },
                    { "line width", to_string(printer.inactiveEdgeWidth) }
                };
                printer.drawEdges( T, options );

                options = {
                    { "vertex", make_optional( to_string(printer.vertexRadius) ) }, // vertex width
                    { "color", make_optional( printer.backgroundColor ) }, // text color
                    { "fill", make_optional( printer.activeVertexColor ) }, // vertex color
                    { "line width", make_optional( to_string(0) ) } // vertex border (same color as text)
                };
                GraphPrinter::OptionsList borderOptions = {
                    { "border", make_optional( to_string(printer.vertexRadius) ) }, // choose shape of vertex
                    { "color", printer.activeEdgeColor }, // additional border color
                    { "line width", to_string(printer.inactiveEdgeWidth) }, // additional border width
                };
                printer.drawVerticesWithInfo( T, options, borderOptions );

                options = { // active edge options
                    { "color", printer.activeEdgeColor },
                    { "line width", to_string(printer.activeEdgeWidth) }
                };
                printer.print( "kpx2010" );





    //
    //
    // END PRINTER NONSENSE
    //
    //














    //************* Step 3 ****************//
    size_tPairSet ePrime;

    // Iterate through vertices in T
    for( auto m=T.finite_vertices_begin(); m!=T.finite_vertices_end(); ++m ) {
        if( printLog ) cout<<"\nm:"<<m->info()<<" ";


//        { // EDGE PRINTING NONSENSE
//            vector< pair<Point,Point> > edgeList;
//            edgeList.reserve( ePrime.size() );
//            // Send resultant graph to output iterator
//            for( size_tPair e : ePrime )
//                edgeList.emplace_back( pointID2VertexHandle.at(e.first)->point(), pointID2VertexHandle.at(e.first)->point() );
//
//            printer.drawEdges( edgeList.begin(), edgeList.end(), options );
//            printer.print( "lw2004" );
//        }



        // Get neighbors of m
        Delaunay::Vertex_circulator N = T.incident_vertices(m);
        if( printLog ) cout<<"N_init:"<<(T.is_infinite(N)?size_t(numeric_limits<size_t>::max):size_t(N->info()))<<" ";
        if( T.is_infinite(N) ) --N;

        Delaunay::Vertex_circulator done(N);

        vector<Delaunay::Vertex_handle> closestInCones( k, v_inf );

        do { // Loop through neighbors and consider forward edges
            if( !T.is_infinite(N) ) {
                // evaluate possible forward edges
                double theta = get_angle<K>(
                    done->point(),
                    m->point(),
                    N->point()
                );
                size_t cone = size_t( (theta-EPSILON) / alpha );


                if( T.is_infinite( closestInCones.at(cone) )
                    || d(m->point(),N->point()) < d(m->point(),closestInCones.at(cone)->point()) )
                {   // If we made it through all that, it's the current closestInCone!
                    closestInCones[cone] = N->handle();
                    if( printLog ) cout<<"s_closest["<<cone<<"]:"<<N->info()<<" ";
                }
            }
        } while( --N != done );
        // We've found all the closest neighbors in each cone

        // Now, let's put the neighbors found so far into a hashed set for quick lookup
        unordered_set<kpx2010::Vertex_handle> selected(k);
        for( auto v : closestInCones )
            if( !T.is_infinite(v) )
                selected.emplace(v);

        // Now we must find every maximal sequence of empty cones
        size_t l = 0; // size of maximal empty sequence
        size_t l_local = 0; // size of current empty sequence
        size_t offset = 0; // offset in case an empty set "wraps" through the end and start of the vector
        size_t startOfSequence = 0; // start of current empty sequence
        unordered_set<size_t> startOfMaximalSequences(k/2);

        for( size_t i=0; i<k+offset; ++i ) {
            if( T.is_infinite( closestInCones.at(i%k) ) ) { // empty
                ++l_local;          // increment
                if( l_local > l ) {  // biggest thus far, clear old starts and update l
                    startOfMaximalSequences.clear();
                    l = l_local;
                }
                if( l_local >= l )  // place the current start in the list
                    startOfMaximalSequences.emplace( startOfSequence );
                if( i+1 == k+offset )   // if we're about to end but on an empty sequence, keep going
                    ++offset;
                if(printLog) cout<<"l_local:"<<l_local<<",";
            } else {                    // filled
                l_local = 0;                 // reset l_local
                startOfSequence = (i+1) % k; // set the start of sequence to the next i
            }
        }
        if( printLog ) {
            cout<<"l:"<<l<<",";
            cout<<"num_seq:"<<startOfMaximalSequences.size()<<",";
        }
        // loop through starts of maximal sequences and add edges for them
        for( auto start : startOfMaximalSequences ) {
            //while( --N != handles.at(start) ); // point N at the start
            if( l > 1 ) {
                // select the first ceil(l/2) unselected edges CCW
                double startAngle = start*alpha;
                size_t remainingToAdd = size_t(rint(ceil(l/2)));

                while( --N != done ); // point N to reference point
                // point N to first neighbor past the empty sequence
                while( T.is_infinite(--N) || get_angle<K>( done->point(), m->point(), N->point() ) < startAngle );
                kpx2010::Vertex_circulator afterSequence(N),
                                          beforeSequence(N);
                ++beforeSequence;

                while( remainingToAdd > 0 && ++N != afterSequence ) {
                    if( !T.is_infinite(N) && !contains( selected, N ) ) {
                        selected.emplace(N);
                        --remainingToAdd;
                    }
                }
                // select the first floor(l/2) unselected edges CW
                double endAngle = startAngle + l*alpha;
                remainingToAdd = size_t(rint(floor(l/2)));

                N = beforeSequence; // move N to the neighbor before the sequence

                while( remainingToAdd > 0 && --N != beforeSequence ) {
                    if( !T.is_infinite(N) && !contains( selected, N ) ) {
                        selected.emplace(N);
                        --remainingToAdd;
                    }
                }
            } else if( l == 1 ) {
                // consider the first CW and CCW edges
                // if one is selected already, add the other
                // otherwise, add the shorter
            }
        }

        bool inserted = false;
        // now add edges from each to the current vertex (u)
        for( auto v : closestInCones ) {
            if( !T.is_infinite(v) ) {
                if( printLog ) cout<<"forward_";
                inserted = createNewEdge( T, ePrime, m, v, n, printLog );
                //degree += size_t(inserted);
                //if( printLog ) cout<<"degree:"<<degree<<",";
            }
        }

    } // END OF STEP 3 LOOP



    // Edge list is only needed for printing. Remove for production.
    vector< pair<Point,Point> > edgeList;
    edgeList.reserve( ePrime.size() );

    // Send resultant graph to output iterator
    for( size_tPair e : ePrime ) {
        // Edge list is only needed for printing. Remove for production.
        edgeList.emplace_back( handles.at(e.first)->point(), handles.at(e.second)->point() );

        *result = make_pair( handles.at(e.first)->point(), handles.at(e.second)->point() );
        ++result;
        *result = make_pair( handles.at(e.second)->point(), handles.at(e.first)->point() );
        ++result;
    }


    //
    //
    // START PRINTER NONSENSE
    //
    //


    if( printLog ) {
        GraphPrinter printer(1);
        GraphPrinter::OptionsList options;

        options = {
            { "color", printer.inactiveEdgeColor },
            { "line width", to_string(printer.inactiveEdgeWidth) }
        };
        printer.drawEdges( T, options );

        options = { // active edge options
            { "color", printer.activeEdgeColor },
            { "line width", to_string(printer.activeEdgeWidth) }
        };
        printer.drawEdges( edgeList.begin(), edgeList.end(), options );


        options = {
            { "vertex", make_optional( to_string(printer.vertexRadius) ) }, // vertex width
            { "color", make_optional( printer.backgroundColor ) }, // text color
            { "fill", make_optional( printer.activeVertexColor ) }, // vertex color
            { "line width", make_optional( to_string(0) ) } // vertex border (same color as text)
        };
        GraphPrinter::OptionsList borderOptions = {
            { "border", make_optional( to_string(printer.vertexRadius) ) }, // choose shape of vertex
            { "color", printer.activeEdgeColor }, // additional border color
            { "line width", to_string(printer.inactiveEdgeWidth) }, // additional border width
        };
        printer.drawVerticesWithInfo( T, options, borderOptions );

        printer.print( "bsx2009" );
        cout<<"\n";
    }




    //
    //
    // END PRINTER NONSENSE
    //
    //

} // function KPX2010

} // namespace gsnunf

#endif // GSNUNF_KPX2010_H
