
#ifndef GSNUNF_KPX2010_H
#define GSNUNF_KPX2010_H

#include <cmath>         // ceil, floor
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

#include <CGAL/algorithm.h>

//#include "GeometricSpannerPrinter.h"
#include "DelaunayGraph.h"
#include "metrics.h"
#include "utilities.h"


namespace gsnunf {

using namespace std;

namespace kpx2010 {

bool selectEdge( const Delaunay_triangulation& T, size_tPairMap &E, const Vertex_handle i, const Vertex_handle j ) {
    assert( T.is_edge( i, j ) );
    //if( printLog ) cout<<"add:("<<i->info()<<","<<j->info()<<") ";

    auto existing = E.begin();
    bool inserted = false;
    tie(existing,inserted) = E.try_emplace( makeNormalizedPair( i->info(), j->info() ), false );
    if(!inserted) existing->second = true;

    return inserted;
}

} // namespace kpx2010

template< typename RandomAccessIterator, typename OutputIterator >
void KPX2010( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, size_t k, bool printLog = false ) {
    using namespace kpx2010;

    // ensure k >= 14
    k = std::max( k, size_t(14) );
    const number_t alpha = 2*PI / k;

    //if(printLog) cout<<"alpha:"<<alpha<<",";

    // Construct Delaunay triangulation

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
    size_tPairMap G_prime; // list of potential edges, value must be true for insertion to result

    // Iterate through vertices in T
    for( auto m=T.finite_vertices_begin(); m!=T.finite_vertices_end(); ++m ) {
        //if( printLog ) cout<<"\n\nm:"<<m->info()<<" ";

        // Get neighbors of m
        Vertex_circulator N = T.incident_vertices(m);
        //if( printLog ) cout<<"N_init:"<<(T.is_infinite(N)?size_t(numeric_limits<size_t>::max):size_t(N->info()))<<" ";
        if( T.is_infinite(N) ) --N;

        Vertex_circulator done(N);
        //if(printLog) cout<<"done:"<<done->info()<<",";

        // closest vertex in each cone
        vector<Vertex_handle> closestInCones( k, v_inf );
        // Now, let's put the neighbors found so far into a hashed set for quick lookup
        unordered_set<Vertex_handle> selected(k);

        do { // Loop through neighbors and consider forward edges
            if( !T.is_infinite(N) ) {
                // evaluate possible forward edges
                number_t theta = get_angle(
                    done->point(),
                    m->point(),
                    N->point()
                );
                auto cone = size_t( theta / alpha );
                //if(printLog) cout<<"N:"<<N->info()<<",theta:"<<theta<<",cone:"<<cone<<",";

                if( T.is_infinite( closestInCones.at(cone) )
                    || distance(m->point(),N->point()) < distance(m->point(),closestInCones.at(cone)->point()) )
                {   // If we made it through all that, it's the current closestInCone!
                    closestInCones[cone] = N;
                    //if( printLog ) cout<<"s_closest["<<cone<<"]:"<<N->info()<<",";
                }
            }
        } while( --N != done );

        // We've found all the closest neighbors in each cone
        // Put them in selected
        for( auto v : closestInCones )
            if( !T.is_infinite(v) )
                selected.emplace(v);

        // Now we must find every maximal sequence of empty cones
        size_t l = 0; // size of maximal empty sequence
        size_t l_local = 0; // size of current empty sequence
        size_t offset = 0; // offset in case an empty set "wraps" through the end and start of the vector
        size_t startOfSequence = 0; // start of current empty sequence
        unordered_set<size_t> startOfMaximalSequences(k/2);

        for( size_t i=0; i<(k+offset); ++i ) {
            //if(printLog) cout<<"i:"<<i<<",";
            if( T.is_infinite( closestInCones.at(i%k) ) ) { // empty cone
                ++l_local;          // increment
                if( l_local > l ) {  // biggest thus far, clear old starts and update l
                    startOfMaximalSequences.clear();
                    l = l_local;
                }
                if( l_local >= l )  // place the current start in the list
                    startOfMaximalSequences.emplace( startOfSequence );
                if( i+1 == k+offset ) {  // if we're about to end but on an empty sequence, keep going
                    ++offset;
                    //if(printLog) cout<<"++offset,";
                }
                //if(printLog) cout<<"l_local:"<<l_local<<",";
            } else {                    // filled cone
                //if(printLog) cout<<"filledby:"<<closestInCones.at(i%k)->info()<<",";
                l_local = 0;                 // reset l_local
                startOfSequence = (i+1) % k; // set the start of sequence to the next i
            }
        }
//        if( printLog ) {
//            cout<<"l:"<<l<<",";
//            cout<<"num_seq:"<<startOfMaximalSequences.size()<<",";
//        }
        // loop through starts of maximal sequences and add edges for them
        for( auto start : startOfMaximalSequences ) {
            number_t startAngle = start*alpha;
//            if( printLog ) cout << "startOfMaximalSeq:"<< start<<",";
//            if( printLog ) cout << "startAngle:"<< startAngle<<",";

            while( --N != done ); // point N to reference point
            // point N to first neighbor past the empty sequence, if it exists
            while( --N != done
              && ( T.is_infinite(N)
                || get_angle( done->point(), m->point(), N->point() ) < startAngle )
            );

            Vertex_circulator afterSequence(N),
                              beforeSequence(N);
            while( T.is_infinite(++beforeSequence) ); // move once CCW, if it's infinite move again

//            if( printLog ) cout << "beforeSeq:"<< beforeSequence->info() <<",";
//            if( printLog ) cout << "afterSeq:"<< afterSequence->info() <<",";

            //return;
            if( l > 1 ) {
                // select the first ceil(l/2) unselected edges CCW
                auto remainingToAdd = size_t(rint(ceil(l/2.0)));
                //if( printLog ) cout << "CCWadds:"<< remainingToAdd<<",";

                while( remainingToAdd > 0 && ++N != afterSequence ) {
                    if( !T.is_infinite(N) && !contains( selected, N ) ) {
                        selected.emplace(N);
                        --remainingToAdd;
                        //if( printLog ) cout << "compensating:"<< N->info() <<",";
                    }
                }

                // select the first floor(l/2) unselected edges CW
                remainingToAdd = size_t(rint(floor(l/2.0)));
                N = beforeSequence; // move N to the neighbor before the sequence
                //if( printLog ) cout << "CWadds:"<< remainingToAdd<<",";

                while( remainingToAdd > 0 && --N != beforeSequence ) {
                    if( !T.is_infinite(N) && !contains( selected, N ) ) {
                        selected.emplace(N);
                        --remainingToAdd;
                        //if( printLog ) cout << "compensating:"<< N->info() <<",";
                    }
                }
            } else if( l == 1 ) {
                //if( printLog ) cout << "addOne,";
                Vertex_handle singleSelection = v_inf;

                // consider the first CW and CCW edges (held in beforeSequence and afterSequence)
                // if one is selected already, add the other
                if( contains(selected,beforeSequence) ^ contains(selected,afterSequence) ) {
                    singleSelection = contains(selected,beforeSequence) ? afterSequence : beforeSequence;
                }
                // otherwise, add the shorter
                else if( !( contains(selected,beforeSequence) || contains(selected,afterSequence) ) ) {
                    singleSelection =
                        CGAL::squared_distance( beforeSequence->point(), m->point() )
                        < CGAL::squared_distance( afterSequence->point(), m->point() )
                            ? beforeSequence : afterSequence;
                }
                // if we found one to select, select it!
                if( !T.is_infinite(singleSelection) ) {
                    selected.emplace( singleSelection );
                    //if( printLog ) cout << "compensating:"<< singleSelection->info() <<",";
                }
            }
        }

        bool inserted = false;
        // now add edges from each to the current vertex (u)
        for( const auto& v : selected ) {
            if( !T.is_infinite(v) ) {
                //if( printLog ) cout<<"forward_";
                inserted = selectEdge( T, G_prime, m, v );
            }
        }

    }

    // Done. Send edges from G_prime with value == true (selected by both endpoints) to output.

    // Edge list is only needed for printing. Remove for production.
//    vector< pair<Point,Point> > edgeList;
//    edgeList.reserve( G_prime.size() );

    // Send resultant graph to output iterator
    for( auto e : G_prime ) {
        if( e.second ) { // e.second holds the bool value of whether both vertices of an edge selected the edge
            // Edge list is only needed for printing. Remove for production.
            //edgeList.emplace_back( handles.at(e.first.first)->point(), handles.at(e.first.second)->point() );

            *result = e.first;
            ++result;
//            *result = reverse_pair(e.first);
//            ++result;
        }
    }


    //
    //
    // START PRINTER NONSENSE
    //
    //

//    if( printLog ) {
//        GraphPrinter printer(1);
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

} // function KPX2010

} // namespace gsnunf

#endif // GSNUNF_KPX2010_H
