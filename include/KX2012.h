
#ifndef GSNUNF_KX2012_H
#define GSNUNF_KX2012_H

#include <cmath>         // ceil, floor
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

#include <boost/functional/hash.hpp> // size_t pair hash

#include <CGAL/algorithm.h> //
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

//#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "metrics.h"
#include "utilities.h"


namespace gsnunf {

    using namespace std;

    namespace kx2012 {

        typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
        typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K>      Vb;
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

        typedef pair<size_t, size_t>                                         size_tPair;
        typedef boost::hash<size_tPair>                                     size_tPairHash;
        typedef unordered_map<size_tPair, bool, size_tPairHash>               size_tPairMap;

        bool selectEdge(const Delaunay& T, size_tPairMap& E, const Vertex_handle i, const Vertex_handle j) {
            assert(T.is_edge(i, j));
            //if( printLog ) cout<<"add:("<<i->info()<<","<<j->info()<<") ";

            auto existing = E.begin();
            bool inserted = false;
            tie(existing, inserted) = E.try_emplace(makeNormalizedPair(i->info(), j->info()), false);
            if (!inserted) existing->second = true;

            return inserted;
        }

    } // namespace kx2012

    template< typename RandomAccessIterator, typename OutputIterator >
    void KX2012(RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false) {
        using namespace kx2012;

        // ensure k >= 14
        size_t k = 14;
        const double alpha = 2 * PI / k;

        //if(printLog) cout<<"alpha:"<<alpha<<",";

        // Construct Delaunay triangulation

        vector<Point> P(pointsBegin, pointsEnd);
        vector<size_t> index;
        spatialSort<K>(P, index);

        //Step 1: Construct Delaunay triangulation
        kx2012::Delaunay T;

        //N is the number of vertices in the delaunay triangulation.
        size_t n = P.size();
        if (n > SIZE_T_MAX - 1) return;

        //Stores all the vertex handles (CGAL's representation of a vertex, its properties, and data).
        vector<kx2012::Vertex_handle> handles(n);

        /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
          (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
        Delaunay::Face_handle hint;
        for (size_t entry : index) {
            auto vh = T.insert(P[entry], hint);
            hint = vh->face();
            vh->info() = entry;
            handles[entry] = vh;
        }

        kx2012::Vertex_handle v_inf = T.infinite_vertex();
<<<<<<< Updated upstream
        size_tPairMap G_prime; // list of potential edges, value must be true for insertion to result

        // Iterate through vertices in T
=======
        size_tPairMap Degree_Eleven_Graph; // list of potential edges, value must be true for insertion to result

        // Iterate through vertices in T
        
>>>>>>> Stashed changes
        for (auto m = T.finite_vertices_begin(); m != T.finite_vertices_end(); ++m) {
            //if( printLog ) cout<<"\n\nm:"<<m->info()<<" ";

            // Get neighbors of m
            kx2012::Vertex_circulator N = T.incident_vertices(m);
<<<<<<< Updated upstream
            //if( printLog ) cout<<"N_init:"<<(T.is_infinite(N)?size_t(numeric_limits<size_t>::max):size_t(N->info()))<<" ";
            if (T.is_infinite(N)) --N;

            kx2012::Vertex_circulator done(N);
            //if(printLog) cout<<"done:"<<done->info()<<",";

            // closest vertex in each cone
            vector<kx2012::Vertex_handle> closestInCones(k, v_inf);
            // Now, let's put the neighbors found so far into a hashed set for quick lookup
            unordered_set<kx2012::Vertex_handle> selected(k);

            do { // Loop through neighbors and consider forward edges
                if (!T.is_infinite(N)) {
                    // evaluate possible forward edges
                    double theta = get_angle<kx2012::K>(
                        done->point(),
                        m->point(),
                        N->point()
                        );
                    size_t cone = size_t(theta / alpha);
                    //if(printLog) cout<<"N:"<<N->info()<<",theta:"<<theta<<",cone:"<<cone<<",";

                    if (T.is_infinite(closestInCones.at(cone))
                        || distance(m->point(), N->point()) < distance(m->point(), closestInCones.at(cone)->point()))
                    {   // If we made it through all that, it's the current closestInCone!
                        closestInCones[cone] = N;
                        //if( printLog ) cout<<"s_closest["<<cone<<"]:"<<N->info()<<",";
                    }
                }
            } while (--N != done);

            // We've found all the closest neighbors in each cone
            // Put them in selected
            for (auto v : closestInCones)
                if (!T.is_infinite(v))
                    selected.emplace(v);

            // Now we must find every maximal sequence of empty cones
            size_t l = 0; // size of maximal empty sequence
            size_t l_local = 0; // size of current empty sequence
            size_t offset = 0; // offset in case an empty set "wraps" through the end and start of the vector
            size_t startOfSequence = 0; // start of current empty sequence
            unordered_set<size_t> startOfMaximalSequences(k / 2);

            for (size_t i = 0; i < (k + offset); ++i) {
                //if(printLog) cout<<"i:"<<i<<",";
                if (T.is_infinite(closestInCones.at(i % k))) { // empty cone
                    ++l_local;          // increment
                    if (l_local > l) {  // biggest thus far, clear old starts and update l
                        startOfMaximalSequences.clear();
                        l = l_local;
                    }
                    if (l_local >= l)  // place the current start in the list
                        startOfMaximalSequences.emplace(startOfSequence);
                    if (i + 1 == k + offset) {  // if we're about to end but on an empty sequence, keep going
                        ++offset;
                        //if(printLog) cout<<"++offset,";
                    }
                    //if(printLog) cout<<"l_local:"<<l_local<<",";
                }
                else {                    // filled cone
                 //if(printLog) cout<<"filledby:"<<closestInCones.at(i%k)->info()<<",";
                    l_local = 0;                 // reset l_local
                    startOfSequence = (i + 1) % k; // set the start of sequence to the next i
                }
            }
            //        if( printLog ) {
            //            cout<<"l:"<<l<<",";
            //            cout<<"num_seq:"<<startOfMaximalSequences.size()<<",";
            //        }
                    // loop through starts of maximal sequences and add edges for them
            for (auto start : startOfMaximalSequences) {
                double startAngle = start * alpha;
                //            if( printLog ) cout << "startOfMaximalSeq:"<< start<<",";
                //            if( printLog ) cout << "startAngle:"<< startAngle<<",";

                while (--N != done); // point N to reference point
                // point N to first neighbor past the empty sequence, if it exists
                while (--N != done
                    && (T.is_infinite(N)
                        || get_angle<K>(done->point(), m->point(), N->point()) < startAngle)
                    );

                kx2012::Vertex_circulator afterSequence(N),
                    beforeSequence(N);
                while (T.is_infinite(++beforeSequence)); // move once CCW, if it's infinite move again

    //            if( printLog ) cout << "beforeSeq:"<< beforeSequence->info() <<",";
    //            if( printLog ) cout << "afterSeq:"<< afterSequence->info() <<",";

                //return;
                if (l > 1) {
                    // select the first ceil(l/2) unselected edges CCW
                    size_t remainingToAdd = size_t(rint(ceil(l / 2.0)));
                    //if( printLog ) cout << "CCWadds:"<< remainingToAdd<<",";

                    while (remainingToAdd > 0 && ++N != afterSequence) {
                        if (!T.is_infinite(N) && !contains(selected, N)) {
                            selected.emplace(N);
                            --remainingToAdd;
                            //if( printLog ) cout << "compensating:"<< N->info() <<",";
                        }
                    }

                    // select the first floor(l/2) unselected edges CW
                    remainingToAdd = size_t(rint(floor(l / 2.0)));
                    N = beforeSequence; // move N to the neighbor before the sequence
                    //if( printLog ) cout << "CWadds:"<< remainingToAdd<<",";

                    while (remainingToAdd > 0 && --N != beforeSequence) {
                        if (!T.is_infinite(N) && !contains(selected, N)) {
                            selected.emplace(N);
                            --remainingToAdd;
                            //if( printLog ) cout << "compensating:"<< N->info() <<",";
                        }
                    }
                }
                else if (l == 1) {
                    //if( printLog ) cout << "addOne,";
                    kpx2010::Vertex_handle singleSelection = v_inf;

                    // consider the first CW and CCW edges (held in beforeSequence and afterSequence)
                    // if one is selected already, add the other
                    if (contains(selected, beforeSequence) ^ contains(selected, afterSequence)) {
                        singleSelection = contains(selected, beforeSequence) ? afterSequence : beforeSequence;
                    }
                    // otherwise, add the shorter
                    else if (!(contains(selected, beforeSequence) || contains(selected, afterSequence))) {
                        singleSelection =
                            CGAL::squared_distance(beforeSequence->point(), m->point())
                            < CGAL::squared_distance(afterSequence->point(), m->point())
                            ? beforeSequence : afterSequence;
                    }
                    // if we found one to select, select it!
                    if (!T.is_infinite(singleSelection)) {
                        selected.emplace(singleSelection);
                        //if( printLog ) cout << "compensating:"<< singleSelection->info() <<",";
                    }
                }
            }

            bool inserted = false;
            // now add edges from each to the current vertex (u)
            for (auto v : selected) {
                if (!T.is_infinite(v)) {
                    //if( printLog ) cout<<"forward_";
                    inserted = selectEdge(T, G_prime, m, v);
                }
            }

        }

        // Done. Send edges from G_prime with value == true (selected by both endpoints) to output.

        // Edge list is only needed for printing. Remove for production.
    //    vector< pair<Point,Point> > edgeList;
    //    edgeList.reserve( G_prime.size() );

        // Send resultant graph to output iterator
        for (auto e : G_prime) {
            if (e.second) { // e.second holds the bool value of whether both vertices of an edge selected the edge
                // Edge list is only needed for printing. Remove for production.
                //edgeList.emplace_back( handles.at(e.first.first)->point(), handles.at(e.first.second)->point() );
=======
            
            if(T.is_infinite(N)) ++N; //Verify N is not starting at infinity
            kx2012::Vertex_circulator done(N); //Artificial end to circulator

            //if( printLog ) cout<<"N_init:"<<(T.is_infinite(N)?size_t(numeric_limits<size_t>::max):size_t(N->info()))<<" ";
            
           
            
            //if(printLog) cout<<"done:"<<done->info()<<",";
            vector<kx2012::Vertex_handle> angle_Set; //Store the 3 verticies in a vector for each angle set
            while(angle_Set.size() < 3)
            {
                if(!T.is_infinite(N)) angle_Set.emplace_back(N);
                N++;
            }
            
            //N is set to 
            Point p = m->point();
            size_t currentPointIndex = m->info();
            double angle_Sum =0;
            set<kx2012::Vertex_handle> WideVertices;
            
            do{ if (T.is_infinite(N)) {N++;}
                
                angle_Sum = get_angle(angle_Set[2]->point(),p,angle_Set[0]->point());

                if(angle_Sum > FOUR_PI_OVER_FIVE)
                {
                    WideVertices.insert(angle_Set.begin(),angle_Set.end());
                    if(printLog)
                    {
                        cout << "Points: {" << angle_Set[2]->info() << ","<< m->info()<< "," << angle_Set[0]->info() << "} make angle: " << angle_Sum << endl;

                    }
                }
                
                angle_Set[0] = angle_Set[1];
                angle_Set[1] = angle_Set[2];
                angle_Set[2] = N++;
                
            }while(angle_Set[0] != done);
            
            for(auto v : WideVertices){
                selectEdge(T,Degree_Eleven_Graph,m,v);
            }
            if (T.is_infinite(N)){N++;}
            done = N;
            Point Cone_Calculation_Point(p.x(),p.y()+1);
            double conalAngle = 0;
            vector<kx2012::Vertex_handle> closestVertexInCone(10,v_inf);
            vector<double> closestPointDistanceInCone(10,INF);
            do{ 
                if (!T.is_infinite(N)&& !contains(WideVertices,N)){

                    conalAngle = get_angle(N->point(),p,Cone_Calculation_Point);   
                    size_t currentCone = conalAngle/PI_OVER_FIVE;
                    double currentDistance = distance(p,N->point());

                    if(currentDistance < closestPointDistanceInCone[currentCone])
                    {
                        closestVertexInCone[currentCone] = N;
                        closestPointDistanceInCone[currentCone] = currentDistance;
                    }

                    
                }
            }while(++N != done);

            for(auto v : closestVertexInCone)
            {
                if(v != v_inf)
                {
                    selectEdge(T,Degree_Eleven_Graph,m,v);
                }
            }
       
        }
       
        // Done. Send edges from G_prime with value == true (selected by both endpoints) to output.

        // Edge list is only needed for printing. Remove for production.
       vector< pair<size_t,size_t> > edgeList;
       edgeList.reserve( Degree_Eleven_Graph.size() );

        // Send resultant graph to output iterator
        for (auto e : Degree_Eleven_Graph) {
            if (e.second) { // e.second holds the bool value of whether both vertices of an edge selected the edge
                // Edge list is only needed for printing. Remove for production.
                edgeList.push_back(e.first);
>>>>>>> Stashed changes

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

<<<<<<< Updated upstream
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
=======
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
           printer.drawEdges( edgeList.begin(), edgeList.end(),P, options );
    
    
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
    
           printer.print( "KX2012" );
           cout<<"\n";
       }
>>>>>>> Stashed changes

        //
        //
        // END PRINTER NONSENSE
        //
        //

    } // function KX2012

} // namespace gsnunf

#endif // GSNUNF_KX2012_H
