//Needs optimizing currently testing.
#ifndef GSNUNF_BHS2017_H
#define GSNUNF_BHS2017_H

//Base libraries.
#include <cmath>         // ceil, floor, isinf
#include <limits>
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

//Boost library
#include <boost/functional/hash.hpp> // size_t pair hash

//CGAL library
#include <CGAL/algorithm.h> //
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Line_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

//Project library
//#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "metrics.h"
#include "utilities.h"


namespace gsnunf {

            using namespace std;

            namespace bhs2017 {

            //CGAL objects
            typedef CGAL::Exact_predicates_inexact_constructions_kernel                         K;
            typedef CGAL::Triangulation_vertex_base_with_info_2<size_t,K>                       Vb;
            typedef CGAL::Triangulation_face_base_2<K>                                          Fb;
            typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                                 Tds;
            typedef CGAL::Delaunay_triangulation_2<K,Tds>                                       Delaunay;
            typedef CGAL::Aff_transformation_2<K>                                               Transformation;
            typedef Delaunay::Vertex_handle                                                     Vertex_handle;
            typedef Delaunay::Vertex_circulator                                                 Vertex_circulator;
            typedef CGAL::Vector_2<K>                                                           Vector_2;
            typedef CGAL::Line_2<K>                                                             Line;
            typedef Delaunay::Point                                                             Point;
            typedef Delaunay::Finite_vertices_iterator                                          Finite_vertices_iterator;
            typedef Delaunay::Finite_edges_iterator                                             Finite_edges_iterator;

            //Project objects
            typedef pair<size_t,size_t>                                                         size_tPair;
            typedef boost::hash<size_tPair>                                                     size_tPairHash;
            typedef unordered_map<size_tPair,bool,size_tPairHash>                               size_tPairMap;
            typedef unordered_map<pair<size_t,size_t>,double,pointPairHash,edgeEquality>        edgeBisectorMap;
            typedef unordered_map<pair<size_t,size_t>,size_t,pointConeHash,pointConeEquality>   pointConeMap;

            const double tan30 = tan(PI / 6);
            //const double cot30 = 1 / tan30;

            const vector<double> bisectorSlopes{ 0, tan30, -1*tan30, 0, tan30, -1*tan30 };
            //const vector<double> orthBisectorSlopes{ INF, -1*cot30, cot30, INF, -1*cot30, cot30 };

            //Finds the cone of p containing vertex q, for this algorithm all vertices have 6 cones (0-5) with an angle of (PI/3).
            inline size_t getSingleCone(const size_t p, const size_t q, const vector<Vertex_handle> &h){
                const double alpha = PI/3;
                const Point refPoint( h.at(p)->point().x() - tan30, h[p]->point().y() + 1 );
                //Point refPoint(h[p]->point().x(), h[p] ->point().y() + 1);

                double theta = get_angle<bhs2017::K>(refPoint, h[p]->point(), h.at(q)->point());

                size_t cone = (theta / alpha);

                return cone;
            }

            // Compute max of getCone(p,q) and (getCone(q,p)+3)%6
            inline size_t getCone( const size_t p, const size_t q, const vector<Vertex_handle> &h ) {
                return max(
                    getSingleCone(p,q,h),
                    (getSingleCone(q,p,h)+3)%6
                );
            }

            //Finds the bisector length of a given edge.
            inline K::FT bisectorLength(const double alpha, const pair<size_t,size_t> &e, const vector<Vertex_handle> &h) {

                size_t cone = getCone(e.first, e.second, h);

                double xCord = h.at(e.first)->point().x();
                double yCord = h[e.first]->point().y() + 1;

                assert(cone<6);
                assert(e.first<h.size());

                if(cone % 3 != 0){
                    xCord = (1/bisectorSlopes.at(cone)) + h[e.first]->point().x();
                }

                Point bisectorPoint(xCord, yCord);

                Line bisectorLine(h[e.first]->point(), bisectorPoint);

                Point intersectionPoint = bisectorLine.projection(h[e.second]->point());

                double bisectorLen = distance(h[e.first]->point(), intersectionPoint);

                return bisectorLen;
            }

            /*
              Step 3: Add incident edges consists of 2 sub steps.
              (3.1) Starts with the empty set E_A.
              (3.2) For each edge in L (sorted in non-decreasing order) Let i be the cone of p containing q if E_A has no edges with endpoint p in the
                    neighborhood of p in cone i and E_A has no edges with endpoint q in the neighborhood of q in cone i+3 then add edge (p,q) to E_A.
            */
            inline void addIncident(vector<pair<size_t,size_t>> &E_A, pointConeMap &AL_E_A, const double alpha,
                const vector<Vertex_handle> &h, const vector<pair<size_t,size_t>> &l){

                //Loops through the entire set L.
                for(auto e : l ){

                    //Separates the edge (p,q) into the vertices p and q.
                    const size_t& p = e.first;
                    const size_t& q = e.second;

                    //Computes the cone of p containing q.
                    size_t p_cone = getCone(p, q, h),
                           q_cone = getCone(q, p, h);
                    if (q_cone != (p_cone + 3) % 6) cout<<"\np_cone:"<<p_cone<<" q_cone:"<<q_cone<<" (p_cone+3)%6:"<<((p_cone+3)%6)<<endl;
                    assert( q_cone == (p_cone + 3) % 6 );

                    /*Evaluates the emptiness of a cone, given a vertex in the set E_A.
                      If a cone is empty then the set E_A does not contain an edge with the
                      given endpoint in the cone calculated above, and the status will be
                      set to true.*/
                    bool p_cone_empty = AL_E_A.find(make_pair(p, p_cone)) == AL_E_A.end();
                    bool q_cone_empty = AL_E_A.find(make_pair(q, q_cone)) == AL_E_A.end();

                    /*Checks that both cone neighborhood are empty, if these are both empty
                    then the condition for step 3 is met and (p,q) is added to E_A. (3.2)*/
                    if( p_cone_empty && q_cone_empty ){
                        E_A.push_back(e);

                        //Adds (p,q) to an adjacency list for future calculation.
                        AL_E_A.emplace(make_pair(p, p_cone), q);
                        AL_E_A.emplace(make_pair(q, q_cone), p);
                    }
                }
            }
            inline void canonicalNeighborhood( vector<size_t>& canNeighbors,
                                               const size_t& p,
                                               const size_t& r,
                                               const size_t cone,
                                               const Delaunay &dt,
                                               const vector<Vertex_handle> &h,
                                               const edgeBisectorMap &b,
                                               bool printLog = false ) {

                pair<size_t,size_t> e = make_pair(p, r);
                /*Vertex circulator oriented to r to find fist and last end vertex. Once r is the circulator is oriented to the first vertex in the cone,
                  that is in the canonical neighborhood. Once found all neighbors are added in clockwise order. For a vertex to be in the canonical
                  neighborhood it must have a bisector length greater than or equal to that of (p,r)*/
                auto N_p = dt.incident_vertices(h[p]);

                while(++N_p != h[r]);

                while(!dt.is_infinite(++N_p)
                      && getCone(p, N_p->info(), h) == cone
                      && (b.at(make_pair(p, N_p->info())) > b.at(e)
                          || abs(b.at(make_pair(p, N_p->info())) - b.at(e)) < EPSILON
                          ));


                while(!dt.is_infinite(--N_p)
                      && getCone(p, N_p->info(), h) == cone
                      && (b.at(make_pair(p, N_p->info())) > b.at(e)
                          || abs(b.at(make_pair(p,N_p->info())) - b.at(e)) < EPSILON
                          //|| N_p->info() == r
                         )){
                    canNeighbors.push_back(N_p->info());
                }

                // print neighborhood of p
                // traverse the neighbors until we reach infinite or the cone isn't 0
//                while(!dt.is_infinite(--N_p) || getCone(p, N_p->info(), h) == 0);
//
//                auto done = N_p;
//                size_t c = 6; // 6 is out of range for cones
//                while( !dt.is_infinite(--N_p) || N_p != done ) {
//                    // if this is a new cone, print cone label
//                    if( getCone( p, N_p->info(), h ) != c ) {
//                        // print canonical neighbors of previous cone
//
//                        c = getCone( p, N_p->info(), h );
//                        cout<<";  C"<<c<<": ";
//                    }
//                    cout<<N_p->info()<<",";
//                }
//                vector<size_t> canNeighborsPrint;
//                canonicalNeighborhood( canNeighborsPrint, p, N_p->info(), cone, dt, h, b );



            }

            /*
              Step 4: Add canonical edges consists of 4 sub steps.
              (4.1) Let r be an element of the canonical neighborhood of p in the cone 0 of p.
              (4.2) Add inner edges if total neigborhood edges is 3 or more.
              (4.3) If r is an end vertex and there is more than one edge in the neighborhood add the edge with endpoint r.
              (4.4) Consider the first and last edge in the canoncial neighborhood. 3 criteria to add.
                (4.4 a) If the edges are in cone 1 or 5 with respect to a and z add.
                (4.4 b) If the edges are in cone 2 or 4 with respect to a and z and cone for has no edge with an end edge point in E_A add.
                (4.4 c) Checks if end edges have a end point a or z in E_A and an edge different from one made with vertex b or y in cone 2 or 4 woth respect
                        to a and z if found the edge (b,c) or (w,y) is added.
            */
            inline void addCanonical(vector<pair<size_t,size_t>> &E_CAN, const size_t p, const size_t r, const double alpha,
                const Delaunay &dt, const vector<Vertex_handle> &h, const edgeBisectorMap &b, pointConeMap& AL_e_a){

                //Creates an edge (p,r)
                pair<size_t,size_t> e = make_pair(p, r);

                //Computes the cone of p containing r.
                size_t p_cone  = getCone(p, r, h);

                //Set of the canonical neighborhood of p in the cone of p containing r. (This cone will be considered as cone 0)
                vector<size_t> canNeighbors;

                canonicalNeighborhood( canNeighbors, p, r, p_cone, dt, h, b );
                assert(canNeighbors.size()>0);
                //Number of edges in the neighborhood.
                int canEdges = canNeighbors.size() - 1;

                //Add inner edges if total neigborhood edges is 3 or more. (4.2)
                for(int i = 1; i < int(canNeighbors.size()) - 2; i++){
                    E_CAN.emplace_back(canNeighbors.at(i), canNeighbors.at(i + 1));
                    assert(dt.is_edge( h.at(canNeighbors.at(i)), h.at(canNeighbors.at(i + 1))));
                }
                //If r is an end vertex and there is more than one edge in the neighborhood add the edge with endpoint r. (4.3)
                if(canNeighbors.front() == r && canEdges > 1){
                    E_CAN.emplace_back(r, canNeighbors.at(1));
                    assert(dt.is_edge( h.at(r), h.at(canNeighbors.at(1))));
                }
                if(canNeighbors.back() == r && canEdges > 1){
                    E_CAN.emplace_back(r, canNeighbors.at(canEdges - 1));
                    assert(dt.is_edge( h.at(r), h.at(canNeighbors.at(canEdges - 1))));
                }

                //First and last edges in the canonical neighborhood are considered and added by 3 criteria. (4.4)
                //Must be at least 1 edge.
                if(canEdges > 0){
                    vector<int> cone(6);
                    for( size_t i=0; i<6; ++i ) {
                        cone[i] = (p_cone+i)%6;
                    }

                    const pair<size_t,size_t> canFirst = make_pair(canNeighbors.front(), canNeighbors.at(1));
                    const pair<size_t,size_t> canLast = make_pair(canNeighbors.at(canEdges - 1), canNeighbors.back());

                    //Iterator to end of map to check if an edge exists.
                    const auto blank = AL_e_a.end();

                    //If the edges are in cone 1 or 5 with respect to a and z add. (4.4 a)
                    if(getCone(canFirst.first, canFirst.second, h) == cone[1]){
                        E_CAN.push_back(canFirst);
                        assert(dt.is_edge( h.at(canFirst.first), h.at(canFirst.second)));
                    }
                    if(getCone(canLast.second, canLast.first, h) == cone[5]){
                        E_CAN.push_back(canLast);
                        assert(dt.is_edge( h.at(canLast.first), h.at(canLast.second)));
                    }
                    const auto u = AL_e_a.find(make_pair(canFirst.first, cone[2]));
                    const auto v = AL_e_a.find(make_pair(canLast.second, cone[4]));
                    //If the edges are in cone 2 or 4 with respect to a and z and cone for has no edge with an end edge point in E_A add. (4.4 b)
                    if( u == blank && getCone(canFirst.first, canFirst.second, h) == cone[2] ){
                        E_CAN.push_back(canFirst);
                        assert(dt.is_edge( h.at(canFirst.first), h.at(canFirst.second)));
                    }
                    if( v == blank && getCone(canLast.second, canLast.first, h) == cone[4] ){
                        E_CAN.push_back(canLast);
                        assert(dt.is_edge( h.at(canLast.first), h.at(canLast.second)));
                    }

                    /*Checks if end edges have an end point a or z in E_A and an edge different from one made with vertex b or y in cone 2 or 4 woth respect
                      to a and z if found the edge (b,c) or (w,y) is added. (4.4 c)*/
                    // (a,b)
                    if( getCone(canFirst.first, canFirst.second, h) == cone[2]
                     && u != blank
                     && u->second != canFirst.second ) {
                        vector<size_t> zCanNeighbors;
                        canonicalNeighborhood( zCanNeighbors, canFirst.first, u->second, cone[2], dt, h, b );
//                        auto maybe_p = find( zCanNeighbors.begin(), zCanNeighbors.end(), p );
//                        if( maybe_p != zCanNeighbors.end() )
//                            zCanNeighbors.erase(maybe_p);
                        auto y = find( zCanNeighbors.begin(), zCanNeighbors.end(), canFirst.second );
                        auto w = y;
                        w += int(y == zCanNeighbors.begin());
                        w -= int(y == zCanNeighbors.end()-1);
//                        if(y==zCanNeighbors.end() || y==zCanNeighbors.end()||y==w) {
//                            cout<<"\np:"<<p<<" w:"<<(*w)<<" y:"<<canFirst.second<<" z:"<<canFirst.first;
//                            cout<<"\nzCanNeighbors1:";
//                            for( auto v : zCanNeighbors ) cout<<v<<" ";
//                            cout<<endl;
//                        }
                        assert( y != zCanNeighbors.end() );
                        assert( y != w );
                        assert(dt.is_edge( h.at(*w), h.at(*y)));
                        E_CAN.emplace_back(*w, *y);
                    }
                    // (y,z)
                    if( getCone(canLast.second, canLast.first, h) == cone[4]
                     && v != blank
                     && v->second != canLast.first ) {
                        vector<size_t> zCanNeighbors;
                        canonicalNeighborhood( zCanNeighbors, canLast.second, v->second, cone[4], dt, h, b );
//                        auto maybe_p = find( zCanNeighbors.begin(), zCanNeighbors.end(), p );
//                        if( maybe_p != zCanNeighbors.end() )
//                            zCanNeighbors.erase(maybe_p);
                        auto y = find( zCanNeighbors.begin(), zCanNeighbors.end(), canLast.first );
                        auto w = y;
                        w += int(y == zCanNeighbors.begin());
                        w -= int(y == zCanNeighbors.end()-1);
//                        if(y==zCanNeighbors.end() || y==zCanNeighbors.end()||y==w) {
//                            cout<<"\np:"<<p<<" w:"<<(*w)<<" y:"<<canLast.first<<" z:"<<canLast.second;
//                            cout<<"\nzCanNeighbors2:";
//                            for( auto v : zCanNeighbors ) cout<<v<<" ";
//                            cout<<endl;
//                        }
                        assert( y != zCanNeighbors.end() );
                        assert( y != w );
                        assert(dt.is_edge( h.at(*w), h.at(*y)));
                        E_CAN.emplace_back(*w, *y);
                    }
                }
            }
        } // namespace BHS2017

    //Main algorithm.
    template<typename RandomAccessIterator, typename OutputIterator>
    void BHS2017(RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false) {

        using namespace bhs2017;

        //Angle of the cones. Results in 6 cones for a given vertex.
        const double alpha = PI / 3;

        //Step 1: Construct Delaunay triangulation
        bhs2017::Delaunay DT(pointsBegin, pointsEnd);

        //N is the number of vertices in the delaunay triangulation.
        size_t n = DT.number_of_vertices();
        if(n > SIZE_T_MAX - 1) return;

        //Stores all the vertex handles (CGAL's representation of a vertex, its properties, and data).
        vector<bhs2017::Vertex_handle> handles(n);

        /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
          (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
        size_t i=0;
        for(auto v = DT.finite_vertices_begin(); v != DT.finite_vertices_end(); ++v) {
            v->info() = i;
            handles[i] = v;
            ++i;
        }

        //Put edges in a vector.
        vector<pair<size_t,size_t>> L;

        for(auto e = DT.finite_edges_begin(); e != DT.finite_edges_end(); ++e) {
            L.emplace_back(e->first->vertex((e->second + 1) % 3)->info(),
                            e->first->vertex((e->second + 2) % 3)->info());
        }

        //Creates a map of edges as keys to its respective bisector length as the value. (Edges are not directional 1-2 is equivilent to 2-1)
        edgeBisectorMap B(L.size());

        for(auto e : L){
            B.emplace(e, bisectorLength(alpha, e, handles));
        }

        //Step 2: Edges in the set L are sorted by their bisector length in non-decreasing order.
        sort(L.begin(), L.end(), [&](const auto& lhs, const auto& rhs) {
            return B.at(lhs) < B.at(rhs);
        });

        //Creates a set which will contain all edges returned by addIncident.
        vector<pair<size_t,size_t>> E_A;

        /*Creates an adjacency list where the inner lists are of size 6 representing the cones. The value stored in a particular inner index
          is the vertex that creates an edge with the outer vertex in the given cone. (i.e. If AL_E_A[10][4] = 5 in cone 4 of vertex 10 there
          exists an edge (10,5).)*/
        pointConeMap AL_E_A;

        //Step 3
        addIncident(E_A, AL_E_A, alpha, handles, L);

        //Add canonical E_CAN
        vector<pair<size_t,size_t>> E_CAN;

        //Step 4
        for(auto  e : E_A ){

            addCanonical(E_CAN, e.first, e.second, alpha, DT, handles, B, AL_E_A);
            addCanonical(E_CAN, e.second, e.first, alpha, DT, handles, B, AL_E_A);

        }

        for(auto e:E_A){
            cout<<e.first<<" "<<e.second<<"\n";
        }
        cout<<endl<<endl;
        for(auto e:E_CAN){
            cout<<e.first<<" "<<e.second<<"\n";
        }
        cout<<endl;

        //Union of sets E_A and E_CAN for final edge set removes duplicates.
        E_A.insert(E_A.end(), E_CAN.begin(), E_CAN.end());
        sort(E_A.begin(), E_A.end(), [](const auto& l, const auto& r) {
            return (min(l.first, l.second) < min(r.first, r.second)
                || (min(l.first, l.second) == min(r.first, r.second) && max(l.first, l.second) < max(r.first, r.second)));
        } );
        E_A.erase(unique(E_A.begin(), E_A.end(), [](const auto& l, const auto& r) {
            return (l.first == r.first && l.second == r.second)
                || (l.first == r.second && l.second == r.first);
        }), E_A.end());

        // Edge list is only needed for printing. Remove for production.
        vector<pair<Point,Point>> edgeList;

        // Send resultant graph to output iterator
            for(auto e : E_A) {
                // Edge list is only needed for printing. Remove for production.
                edgeList.emplace_back(handles.at(e.first)->point(), handles.at(e.second)->point());

                *result = make_pair(handles.at(e.first)->point(), handles.at(e.second)->point());
                ++result;
                *result = make_pair(handles.at(e.second)->point(), handles.at(e.first)->point());
                ++result;
            }

            // START PRINTER NONSENSE
            if(printLog) {
                GraphPrinter printer(1);
                GraphPrinter::OptionsList options;

                options = {
                    {"color", printer.inactiveEdgeColor},
                    {"line width", to_string(printer.inactiveEdgeWidth)}
                };
                printer.drawEdges(DT, options);

                options = { // active edge options
                    {"color", printer.activeEdgeColor},
                    {"line width", to_string(printer.activeEdgeWidth)}
                };
                printer.drawEdges(edgeList.begin(), edgeList.end(), options);


                options = {
                    {"vertex", make_optional(to_string(printer.vertexRadius))}, // vertex width
                    {"color", make_optional(printer.backgroundColor)}, // text color
                    {"fill", make_optional(printer.activeVertexColor)}, // vertex color
                    {"line width", make_optional(to_string(0))} // vertex border (same color as text)
                };
                GraphPrinter::OptionsList borderOptions = {
                    {"border", make_optional(to_string(printer.vertexRadius))}, // choose shape of vertex
                    {"color", printer.activeEdgeColor}, // additional border color
                    {"line width", to_string(printer.inactiveEdgeWidth)}, // additional border width
                };
                printer.drawVerticesWithInfo(DT, options, borderOptions);

                printer.print("BHS2017");
                cout << "\n";
            }
            // END PRINTER NONSENSE

    } // function BHS2017

} // namespace gsnunf

#endif // GSNUNF_BHS2017_H
