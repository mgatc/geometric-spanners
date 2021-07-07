#ifndef GSNUNF_BCC2012_H
#define GSNUNF_BCC2012_H

#include <bitset>
#include <cmath>         // ceil, floor, isinf
#include <functional>
#include <limits>
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles


#include <CGAL/algorithm.h> //
#include <CGAL/circulator.h>

#include "DelaunayGraph.h"
#include "GeometricSpannerPrinter.h"
#include "metrics.h"
#include "utilities.h"


namespace unf_planespanners {

    using namespace std;

    namespace bcc2012 {

        enum Q_primePosition {
            between_j_i = 0, between_i_k = 1, not_set = 2
        };

        struct WedgeParameters {
            index_t p;
            index_t q;
            cone_t cone;
        };

        inline cone_t getCone(const vector<VertexHandle> &handles,
                              const vector<index_t> &closest,
                              const size_t p,
                              const size_t q,
                              const number_t alpha) {
            return cone_t(getAngle(handles[closest[p]]->point(),
                                   handles[p]->point(),
                                   handles[q]->point()) / alpha);
        }

        inline cone_t getPreviousCone(const cone_t cone, const cone_t numCones) {
            return (cone + numCones - 1) % numCones;
        }

        template<cone_t DEGREE, cone_t NUM_CONES = DEGREE + 1>
        inline bool vertexAgreesOnEdge(const vector<VertexHandle> &handles,
                                       vector<index_t> &closest,
                                       const vector<bitset<NUM_CONES>> &filled,
                                       const index_t p,
                                       const index_t q,
                                       cone_t &cone,
                                       cone_t &conePrev,
                                       bool &qOnBoundary) {

            if (closest.at(p) == SIZE_T_MAX) { // First, make sure the closest vertex is set
                closest.at(p) = q;
                qOnBoundary = true; // the only vertex that falls on a boundary should be the closest
            }
            const number_t ALPHA = 2 * PI / NUM_CONES;
            cone = getCone(handles, closest, p, q, ALPHA);
            conePrev = getPreviousCone(cone, NUM_CONES);

            bool pGivenConeFilled = filled.at(p)[cone],
                    pPrevConeFilled = filled.at(p)[conePrev];

            return (!qOnBoundary && !pGivenConeFilled)
                   || (qOnBoundary && (!pGivenConeFilled || !pPrevConeFilled));
        }

        template<cone_t DEGREE, cone_t NUM_CONES = DEGREE + 1>
        inline void updateVertexConeStatus(vector<bitset<NUM_CONES>> &filled,
                                           vector<WedgeParameters> &wedge,
                                           const index_t p,
                                           const index_t q,
                                           const bool qOnBoundary,
                                           const cone_t cone,
                                           const cone_t conePrev) {
            if (qOnBoundary) {
                if (!filled.at(p)[conePrev])
                    wedge.push_back({p, q, conePrev});
                filled.at(p)[conePrev] = true;
            }
            if (!filled.at(p)[cone])
                wedge.push_back({p, q, cone});
            filled.at(p)[cone] = true;
        }


/*
 *  Performs algorithm "wedge" with the exception of line 1. The
 *  for loop must be accomplished outside of the function.
 *  Degree 6 and 7 wedge algorithms are implemented as templates.
 *  This is the primary template, but will never be used. Instead,
 *  the specialized templates are used, defined below.
 */
        template<cone_t DEGREE, cone_t NUM_CONES = DEGREE + 1>
        inline void wedge(const DelaunayTriangulation &DT,
                          const vector<VertexHandle> &handles,
                          const vector<index_t> &closest,
                          const WedgeParameters &params,
                          vector<Edge> &addToE_star,
                          const VertexCirculator &q_i) {
            assert(DEGREE == 6 || DEGREE == 7);
        }

        template<>
        inline void wedge<7>(const DelaunayTriangulation &DT,
                             const vector<VertexHandle> &handles,
                             const vector<index_t> &closest,
                             const WedgeParameters &params,
                             vector<Edge> &addToE_star,
                             const VertexCirculator &q_i) {

            const number_t ALPHA = 2 * PI / 8;

            const index_t &p = params.p;
            const cone_t &cone = params.cone;

            // q_m[i] holds the circulator for q_{m-i}
            vector<VertexCirculator> q_m(3);

            // Setup function objects for wedge
            vector<std::function<VertexCirculator(VertexCirculator &)>> step;
            step.emplace_back([](VertexCirculator &c) { return ++c; });
            step.emplace_back([](VertexCirculator &c) { return --c; });

            for (size_t i = 0; i < step.size(); i++) {
                fill(q_m.begin(), q_m.end(), q_i);

                // set and increment q_m and q_{m-1} so the sequence will be correct at the start of the loop
                for (size_t j = 0; j < q_m.size() - 1; ++j)
                    step[i](q_m[j]);

                index_t a = p;
                //     b = q_i
                index_t c = p;

                if (i == 0)
                    c = q_m[1]->info();
                else
                    a = q_m[1]->info();

                // Get the first and last vertex in cone_p, called q_j and q_k
                // Rotate q_j CCW until we leave the cone
                while (!DT.is_infinite(q_m[0])
                       && getCone(handles, closest, p, q_m[0]->info(), ALPHA) == cone
                       && !DT.is_infinite(step[i](q_m[0]))
                       && getCone(handles, closest, p, q_m[0]->info(), ALPHA) == cone) {
                    if (q_m[2] != q_i
                        || (q_m[2] == q_i &&
                            getAngle(handles.at(a)->point(), q_i->point(), handles.at(c)->point()) > PI_OVER_TWO)) {
                        addToE_star.emplace_back(q_m[1]->info(), q_m[2]->info());
                    }
                    step[i](q_m[2]);
                    step[i](q_m[1]);
                }
            }
        }

        template<>
        inline void wedge<6>(const DelaunayTriangulation &DT,
                             const vector<VertexHandle> &handles,
                             const vector<index_t> &closest,
                             const WedgeParameters &params,
                             vector<Edge> &addToE_star,
                             const VertexCirculator &q_i) {

            const number_t ALPHA = 2 * PI / 7;

            const index_t &p = params.p;
            const index_t &q = params.q;
            const cone_t &cone = params.cone;

            // Line 2: Build Q
            auto N_p = DT.incident_vertices(handles[p]);

            while (++N_p != handles[q]); // point to q aka q_i

            while (!DT.is_infinite(++N_p) // move CCW until we leave the cone
                   && (getCone(handles, closest, p, N_p->info(), ALPHA) == cone
                       || N_p->info() == q));

            vector<index_t> Q; // the ordered neighbors in the current cone

            while (!DT.is_infinite(--N_p) // move CW until we leave the cone, adding each to Q
                   && (getCone(handles, closest, p, N_p->info(), ALPHA) == cone
                       || N_p->info() == q)) { // handles the case where q is on a boundary
                Q.push_back(N_p->info());
            }
            assert(!Q.empty());

            // Line 3: Build Q'
            unordered_set<index_t> Q_prime; // select elements from Q (line 3)

            size_t j = 0, // the first point in the ordered list of neighbors is q_j
            i = 0, // the opposite end of the edge
            k = Q.size() - 1; // the last point in the ordered list is q_k

            // find the index of q_i in Q
            for (size_t n = 0; n < Q.size(); ++n)
                if (Q[n] == q)
                    i = n;


            Q_primePosition Q_primePos = not_set;

            for (size_t n = j + 1; n < k; ++n) {
                if (n != i // q_i
                    && getAngle(handles[Q[n + 1]]->point(),
                                handles[Q[n]]->point(),
                                handles[Q[n - 1]]->point()) < SIX_PI_OVER_SEVEN) {
                    Q_prime.insert(Q[n]);
                    Q_primePos = static_cast<Q_primePosition>(n > i); // if we have passed i, n>i will be true == 1 == between_i_k
                }
//        if(printLog)cout<<getAngle<K>( handles.at(Q.at(n+1))->point(),
//                             handles.at(Q.at(n))->point(),
//                             handles.at(Q.at(n-1))->point() )<<" <>? "<<SIX_PI_OVER_SEVEN<<"\n";
            }
//    if(printLog) {
//        cout<<"p:"<<p<<" q:"<<q<<"\n";
//        cout<<"Q:";
//        for( auto v : Q ) cout<<v<<" ";
//
//        cout<<"\n";
//        cout<<"Q':";
//        for( auto v : Q_prime ) cout<<v<<" ";
//
//        cout<<"\n";
//
//        cout<<"j:"<<j<<" i:"<<i<<" k:"<<k<<endl;
//        cout<<"j+1:"<<int(j+1)<<" i-2:"<<int(i-2)<<" k:"<<k<<endl;
//    }

            // Line 4: Add select edges
            if (i > 1)
                for (size_t n = j + 1; n < i - 2; ++n)
                    if (!contains(Q_prime, Q[n]) && !contains(Q_prime, Q[n + 1]))
                        addToE_star.emplace_back(Q[n], Q[n + 1]);

//    if(printLog)cout<<"add to E_star:";
//    if(printLog)for( auto e : addToE_star ) cout<<e.first<<" "<<e.second<<" - ";

            if (k > 1)
                for (size_t n = i + 1; n < k - 2; ++n)
                    if (!contains(Q_prime, Q[n]) && !contains(Q_prime, Q[n + 1]))
                        addToE_star.emplace_back(Q[n], Q[n + 1]);

//    if(printLog)cout<<"add to E_star:";
//    if(printLog)for( auto e : addToE_star ) cout<<e.first<<" "<<e.second<<" - ";

            size_t f = i, // will hold the index in Q of the first point in Q_prime
            a = i;

            // Line 5:
            switch (Q_primePos) {
                case not_set:
                    // skip
                    break;
                case between_i_k:
                    // Line 6-7
                    if (i != j
                        && i - 1 != j
                        && getAngle(handles[p]->point(),
                                    handles[Q[i]]->point(),
                                    handles[Q[i - 1]]->point()) > FOUR_PI_OVER_SEVEN)
                        addToE_star.emplace_back(Q[i], Q[i - 1]);

                    //if(printLog&&i!=j&&i-1!=j)cout<<getAngle<K>(handles.at(p)->point(), handles.at(Q.at(i))->point(), handles.at(Q.at(i-1))->point())<<"\n";

                    while (++f < Q.size() && !contains(Q_prime, Q[f]));

                    a = f - int(f == Q.size()) * Q.size(); // the first point after f not in Q_prime, avoid vector overflow
                    while (++a < Q.size() && contains(Q_prime, Q[a]));

                    //if(printLog)cout<<"f:"<<f<<" a:"<<a<<"\n";

                    if (f == i + 1) {
                        if (a != k && getAngle(handles[Q[i + 1]]->point(),
                                               handles[Q[i]]->point(),
                                               handles[p]->point()) < FOUR_PI_OVER_SEVEN) {
                            addToE_star.emplace_back(Q[f], Q[a]);
                        }
//            if(printLog)cout<<getAngle<K>( handles.at(Q.at(i+1))->point(),
//                                         handles.at(Q.at(i))->point(),
//                                         handles.at(p)->point() )<<"\n";
                        if (i != j && i != k && f + 1 != k && getAngle(handles[Q[i + 1]]->point(),
                                                                       handles[Q[i]]->point(),
                                                                       handles[p]->point()) > FOUR_PI_OVER_SEVEN) {
                            addToE_star.emplace_back(Q[i], Q[f + 1]);
                        }
                    } else {
                        size_t l = Q.size(); // the last point in Q_prime
                        while (--l > 0 && !contains(Q_prime, Q[l]));
                        if (l > 0) {
                            size_t b = l; // set to max n that is less than l in Q but not Q_prime
                            while (--b > 0 && contains(Q_prime, Q[b]));

//                if(printLog)cout<<"l:"<<l<<" b:"<<b<<endl;
                            if (l == k - 1) {
                                addToE_star.emplace_back(Q[l], Q[b]);
                            } else {
                                addToE_star.emplace_back(Q[b], Q[l + 1]);
                                if (contains(Q_prime, Q[l - 1]))
                                    addToE_star.emplace_back(Q[l], Q[l - 1]);
                            }
                        }
                    }
                    break;
                case between_j_i:
                    // Line 6-7
                    if (i != k
                        && i + 1 != k
                        && getAngle(handles[Q[i + 1]]->point(),
                                    handles[Q[i]]->point(),
                                    handles[p]->point()) > FOUR_PI_OVER_SEVEN) {
                        addToE_star.emplace_back(Q[i], Q[i + 1]);
                    }

//        if(printLog)if(i!=k&&i+1!=k)cout<<getAngle<K>(handles.at(Q.at(i+1))->point(), handles.at(Q.at(i))->point(), handles.at(p)->point())<<"\n";

                    while (--f > 0 && !contains(Q_prime, Q[f]));

                    a = f + int(f == 0) *
                            Q.size(); // the first point after f not in Q_prime, protect against unsigned underflow
                    while (--a > 0 && contains(Q_prime, Q[a]));

//        if(printLog)cout<<"f:"<<f<<" a:"<<a<<"\n";

                    if (f == i - 1) {
//            if(printLog)cout<<getAngle<K>( handles.at(p)->point(),
//                                         handles.at(Q.at(i))->point(),
//                                         handles.at(Q.at(i-1))->point() )<<"\n";
                        if (a != j && getAngle(handles[p]->point(),
                                               handles[Q[i]]->point(),
                                               handles[Q[i - 1]]->point()) < FOUR_PI_OVER_SEVEN) {
                            addToE_star.emplace_back(Q[f], Q[a]);
                        }
                        if (i != j && i != k && f - 1 != j && getAngle(handles[p]->point(),
                                                                       handles[Q[i]]->point(),
                                                                       handles[Q[i - 1]]->point()) >
                                                              FOUR_PI_OVER_SEVEN) {
                            addToE_star.emplace_back(Q[i], Q[f - 1]);
                        }
                    } else {
                        size_t l = 0; // the last point in Q_prime
                        while (++l < Q.size() && !contains(Q_prime, Q[l]));
                        if (l > 0) {
                            size_t b = l; // set to max n that is less than l in Q but not Q_prime
                            while (++b < Q.size() && contains(Q_prime, Q[b]));

//                cout<<"l:"<<l<<" b:"<<b<<endl;
                            if (l == j + 1) {
                                addToE_star.emplace_back(Q[l], Q[b]);
                            } else {
                                addToE_star.emplace_back(Q[b], Q[l - 1]);
                                if (contains(Q_prime, Q[l + 1]))
                                    addToE_star.emplace_back(Q[l], Q[l + 1]);
                            }
                        }
                    }
                    break;
            }

//    if(printLog)cout<<"add to E_star:";
//    if(printLog)for( auto e : addToE_star ) cout<<e.first<<" "<<e.second<<" - ";
//
//    if(printLog)cout<<"\n\n";
        }
//#pragma GCC diagnostic pop


    } // namespace bcc2012

    template<size_t DEGREE = 7, size_t NUM_CONES = DEGREE + 1,
            typename RandomAccessIterator, typename OutputIterator>
    void BCC2012(RandomAccessIterator pointsBegin,
                 RandomAccessIterator pointsEnd,
                 OutputIterator result,
                 bool printLog = false) {
        using namespace bcc2012;

        assert(DEGREE == 7 || DEGREE == 6);

//    if(printLog) cout<<"\nnumCones:"<<NUM_CONES<<"\n";
//    if(printLog) cout<<"ALPHA:"<<ALPHA<<"\n";

        // Construct Delaunay triangulation
        vector<Point> P(pointsBegin, pointsEnd);
        vector<index_t> index;
        spatialSort<K>(P, index);

        //Step 1: Construct Delaunay triangulation
        DelaunayTriangulation DT;

        //N is the number of vertices in the delaunay triangulation.
        index_t n = P.size();
        if (n > SIZE_T_MAX - 1) return;

        //Stores all the vertex handles (CGAL's representation of a vertex, its properties, and data).
        vector<VertexHandle> handles(n);

        /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
          (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
        FaceHandle hint;
        for (size_t entry : index) {
            auto vh = DT.insert(P[entry], hint);
            hint = vh->face();
            vh->info() = entry;
            handles[entry] = vh;
        }

//    if(printLog) cout<<"n:"<<n<<"\n";


        // Put edges in a vector, then sort on weight
        vector<Edge> L;

        for (auto e = DT.finite_edges_begin(); e != DT.finite_edges_end(); ++e) {
            L.emplace_back(e->first->vertex((e->second + 1) % 3)->info(),
                           e->first->vertex((e->second + 2) % 3)->info());
        }
        sort(L.begin(), L.end(), [&](const auto &lhs, const auto &rhs) {
            return getEdgeLength(lhs, P) < getEdgeLength(rhs, P);
        });

        vector<index_t> closest(n, SIZE_T_MAX); // id of closest vertex for orienting cones
        vector<bitset<NUM_CONES>> filled(n); // status of each cone
        vector<Edge> E; // output edge list
        vector<Edge> E_star; // edges added from "Wedge"

        for (auto pq : L) {
            size_t p = pq.first,
                   q = pq.second;

//        if(printLog) cout<<"p-q:"<<p<<" - "<<q<<"\n";
//        if(printLog) cout<<"p-q:"<<handles.at(p)->point()<<" - "<<handles.at(q)->point()<<"\n";
//        if(printLog) cout<<"  p_filled:"<<filled.at(p)<<"\n";
//        if(printLog) cout<<"  q_filled:"<<filled.at(q)<<"\n";

            // If either p or q's cone is filled, don't even bother
            if (filled.at(p).count() == NUM_CONES || filled.at(q).count() == NUM_CONES)
                continue;

            // Politely ask p if it wants an edge to q
            cone_t cone_p = 0,
                   cone_pPrev = 0;
            bool qOnBoundary = false;
            bool pAbides = vertexAgreesOnEdge<DEGREE>(handles, closest, filled, p, q,
                                                      cone_p, cone_pPrev, qOnBoundary);

//        if(printLog) cout<<"  cone_p:"<<cone_p<<"\n";
//        if(printLog && qOnBoundary) cout<<"  qOnBoundary\n";
//        if(printLog && pAbides ) cout<<"  p abides!\n";


            // Politely ask q if it wants an edge to p
            cone_t cone_q = 0,
                   cone_qPrev = 0;
            bool pOnBoundary = false;
            bool qAbides = vertexAgreesOnEdge<DEGREE>(handles, closest, filled, q, p,
                                                      cone_q, cone_qPrev, pOnBoundary);

//        if(printLog) cout<<"  cone_q:"<<cone_q<<"\n";
//        if(printLog && pOnBoundary) cout<<"  pOnBoundary\n";
//        if(printLog && qAbides ) cout<<"  q abides!\n";

            // Only continue if p and q both consent to add the edge
            if (pAbides && qAbides) {
                E.emplace_back(p, q); // Place the edge

                // Wedge on each cone of pq and qp
                // There will be at least one for each, but there could
                // be two cones for one or both pq and qp if the edgeat(0)
                // falls on the boundary of a cone and the cone is not already filled
                vector<WedgeParameters> W; // holds the parameters for each call to wedge

                // Bookkeeping for p
                updateVertexConeStatus<DEGREE>(filled, W, p, q, qOnBoundary, cone_p, cone_pPrev);

                // Bookkeeping for q
                updateVertexConeStatus<DEGREE>(filled, W, q, p, pOnBoundary, cone_q, cone_qPrev);

                vector<Edge> addToE_star;

                // Wedge on p, q
                for (auto params : W) {
                    // find q
                    auto q_z = DT.incident_vertices(handles.at(params.p));
                    while (++q_z != handles.at(params.q)); // point to q
                    const auto q_i(q_z);

                    wedge<DEGREE>(DT, handles, closest, params, addToE_star, q_i);
                }
                E_star.insert(E_star.end(), addToE_star.begin(), addToE_star.end());
            }
        }

        // Combine E and E_star, remove duplicates
        E.insert(E.end(), E_star.begin(), E_star.end());
        sort(E.begin(), E.end());
        E.erase(unique(E.begin(), E.end(), [](const auto &l, const auto &r) {
            return (l.first == r.first && l.second == r.second)
                   || (l.first == r.second && l.second == r.first);
        }), E.end());


        // Edge list is only needed for printing. Remove for production.
//    vector< pair<Point,Point> > edgeList;
//    edgeList.reserve( E.size() );

        // Send resultant graph to output iterator
        for (auto e : E) {
            // Edge list is only needed for printing. Remove for production.
            //edgeList.emplace_back( handles.at(e.first)->point(), handles.at(e.second)->point() );

            *result = e;
            ++result;
        }

        //
        //
        // START PRINTER NONSENSE
        //
        //

//    if( printLog && n <= 200 ) {
//        GraphPrinter printer(1);
//        GraphPrinter::OptionsList options;
//
//        options = {
//            { "color", printer.inactiveEdgeColor },
//            { "line width", to_string(printer.inactiveEdgeWidth) }
//        };
//        printer.drawEdges( DT, options );
//
//        options = { // active edge options
//            { "color", printer.activeEdgeColor },
//            { "line width", to_string(printer.activeEdgeWidth) }
//        };
//        printer.drawEdges( E.begin(), E.end(), P, options );
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
//        printer.drawVerticesWithInfo( DT, options, borderOptions );
//
//        printer.print( "bcc2012" );
//        cout<<"\n";
//    }

        //
        //
        // END PRINTER NONSENSE
        //
        //

    } // function BCC2012



} // namespace unf_planespanners

#endif // GSNUNF_BCC2012_H

