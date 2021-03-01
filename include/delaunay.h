#ifndef GSNUNF_DELAUNAY_H
#define GSNUNF_DELAUNAY_H

#include <algorithm> // min, max
#include <cmath> // ceil
#include <unordered_set> // hashed adjacency list
#include <vector> // vertex containers

#include <boost/functional/hash.hpp> // hashing pairs
#include <boost/heap/fibonacci_heap.hpp> // ordering

//#include <CGAL/algorithm.h>
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/DelaunayTriangulationTraits_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/utils.h> // min, max
#include <CGAL/Vector_2.h>

#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "utilities.h"

namespace gsnunf {

using namespace std;

namespace delaunay {

    using namespace CGAL;

    typedef CGAL::Exact_predicates_inexact_constructions_kernel         Epick;

    template< class Geom_traits >
    class StandardDelaunayTriangulationTraits_2 : public Geom_traits {
      public:
        typedef typename Geom_traits::Point_2       Point_2;
        typedef typename Geom_traits::FT            FT;

        struct Side_of_oriented_circle_2 {
            Oriented_side operator()( const Point_2& p,
                                      const Point_2& q,
                                      const Point_2& r,
                                      const Point_2& s ) {

                /* IF USING THIS CLASS AS A MODEL TO IMPLEMENT
                 * AN EXOTIC DELAUNAY CONDITION, IMPLEMENT YOUR
                 * PREDICATE HERE AND RETURN THE APPROPRIATE
                 * Oriented_side
                 */

                FT px = p.x();
                FT py = p.y();
                FT qx = q.x();
                FT qy = q.y();
                FT rx = r.x();
                FT ry = r.y();
                FT tx = s.x();
                FT ty = s.y();
                //  sign_of_determinant(px, py, px*px + py*py, 1,
                //                      qx, qy, qx*qx + qy*qy, 1,
                //                      rx, ry, rx*rx + ry*ry, 1,
                //                      tx, ty, tx*tx + ty*ty, 1);
                // We first translate so that p is the new origin.
                FT qpx = qx-px;
                FT qpy = qy-py;
                FT rpx = rx-px;
                FT rpy = ry-py;
                FT tpx = tx-px;
                FT tpy = ty-py;
                // The usual 3x3 formula can be simplified a little bit to a 2x2.
                //         - sign_of_determinant(qpx, qpy, square(qpx) + square(qpy),
                //                                  rpx, rpy, square(rpx) + square(rpy),
                //                                  tpx, tpy, square(tpx) + square(tpy)));
                return sign_of_determinant<FT>( qpx*tpy - qpy*tpx, tpx*(tx-qx) + tpy*(ty-qy),
                                              qpx*rpy - qpy*rpx, rpx*(rx-qx) + rpy*(ry-qy));
            }
        };
        Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const {
            return Side_of_oriented_circle_2();
        }
    };

    template< class Geom_traits >
    class TriangularDistanceDelaunayTriangulationTraits_2 : public Geom_traits {
      public:
        typedef typename Geom_traits::Point_2       Point_2;
        typedef typename Geom_traits::FT            FT;
        typedef typename Geom_traits::Line_2        Line_2;



        struct Side_of_oriented_circle_2 {
            Oriented_side operator()( const Point_2& p,
                                      const Point_2& q,
                                      const Point_2& r,
                                      const Point_2& s ) {

                /* IF USING THIS CLASS AS A MODEL TO IMPLEMENT
                 * AN EXOTIC DELAUNAY CONDITION, IMPLEMENT YOUR
                 * PREDICATE HERE AND RETURN THE APPROPRIATE
                 * Oriented_side
                 */

//                FT px = p.x();
//                FT py = p.y();
//                FT qx = q.x();
//                FT qy = q.y();
//                FT rx = r.x();
//                FT ry = r.y();
//                FT sx = s.x();
//                FT sy = s.y();

                vector<pair<Point_2, Point_2>> triEdges = {{p,q}, {p,r}, {q,r}};
                const vector<double> orthBisectorSlopes{ 0, -1*COT30, COT30, 0, -1*COT30, COT30 };

                //Finds the cone of p containing vertex q, for this algorithm all vertices have 6 cones (0-5) with an angle of (PI/3).
                const double alpha = PI/3; // 60 Degree for equalateral triangle.
                Point_2 focus, opPoint, refPoint;
                double theta;
                size_t cone, coneS;

                for(auto e: triEdges){
                    // a w.r.t. b
                    focus = e.first;
                    opPoint = e.second;
                    refPoint = Point_2( e.first.x() - TAN30, e.first.y() + 1 );
                    theta = get_angle<Geom_traits>(refPoint, e.first, e.second);
                    cone = (theta / alpha);

                    if(cone % 2){
                        // b W.r.t a
                        focus = e.second;
                        opPoint = e.first;
                        refPoint = Point_2( e.second.x() - TAN30, e.second.y() + 1 );
                        theta = get_angle<Geom_traits>(refPoint, e.second, e.first);
                        cone = (theta / alpha);

                    }

                    //Find where s is
                    theta = get_angle<Geom_traits>(refPoint, focus, s);
                    coneS = (theta / alpha);

                    if(cone == coneS){

                        double bisectorLenOP = bisectorLength({focus, opPoint}, orthBisectorSlopes[cone]);
                        double bisectorLenS = bisectorLength({focus, s}, orthBisectorSlopes[cone]);

                        if( bisectorLenS < bisectorLenOP){
                            return ON_POSITIVE_SIDE;
                        }


                    }


                }

                return ON_NEGATIVE_SIDE;

            }

            //Finds the bisector length of a given edge.
            inline typename Geom_traits::FT bisectorLength(const pair<Point_2,Point_2> &e, double slope) {

                typename Geom_traits::FT xCord = e.first.x();
                typename Geom_traits::FT yCord = e.first.y() + 1;

                xCord = xCord - slope;

                Point_2 bisectorPoint(xCord, yCord);

                Line_2 bisectorLine(e.first, bisectorPoint);

                Point_2 intersectionPoint = bisectorLine.projection(e.second);

                typename Geom_traits::FT bisectorLen = gsnunf::distance(e.first, intersectionPoint);

                return bisectorLen;
            }
        };
        Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const {
            return Side_of_oriented_circle_2();
        }
    };



    typedef TriangularDistanceDelaunayTriangulationTraits_2<Epick>                K;
    //typedef Epick                                                       K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K>      Vb;
    typedef CGAL::Triangulation_face_base_2<K>                          Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
    typedef CGAL::Aff_transformation_2<K>                               Transformation;
    typedef Delaunay::Vertex_handle                                     Vertex_handle;
    typedef Delaunay::Vertex_circulator                                 Vertex_circulator;
    typedef CGAL::Vector_2<K>                                           Vector_2;
    typedef Delaunay::Finite_vertices_iterator                          Finite_vertices_iterator;
    typedef Delaunay::Finite_edges_iterator                             Finite_edges_iterator;
    typedef K::Line_2                                            Line;
    typedef K::Point_2                                           Point;

}

// alpha is set to pi/2
template< typename RandomAccessIterator, typename OutputIterator >
void delaunay_testing( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result ) {
    using namespace delaunay;

    Delaunay T( pointsBegin, pointsEnd );







    //
    //
    // START PRINTER NONSENSE
    //
    //


        GraphPrinter printer(5);
        GraphPrinter::OptionsList options;

        options = {
            { "color", printer.inactiveEdgeColor },
            { "line width", to_string(printer.inactiveEdgeWidth) }
        };
        printer.drawEdges( T, options );

//        options = { // active edge options
//            { "color", printer.activeEdgeColor },
//            { "line width", to_string(printer.activeEdgeWidth) }
//        };
//        printer.drawEdges( edgeList.begin(), edgeList.end(), options );


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

        printer.print( "delaunay2" );
        cout<<"\n";




    //
    //
    // END PRINTER NONSENSE
    //
    //

} // function delaunay_testing

} // namespace gsnunf

#endif // GSNUNF_DELAUNAY_H

