//Needs optimizing currently testing.
#ifndef GSNUNF_BGHP2010_H
#define GSNUNF_BGHP2010_H

//Base libraries.
#include <cmath>         // ceil, floor, isinf
#include <limits>
#include <unordered_set> // selected
#include <unordered_map> // G_prime
#include <vector>        // handles

//Boost library
#include <boost/functional/hash.hpp> // size_t pair hash

//CGAL library
#include <CGAL/algorithm.h>
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Line_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

//Project library
#include "GeometricSpannerPrinter.h"
//#include "GraphAlgoTV.h"
#include "metrics.h"
#include "TDDelaunay.h"
#include "utilities.h"


namespace gsnunf {

using namespace std;

namespace bghp2010 {

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
typedef Epick                     K;
typedef K::Point_2                Point_2;
typedef K::FT                     FT;

typedef HalfThetaTriangulation<K> TD_Delaunay_2;
typedef TD_Delaunay_2::Vertex_descriptor Vertex_descriptor;

//Finds the bisector length of a given edge.
inline K::FT bisectorLength( const pair<size_t,size_t> &e, const vector<Point_2> &h )
{
    size_t cone = getCone(e.first, e.second, h);

    double xCord = h.at(e.first).x();
    double yCord = h[e.first].y() + 1;

    assert(cone<6);
    assert(e.first<h.size());

    xCord = h[e.first].x() - orthBisectorSlopes.at(cone);

    Point_2 bisectorPoint(xCord, yCord);

    K::Line_2 bisectorLine(h[e.first], bisectorPoint);

    Point_2 intersectionPoint = bisectorLine.projection(h[e.second]);

    double bisectorLen = distance(h[e.first], intersectionPoint);

    return bisectorLen;
}



} // namespace BGHP2010


// Main algorithm.
template<typename RandomAccessIterator, typename OutputIterator>
void BGHP2010(RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, bool printLog = false)
{
    using namespace bghp2010;

    vector<Point_2> P( pointsBegin, pointsEnd );

    TD_Delaunay_2 D( P.begin(), P.end() );
    {
        //Timer tim;
        size_t n = D.number_of_vertices();
















        // Send resultant graph to output iterator
//        for(auto e : S)
//        {
//            *result = e;
//            ++result;
////            *result = reverse_pair(e);
////            ++result;
//        }





        // START PRINTER NONSENSE
        if(printLog)
        {
            vector<pair<Point_2,Point_2>> edgeList;

            /*for(auto e : S)
            {
                edgeList.emplace_back(P.at(e.first), P.at(e.second));
            }*/

            GraphPrinter printer(0.01);
            GraphPrinter::OptionsList options;

            options = {
                {"color", printer.inactiveEdgeColor},
                {"line width", to_string(printer.inactiveEdgeWidth)}
            };
            printer.drawEdgesOfHalfTheta(D, options);

            options = { // active edge options
                {"color", printer.activeEdgeColor},
                {"line width", to_string(printer.activeEdgeWidth)}
            };
            //printer.drawEdges(edgeList.begin(), edgeList.end(), options);

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
            printer.drawVerticesWithInfo(D.points_begin(), D.points_end(), options, borderOptions);

            printer.print("BGHP2010");
            cout << "\n";
        }
        // END PRINTER NONSENSE
    }
} // function BGHP2010

} // namespace gsnunf

#endif // GSNUNF_BGHP2010_H
