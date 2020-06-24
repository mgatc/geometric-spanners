#ifndef GSNUNF_CGALCOMPONENTS_H
#define GSNUNF_CGALCOMPONENTS_H

/*
 * Include headers in the following order:
 *      Related header,
 *      C system headers,
 *      C++ standard library headers,
 *      other libraries' headers,
 *      your project's headers.
 */

/**
 *  Global stl includes
 */
#include <string>
#include <unordered_set>
#include <unordered_map>

/**
 *  Global CGAL includes
 */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // CGAL
#include <CGAL/Triangulation_vertex_base_with_info_2.h>         // Triangulations
#include <CGAL/Delaunay_triangulation_2.h>                      // Triangulations
#include <CGAL/point_generators_2.h>                            // Random point generation, testing
#include <CGAL/circulator.h>                                    // Vertex neighbor searches
#include <CGAL/algorithm.h>

/**
 *  Other global includes
 */
#include "Vertex_info.h"

using namespace std;

/**
 *  Global types
 */
typedef CGAL::Exact_predicates_inexact_constructions_kernel                        K;
    typedef K::Point_2                                                         Point;

typedef CGAL::Triangulation_vertex_base_with_info_2< gsnunf::Vertex_info, K >     Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                                 Tds;

typedef CGAL::Delaunay_triangulation_2< K, Tds >               DelaunayTriangulation;
    typedef DelaunayTriangulation::Finite_vertices_iterator Finite_vertices_iterator;
    typedef DelaunayTriangulation::Finite_edges_iterator       Finite_edges_iterator;
    typedef DelaunayTriangulation::Face_handle                           Face_handle;
    typedef DelaunayTriangulation::Vertex_handle                       Vertex_handle;
    typedef DelaunayTriangulation::Vertex_circulator               Vertex_circulator;

typedef CGAL::Creator_uniform_2<double,Point>                                Creator;
typedef CGAL::Container_from_circulator<Vertex_circulator>          Vertex_container;

// For sorting containers of Points, as with std::sort()
// Not used by anything right now, but useful enough to keep around
struct LTR { // comparator for Point type
    bool operator () ( const Point &lhs, const Point &rhs ) {
        return lhs.x() < rhs.x()
            || ( lhs.x() == rhs.x() && lhs.y() < rhs.y() );
    }
};

#endif // GSNUNF_CGALCOMPONENTS_H
