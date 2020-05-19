#ifndef CGALCOMPONENTS_H
#define CGALCOMPONENTS_H

#include <vector>
#include <list>
#include <algorithm>
#include <ctime>
#include <unordered_set>
#include <fstream>
#include <iterator>
#include <string>
#include <random>
#include <chrono>

//#include <boost/thread/thread.hpp>
/**********************CGAL headers files**************************/
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/random_selection.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2.h>

#include <CGAL/Circular_kernel_intersections.h>
#include <CGAL/Exact_circular_kernel_2.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>

#include <CGAL/Iso_rectangle_2.h>

//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<double> K;

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//typedef CGAL::Exact_predicates_exact_constructions_kernel K;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//
//#include <CGAL/Cartesian.h>
//typedef CGAL::Cartesian<double> K; // throws segmentation fault from LiuLu
/******************************************************************/
using namespace CGAL;
using namespace std;

typedef Delaunay_triangulation_2<K>  DelaunayTriangulation;
typedef K::Point_2                       Point;
typedef K::Segment_2                     Segment;
typedef double                           NumberType;
typedef Creator_uniform_2<NumberType,Point>  Creator;
typedef DelaunayTriangulation::Vertex_handle Vertex_handle;
typedef Circle_2<K> Circle;
typedef Iso_rectangle_2<K> Rectangle;

typedef CGAL::Exact_circular_kernel_2             Circular_k;
typedef CGAL::Point_2<Circular_k>                 PointInCircularKernel;
typedef CGAL::Circle_2<Circular_k>                CircleInCircularKernel;
typedef CGAL::Circular_arc_2<Circular_k>          Circular_arc;
typedef CGAL::CK2_Intersection_traits<Circular_k, CircleInCircularKernel, CircleInCircularKernel>::type Intersection_result;

typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

#endif

/*****************************/
/** Install full boost from synaptic **/
/** Do not forget the flags for compiler and linker **/
