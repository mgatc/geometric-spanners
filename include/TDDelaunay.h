#ifndef TD_DELAUNAY_H
#define TD_DELAUNAY_H


#include <CGAL/Compact_container.h>
#include <CGAL/Construct_theta_graph_2.h>
#include <CGAL/property_map.h>

#include <boost/functional/hash.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include "utilities.h"

//namespace std {
//    template < class Gt >
//    struct hash<Vertex_base<Gt>> {
//        size_t operator()(const Vertex_base<Gt>& v) const {
//            size_t seed = 0;
//
//            boost::hash_combine(seed, boost::hash_value(v->point()));
//
//            return seed;
//        }
//    };
//}

namespace gsnunf {

//template< class Gt >
//class Vertex_base {
//  public:
//    typedef typename Gt::Point_2 Point_2;
//    Vertex_base() {}
//    Vertex_base(const Point_2& p) : _p(p) {}
//    void set_point(const Point_2& p) { _p = p; }
//    const Point_2&  point() const { return _p; }
//    friend bool operator<(const Vertex_base<Gt>& l, const Vertex_base<Gt>& r ) {
//        return l->point() < r->point();
//    }
//    bool operator==(const Vertex_base<Gt>& other) const {
//        return _p == other->point();
//    }
//  private:
//    Point_2 _p;
//};
//class Vertex_circulator {
//
//};
//Cone angles.
const number_t tan30 = TAN30;
const number_t cot30 = 1 / tan30;

const number_t alpha = PI/3;

//Slopes of the cone boundary lines.
const vector<number_t> bisectorSlopes{ INF, tan30, -1*tan30, INF, tan30, -1*tan30 };
const vector<number_t> orthBisectorSlopes{ 0, -1*cot30, cot30, 0, -1*cot30, cot30 };

//Finds the cone of p containing vertex q, for this algorithm all vertices have 6 cones (0-5) with an angle of (PI/3).
template<class Point>
inline size_t getSingleCone(const size_t p, const size_t q, const vector<Point> &h)
{
    if( CGAL::compare_y(h[p], h[q]) == CGAL::EQUAL ) {
        return 1 + 3*int(CGAL::compare_x(h[p], h[q]) == CGAL::LARGER);
    }
    const Point refPoint( h.at(p).x() - tan30, h[p].y() + 1 );
    //Point refPoint(h[p]->point().x(), h[p] ->point().y() + 1);

    number_t theta = get_angle(refPoint, h[p], h.at(q));

    size_t cone = (theta / alpha);

    return cone;
}

//Compute max of getCone(p,q) and (getCone(q,p)+3)%6, is used to make sure cones are calculated correctly.
template<class Point>
inline size_t getCone( const size_t p, const size_t q, const vector<Point> &h )
{
    return p < q ?
        getSingleCone(p,q,h)
        : ( getSingleCone(q,p,h)+3 ) % 6;
}

template< class Gt >
class HalfThetaTriangulation {
  public:

    // select the kernel type
    typedef Gt                                              Geom_traits;
    typedef typename Geom_traits::Point_2                   Point_2;
    typedef typename Geom_traits::Direction_2               Direction_2;
   // typedef Vertex_base<Geom_traits>                        Vertex_handle;
    typedef typename std::vector<Point_2>::iterator         Point_iterator;
    //typedef typename std::vector<Vertex_handle>::iterator   Vertex_iterator;
    typedef boost::adjacency_list<boost::vecS,
                              boost::vecS,
                              boost::bidirectionalS,
                              Point_2,
                              size_t>                       Graph;
    typedef typename boost::graph_traits<Graph>::edge_iterator     Edge_iterator;
    typedef typename boost::graph_traits<Graph>::in_edge_iterator  In_edge_iterator;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator Out_edge_iterator;
    typedef typename boost::graph_traits<Graph>::edge_descriptor   Edge_descriptor;
    typedef typename boost::graph_traits<Graph>::vertex_iterator   Vertex_iterator;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex_descriptor;

    // define the Graph

    Graph _G;
    CGAL::Construct_theta_graph_2<Geom_traits, Graph> halfTheta;

    template <typename InputIt>
    HalfThetaTriangulation(InputIt first, InputIt last, const Geom_traits& gt = Geom_traits()) :
                halfTheta(6, Direction_2(1,0), CGAL::ODD_CONES),
                _gt(gt)
    {
        insert(first,last);
    }

    HalfThetaTriangulation(const Geom_traits& gt = Geom_traits()) :
                halfTheta(6, Direction_2(1,0), CGAL::ODD_CONES),
                _gt(gt)
    {

    }

    template < class InputIt>
    std::ptrdiff_t insert(InputIt first, InputIt last) {
        auto n = number_of_vertices();

        _P.insert( _P.end(), first, last );
        size_t id=0;
        for(auto pit=_P.begin(); pit!=_P.end(); pit++) {
            info.emplace(pit, id++);
        }
        halfTheta( _P.begin(), _P.end(), _G);          // Construct the half-theta graph


        return number_of_vertices() - n;
    }

    const Geom_traits& geom_traits() const { return _gt;}

    std::ptrdiff_t number_of_vertices() {
        return boost::num_vertices(_G);
    }
    Vertex_iterator finite_vertices_begin() const {
        return boost::vertices(_G).first;
    }
    // finite_vertices_end
    Vertex_iterator finite_vertices_end() const {
        return boost::vertices(_G).second;
    }
    // finite_points_begin
    Point_iterator points_begin() {
        return _P.begin();
    }
    // finite_points_end
    Point_iterator points_end() {
        return _P.end();
    }
    // finite_edges_begin
    Edge_iterator edges_begin() const {
        return boost::edges(_G).first;
    }
    // finite_edges_end
    Edge_iterator edges_end() const {
        return boost::edges(_G).second;
    }
    Point_2& point(size_t index) {
        return _G[index];
    }
    // get negative cone edges
    In_edge_iterator negative_cone_edges_begin(const Vertex_descriptor& u) const {
        return in_edges(u, _G).first;
    }
    In_edge_iterator negative_cone_edges_end(const Vertex_descriptor& u) const {
        return in_edges(u, _G).second;
    }
    // get positive cone edges
    Out_edge_iterator positive_cone_edges_begin(const Vertex_descriptor& u) const {
        return out_edges(u, _G).first;
    }
    Out_edge_iterator positive_cone_edges_end(const Vertex_descriptor& u) const {
        return out_edges(u, _G).second;
    }

    //Compute max of getCone(p,q) and (getCone(q,p)+3)%6, is used to make sure cones are calculated correctly.
    inline size_t getCone( const Vertex_descriptor &p, const Vertex_descriptor &q ) const {
        return gsnunf::getCone(p,q,_P);
    }
    inline bool edgeExists(const pair<Vertex_descriptor,Vertex_descriptor>& e) const {
        return boost::edge(e.first, e.second,_G).second;
    }
    inline pair<pair<Vertex_descriptor,Vertex_descriptor>, bool>
    eitherEdge(const Vertex_descriptor& u, const Vertex_descriptor& v) const
    {
        auto orientedEdge = make_pair(u,v),
             reversedEdge = reverse_pair(orientedEdge);
        bool reversedEdgeExists = edgeExists(reversedEdge);
        return make_pair(
            reversedEdgeExists ? reversedEdge : orientedEdge,
            reversedEdgeExists || edgeExists(orientedEdge)
        );
    }
    template< class EdgeList >
    void fan_of_cone(const Vertex_descriptor& v, const size_t cone, EdgeList &fan) const
    {
        // get in_edges of v
        for( auto neit=negative_cone_edges_begin(v);
             neit!=negative_cone_edges_end(v); ++neit )
        {
            auto e = *neit;
            if( getCone( target(e), source(e) ) == cone )
            {
                fan.insert( fan.end(), e );
            }
        }
        // sort edges in counter-clockwise order
        sort( fan.begin(), fan.end(),
            [&] ( const auto &lhs, const auto &rhs )
            {
                assert( target(lhs) == target(rhs) );
                const Point_2 refPoint( _P[target(lhs)].x()-TAN30,
                                        _P[target(lhs)].y()+1 );
                return get_angle( refPoint,
                                  _P[target(lhs)],
                                  _P[source(lhs)])
                     > get_angle( refPoint,
                                  _P[target(rhs)],
                                  _P[source(rhs)]);
            }
        );
//        cout<< "Fan "<<cone<< " of "<<v<<": ";
//        for( auto e : fan )
//        {
//            cout<<source(e)<<"-"<<target(e)<<" ";
//        }
//        cout<<"\n";
    }
    Vertex_descriptor source( const Edge_descriptor& e ) const {
        return boost::source(e,_G);
    }
    Vertex_descriptor target( const Edge_descriptor& e ) const {
        return boost::target(e,_G);
    }
    Vertex_descriptor parent( const Vertex_descriptor child, const size_t i ) const {
        for( auto it=positive_cone_edges_begin(child);
             it!=positive_cone_edges_end(child); ++it )
        {
            Vertex_descriptor v = target(*it);
            if( getCone(child,v) == i )
                return v;
        }
        return SIZE_T_MAX;
    }


  protected:
    Geom_traits _gt;
    //Vertex_handle _infinite_vertex;
    std::vector<Point_2> _P;
    std::vector<Vertex_descriptor> _V;
    std::map<Point_iterator,size_t> info;

  private:
    const double alpha = PI / 3;
};


}

#endif //TD_DELAUNAY_H
