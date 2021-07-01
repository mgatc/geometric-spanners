#ifndef GSNUNF_UTILITIES_H
#define GSNUNF_UTILITIES_H

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp> // hashing pairs
#include <boost/heap/fibonacci_heap.hpp> // ordering

#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Vector_2.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>



namespace gsnunf {

using namespace std;


typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;

typedef K::Point_2                                              Point;
typedef K::Vector_2                                             Vector_2;
typedef K::FT                                                   number_t;

typedef size_t index_t;
typedef size_t cone_t;

typedef pair<size_t,size_t>                             size_tPair;
typedef boost::hash<size_tPair>                         size_tPairHash;
typedef unordered_set<size_tPair,size_tPairHash>        size_tPairSet;
typedef unordered_map<size_tPair,bool,size_tPairHash>   size_tPairMap;

const double EPSILON = 0.000001;
const double INF = std::numeric_limits<double>::infinity();
const size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

const double PI = M_PI;
const double PI_OVER_TWO = PI / 2;
const double SIX_PI_OVER_SEVEN = 6*PI / 7;
const double FOUR_PI_OVER_SEVEN = 4*PI / 7;
const double TAN30 = tan(PI / 6);
const double COS30 = cos(PI / 6);
const double COT30 = 1/TAN30;
const double PI_OVER_FIVE = PI / 5;
const double FOUR_PI_OVER_FIVE = 4*PI / 5;



template< class T >
bool contains( const T& V, const typename T::key_type& v ) {
    return V.find(v) != V.end();
}

template< typename first_t, typename second_t >
std::pair<second_t,first_t> reverse_pair(const std::pair<first_t,second_t>& in)
{
    return std::make_pair( in.second, in.first );
}

/* If V contains v, remove v.
 * If V does not contain v, add it.
 * Return new value.
 */
template< class T >
bool toggle( T& V, const typename T::key_type& v ) {
    bool found = contains( V, v );
    if( found ) V.erase( v );
    else V.insert( v );
    return !found;
}

template< typename T >
inline std::pair<T,T> makeNormalizedPair( const T& i, const T& j ) {
    return make_pair(
        CGAL::min(i,j),
        CGAL::max(i,j)
    );
}

template< typename Point >
inline number_t distance( Point p, Point q ) {
    return p == q ? 0 : CGAL::sqrt( CGAL::squared_distance(p,q) );
}

template< class Point_2 >
inline number_t get_angle( const Point_2& p, const Point_2& q, const Point_2& r ) {
    //CGAL::Vector_2<K> pq( p, q );
    //CGAL::Vector_2<K> rq( r, q );
    auto pq = q - p,
         rq = q - r;

    number_t result = atan2( pq.y(), pq.x() ) - atan2( rq.y(), rq.x() );

    // atan() returns a value between -PI and PI. From zero ("up"), CCW rotation is negative and CW is positive.
    // Our zero is also "up," but we only want positive values between 0 and 2*PI:

    result = fmod( result+2*PI, 2*PI );
    //result = CGAL::min( result, 2*PI );
    //cout<<"angle("<<p<<","<<q<<","<<r<<")="<<result*180/PI<<" ";
    return result;
}
template< class Container >
inline double get_angle( const size_t p, const size_t q, const size_t r, const Container &P )
{
    return get_angle( P[p], P[q], P[r] );
}
template< class K >
inline double get_angle( const typename K::Point_2& p, const typename K::Point_2& q, const typename K::Point_2& r ) {
    return get_angle(p,q,r);
}

template< typename TriWithInfo, typename OutputIterator >
void getVertexInfo( const TriWithInfo& Triangulation, OutputIterator res ) {
    for( auto v=Triangulation.finite_vertices_begin(); v!=Triangulation.finite_vertices_end(); ++v ) {
        *res = to_string( v->info() );
       ++res;
    }
}

struct pointPairHash{
    size_t operator()(const pair<size_t, size_t>& edge) const noexcept{
        size_t seed = 31;
        boost::hash_combine(seed, CGAL::min(edge.first, edge.second));
        boost::hash_combine(seed, CGAL::max(edge.first, edge.second));
        return seed;
    }
};

struct edgeEquality{
    bool operator() (const pair<size_t, size_t> &edgeA, const pair<size_t, size_t> &edgeB) const noexcept{
        return (edgeA.first == edgeB.first && edgeA.second == edgeB.second) || (edgeA.first == edgeB.second && edgeA.second == edgeB.first);
    }
};

struct pointConeHash{
    size_t operator()(const pair<size_t,size_t> &PC) const noexcept{
        size_t seed = 31;
        boost::hash_combine(seed, PC.first);
        boost::hash_combine(seed, PC.second);
        return seed;
    }
};

struct pointConeEquality{
    bool operator()(const pair<size_t,size_t> &PCA, const pair<size_t,size_t> &PCB) const noexcept{
        return PCA.first == PCB.first && PCA.second == PCB.second;
    }
};

template< class K>
void spatialSort(vector<typename K::Point_2> &P, vector<size_t> &index)
{
    typedef CGAL::Spatial_sort_traits_adapter_2<K,
          typename CGAL::Pointer_property_map< typename K::Point_2>::type > Search_traits_2;

    index.clear();
    index.reserve(P.size());

    std::copy( boost::counting_iterator<std::size_t>(0),
               boost::counting_iterator<std::size_t>(P.size()),
               std::back_inserter(index));

    CGAL::spatial_sort( index.begin(),
                        index.end(),
                        Search_traits_2(CGAL::make_property_map(P)) );
    //cout<<"done sorting"<<endl;
}




struct comparatorForMinHeap {
    bool operator()(const size_tPair &n1, const size_tPair &n2) const {
        return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
    }
};
template< class Graph, class OutputIterator>
void reverseLowDegreeOrdering( const Graph &T, OutputIterator out )
{


    typedef boost::heap::fibonacci_heap<size_tPair,boost::heap::compare<comparatorForMinHeap>> Heap;
    typedef Heap::handle_type handle;
    typedef typename Graph::Vertex_circulator Vertex_circulator;

    const size_t n = T.number_of_vertices();

    Heap H;
    vector<handle> handleToHeap(n);
    //vector<size_t> piIndexedByV(n);
    vector<size_t> piIndexedByPiU(n);
    vector<unordered_set<size_t>> currentNeighbors(n);

    // Initialize the vector currentNeighbors with appropriate neighbors for every vertex
    for( auto it=T.finite_vertices_begin();
         it!=T.finite_vertices_end(); ++it )
    {
        Vertex_circulator N = T.incident_vertices(it),
            done(N);
        do {
            if( !T.is_infinite(N) )
                currentNeighbors.at(it->info()).insert( N->info() );
        } while( ++N != done );

        size_t degree = currentNeighbors.at(it->info()).size();
        handleToHeap[it->info()] = H.emplace( degree, it->info() );
    }

    // Use a heap to walk through G_0 to G_{n-1} and set up the Pi for every vertex
    size_t i = n-1; // start at the last valid index

    while( !H.empty() ) {
        size_tPair p = H.top();
        H.pop();
        // make sure our math is correct, e.g., degree from heap key matches neighbor container size
        assert( p.first == currentNeighbors.at( p.second ).size() );
        assert( 0 <= p.first && p.first <= 5 ); // Lemma 1

        // Erase this vertex from incidence list of neighbors and update the neighbors' key in the heap
        for( size_t neighbor : currentNeighbors.at( p.second ) ) {
            currentNeighbors.at(neighbor).erase(p.second);
            handle h = handleToHeap.at(neighbor);
            size_tPair q = make_pair( currentNeighbors.at( neighbor ).size(), neighbor );
            H.update(h,q);
            H.update(h);
        }
        currentNeighbors.at(p.second).clear();
        //piIndexedByV[p.second] = i;
        piIndexedByPiU[i] = p.second;
        --i;
    }

    for( auto v : piIndexedByPiU )
    {
        *out = v;
        ++out;
    }
}


} // namespace gsnunf

#endif // GSNUNF_UTILITIES_H
