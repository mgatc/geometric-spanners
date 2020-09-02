#ifndef GSNUNF_LW2004_3_H
#define GSNUNF_LW2004_3_H

#include <cmath>
#include <float.h>
#include <forward_list>
#include <fstream>
#include <list>
#include <queue>
#include <unordered_set>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include <CGAL/Aff_transformation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <CGAL/circulator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Vector_2.h>




namespace gsnunf {

using namespace std;

namespace lw2004_3 {

using namespace CGAL;

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K>    Vb;
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

typedef pair<size_t,size_t>                                         size_tPair;
typedef boost::hash<size_tPair>                                     size_tPairHash;
typedef unordered_set<size_tPair,size_tPairHash>                    size_tPairSet;

template<class ArgumentType, class ResultType>
struct unary_funct {
    typedef ArgumentType argument_type;
    typedef ResultType result_type;
};

struct AutoCount : public unary_funct<const Point&,std::pair<Point,size_t> > {
    mutable size_t i;
    AutoCount() : i(0) {}
    pair<Point,size_t> operator()(const Point& p) const {
        return make_pair(p,i++);
    }
};

struct comparatorForMinHeap {
    bool operator()(const size_tPair &n1, const size_tPair &n2) const {
        return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
    }
};

typedef boost::heap::fibonacci_heap<size_tPair,boost::heap::compare<comparatorForMinHeap>> Heap;
typedef Heap::handle_type handle;

inline void createNewEdge( const Delaunay& T, const vector<Delaunay::Vertex_handle>& handles, size_tPairSet &E, const size_t i, const size_t j, const size_t n ) {
    // need access to T and pointID2VertexHandle
    assert( T.is_edge( handles.at(i), handles.at(j) ) );

    if( std::max(i,j) > n-1)
        cout <<  "Ooops! out-of-range pointID found! -> " << std::max(i,j) << endl;
    E.insert(make_pair(std::min(i,j), std::max(i,j) ));
}

template< class DT >
double get_angle( const DT& T, typename DT::Vertex_handle &p, const typename DT::Vertex_handle &q, const typename DT::Vertex_handle &r ) {
    assert( !T.is_infinite(p) );
    assert( !T.is_infinite(q) );
    assert( !T.is_infinite(r) );

    Vector_2 pq( p->point(), q->point() );
    Vector_2 rq( r->point(), q->point() );

    double result = atan2( rq.y(), rq.x() ) - atan2( pq.y(), pq.x() );

    // atan() returns a value between -PI and PI. From zero ("up"), CCW rotation is negative and CW is positive.
    // Our zero is also "up," but we only want positive values between 0 and 2*PI:

    result *= -1; // First, invert the result. This will associate CW rotation with positive values.
    if( result < EPSILON ) // Then, if the result is less than 0 (or epsilon for floats) add 2*PI.
        result += 2*PI;

    return CGAL::min( result, 2*PI );
}

} // namespace lw2004_2

// alpha is set to pi/2
template< typename RandomAccessIterator, typename OutputIterator >
void LW2004_3( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, double alpha = PI/2 ) {
    using namespace lw2004_3;

    // ensure valid alpha
    alpha = CGAL::max( EPSILON, CGAL::min( alpha, PI/2 ) );

    // add points to vector and remove any duplicates or AutoCount() will fail
    vector<Point> points( pointsBegin, pointsEnd );
    std::sort( points.begin(), points.end() );
    auto last = std::unique( points.begin(), points.end() );
    points.erase( last, points.end() );

    const size_t n = points.size();
//    cout<<n<<",";
    //cout << "Step 1 starts...\n";
    Delaunay T;

    T.insert( boost::make_transform_iterator(points.begin(),AutoCount()),
              boost::make_transform_iterator(points.end(),  AutoCount()));

    Delaunay::Vertex_handle v_inf = T.infinite_vertex();

    //cout << "Step 1 is over...\n";
    // TriangulationPrinter tp(T);
    // tp.draw("del");
    //************* Step 2 ****************//
    vector<Delaunay::Vertex_handle> pointID2VertexHandle(n, v_inf);
    for( auto vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit ) {
        //assert( !T.is_infinite(vit) );
        pointID2VertexHandle[ vit->info() ] = vit;
    }

    Heap H;
    vector<handle> handleToHeap(n);
    vector<size_t> piIndexedByV(n), piIndexedByPiU(n);
    vector<unordered_set<size_t>> currentNeighbours(n);

   // size_t maxDegree = 0;
    // Initialize the vector currentNeighbours with appropriate neighbours for every vertex
    for( size_t it = 0; it < n; it++ ) {
        Delaunay::Vertex_circulator N = T.incident_vertices( pointID2VertexHandle.at(it) ),
            done(N);
        //if (vc != 0) {
        do {
            if( !T.is_infinite(N) ) {
                //degree++;
                currentNeighbours.at(it).insert( N->info() );
            }
        } while( ++N != done );
        //}
       // if(degree > maxDegree)
        //    maxDegree = degree;

        size_t degree = currentNeighbours.at(it).size();
        handleToHeap[it] = H.emplace( degree,it );
    }

    //cout << "Maximum degree in the Delaunay Triangulation: " << maxDegree << endl;
    // Use a heap to walk through G_0 to G_{n-1} and set up the Pi for every vertex
    size_t i = n-1; // start at the last valid index
    while(!H.empty()) {
        size_tPair p = H.top();
        H.pop();

        for( size_t neighbour : currentNeighbours.at( p.second ) ) {
            currentNeighbours.at(neighbour).erase(p.second);

            handle h = handleToHeap.at(neighbour);
            size_t degree = (*h).first;
            assert( degree > 0 );
            size_tPair q = make_pair( degree-1, neighbour );
            H.update(h,q);
            H.update(h);
        }
        currentNeighbours.at(p.second).clear();
        piIndexedByV[p.second] = i;
        piIndexedByPiU[i] = p.second;
        --i;
    }

    handleToHeap.clear();
    H.clear();
    currentNeighbours.clear();
    //cout << "Step 2 is over...\n";

    //************* Step 3 ****************//
    // In this step we assume alpha = pi/2 in order to minimize the degree
    size_tPairSet ePrime; // without set duplicate edges could be inserted (use the example down below)
    vector<bool> isProcessed(n, false);
    Delaunay::Vertex_handle u_handle = v_inf;

    //cout<<piIndexedByPiU.size()<<"\n";

    // Iterate through vertices by pi ordering
    for( size_t u : piIndexedByPiU ) {
        u_handle = pointID2VertexHandle.at(u);
        assert( !T.is_infinite(u_handle) );
        //cout<<u_handle->point()<<"\n";
        // Get neighbors of u
        Delaunay::Vertex_circulator N = T.incident_vertices( u_handle );
        while( T.is_infinite(N) ) ++N; // Make sure N isn't infinite to start with
        Delaunay::Vertex_circulator done(N);
        // Rotate N until reaching a processed vertex or the original vertex
        while( ( T.is_infinite( ++N ) || !isProcessed.at( N->info() ) ) && N != done );
        // Find and store sector boundaries, start with N
        vector<Delaunay::Vertex_handle> sectorBoundaries{ N->handle() };
        while( ++N != sectorBoundaries.front() ) {
            if( !T.is_infinite(N) && isProcessed.at(N->info()) )
                sectorBoundaries.push_back( N->handle() );
        }
        assert( sectorBoundaries.size() <= 5 );

        // Now, compute the angles of the sectors, the number of cones in each sector,
        // and the actual angles
        vector<double> alphaReal( sectorBoundaries.size() );
        vector< vector<Delaunay::Vertex_handle> > closest( sectorBoundaries.size() );

        for( size_t i=0; i<sectorBoundaries.size(); ++i ) {
            double sectorAngle = get_angle(
                T,
                sectorBoundaries.at(i),
                u_handle,
                sectorBoundaries.at( (i+1)%sectorBoundaries.size() )
            );
            size_t numCones = rint( ceil( sectorAngle / alpha ) );
            assert( numCones > 0 ); // guard against /0
            alphaReal[i] = sectorAngle / numCones;
            closest.at(i).resize( numCones, v_inf );
        }

        Delaunay::Vertex_handle lastN = v_inf;
        if( isProcessed.at( N->info() ) ) ++N; // if current neighbor is processed, step
        size_t sector = 0;

        do { // Loop through neighbors and add appropriate edges
            if( !T.is_infinite(N) ) {
                if( isProcessed.at( N->info() ) ) {
                    ++sector;
                } else {
                    // evaluate possible forward edges
                    double theta = get_angle(
                        T,
                        sectorBoundaries.at(sector),
                        u_handle,
                        N->handle()
                    );
                    // Due to floating point inaccuracies, a sector of size 360
                    // can result in an out-of-range cone.
                    size_t cone = CGAL::min(
                        size_t( theta / alphaReal.at(sector) ),
                        closest.at(sector).size()-1
                    );
                    // Store value until after all neighbors are processed, then add
                    if( T.is_infinite( closest.at(sector).at(cone) )
                      || Vector_2( u_handle->point(), N->point() ).squared_length()
                       < Vector_2( u_handle->point(), closest.at(sector).at(cone)->point() ).squared_length() )
                            closest.at(sector).at(cone) = N->handle();   // if the saved vertex is infinite or longer than the current one, update

                    // cross edges
                    if( !T.is_infinite( lastN ) && !isProcessed.at( lastN->info() ) )
                        createNewEdge( T, pointID2VertexHandle, ePrime, lastN->info(), N->info(), n );
                }
            }
            lastN = N->handle();
        } while( ++N != done );

        // Add edges in closest
        for( auto segment : closest )
            for( auto v : segment )
                if( !T.is_infinite(v) )
                    createNewEdge( T, pointID2VertexHandle, ePrime, u, v->info(), n );
    }
    for( size_tPair e : ePrime ) {
        *result = make_pair( points[e.first], points[e.second] );
        ++result;
        *result = make_pair( points[e.second], points[e.first] );
        ++result;
    }
    //cout << "------------------------\nEdges: " << ePrime.size() << endl;
    //GraphPrinter gp(points, ePrime);
    //gp.draw();

}

} // namespace gsnunf

#endif // GSNUNF_LW2004_3_H
