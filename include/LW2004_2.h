#ifndef GSNUNF_LW2004_2_H
#define GSNUNF_LW2004_2_H

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
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/random_selection.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>



namespace gsnunf {

using namespace std;

namespace lw2004_2 {

using namespace CGAL;

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K>    Vb;
typedef CGAL::Triangulation_face_base_2<K>                          Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef CGAL::Aff_transformation_2<K>                               Transformation;
typedef CGAL::Vector_2<K> Vector;
typedef Delaunay::Point                                             Point;
typedef Delaunay::Finite_vertices_iterator                          Finite_vertices_iterator;
typedef Delaunay::Finite_edges_iterator                             Finite_edges_iterator;

typedef pair<size_t,size_t>                                         size_tPair;
typedef boost::hash<size_tPair>                                     size_tPairHash;
typedef unordered_set<size_tPair,size_tPairHash>                    size_tPairSet;

class GraphPrinter {
    vector<Point> points;
    size_tPairSet E;
    double radiusOfPoints;

    public:
        GraphPrinter(const vector<Point>& points, const size_tPairSet& E, double radiusOfPoints = 0.09)
            : points(points), E(E), radiusOfPoints(radiusOfPoints) {}

        void draw(string fName = "temp") {
            FILE *fp = fopen(  fName.c_str(),"w");
            fprintf(fp,"\\documentclass{standalone} \n\\usepackage{tikz} \n \n\n\\begin{document}\n");
            fprintf(fp,"\n\n\n\\begin{tikzpicture}\n\n");

            for(pair<size_t,size_t> p : E) {
                double x1 = points[p.first].x();
                double y1 = points[p.first].y();
                double x2 = points[p.second].x();
                double y2 = points[p.second].y();
                fprintf(fp, "\\draw [thick] (%f,%f) -- (%f,%f);\n", x1,y1, x2,y2);
            }

            for(size_t i = 0; i < points.size(); i++)
                fprintf(fp,"\\draw [fill=red,stroke=red] (%f,%f) circle [radius=%f] node[label=above:%zu] {};\n",points[i].x(),points[i].y(),radiusOfPoints,i);

            //for(size_t i = 0; i < points.size(); i++)
              //  fprintf(fp,"\\draw [fill=red,stroke=red] (%f,%f) circle [radius=%f]  {};\n",points[i].x(),points[i].y(),radiusOfPoints);


            fprintf(fp,"\n\n\\end{tikzpicture}");
            fprintf(fp,"\n\n\\end{document}");
            fclose(fp);

            cout << "\nOutput PDF generation started...\n";
            string command = "pdflatex " + fName + " > /dev/null";
            system(command.c_str());
            cout << "PDF generation terminated...\n";

            command = "evince " + fName + ".pdf &";
            system(command.c_str());
        }
};

class TriangulationPrinter {
    double radiusOfPoints;
    Delaunay T;
public:
    TriangulationPrinter(Delaunay T, double radiusOfPoints = 0.09) {
        this-> T = T;
        this->radiusOfPoints = radiusOfPoints;
    }
    void draw(string fName = "temp") {
        FILE *fp = fopen(  fName.c_str(),"w");
        fprintf(fp,"\\documentclass{standalone} \n\\usepackage{tikz} \n \n\n\\begin{document}\n");
        fprintf(fp,"\n\n\n\\begin{tikzpicture}\n\n");

        for(Finite_edges_iterator eit = T.finite_edges_begin();  eit != T.finite_edges_end(); ++eit) {
            Delaunay::Edge e = *eit;
            double x1 = e.first->vertex( (e.second+1)%3 )->point().x();
            double y1 = e.first->vertex( (e.second+1)%3 )->point().y();
            double x2 = e.first->vertex( (e.second+2)%3 )->point().x();
            double y2 = e.first->vertex( (e.second+2)%3 )->point().y();
            fprintf(fp, "\\draw [thick] (%f,%f) -- (%f,%f);\n", x1,y1, x2,y2);
        }

        for(Finite_vertices_iterator it = T.finite_vertices_begin();  it != T.finite_vertices_end(); ++it)
            fprintf(fp,"\\draw [fill=red,stroke=red] (%f,%f) circle [radius=%f] node[label=above:%zu] {};\n",it->point().x(),it->point().y(),radiusOfPoints,it->info());


        fprintf(fp,"\n\n\\end{tikzpicture}");
        fprintf(fp,"\n\n\\end{document}");
        fclose(fp);

        cout << "\nOutput PDF generation started...\n";
        string command = "pdflatex " + fName + " > /dev/null";
        system(command.c_str());
        cout << "PDF generation terminated...\n";

        command = "evince " + fName + ".pdf &";
        system(command.c_str());

    }
};

// Use high degree delaunay for testing
// See if connected
// Find degree
// Find stretch factor
// check planar

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

void generateRandomPoints(vector<Point> &points, const size_t n, const string outputFileName = "") {
    typedef Creator_uniform_2<double,Point>  Creator;
    //Random_points_in_disc_2<Point,Creator> g(10);
    Random_points_in_square_2<Point,Creator> g(10);

    //random_convex_set_2(n,std::back_inserter(points), g);
    // Random_points_on_circle_2<Point,Creator> gen(50);
    // std::copy_n( gen, n, std::back_inserter(points));

    std::copy_n( g, n, std::back_inserter(points));

    if( outputFileName != "") {
        ofstream out;
        out.open(outputFileName);
        for(Point p : points)
            out << p << endl;
        out.close();
    }
}

void readPointsFromFile( vector<Point> &points, const string outputFileName ) {
    ifstream in (outputFileName);
    if (in.is_open()) {
        long double x,y;
        while ( in >> x >> y ) {
            points.push_back(Point(x,y));
        }
        in.close();
    }
}

struct comparatorForMinHeap {
    bool operator()(const size_tPair &n1, const size_tPair &n2) const {
        return (n1.first > n2.first) || ((n1.first == n2.first) && (n1.second > n2.second));
    }
};

typedef boost::heap::fibonacci_heap<size_tPair,boost::heap::compare<comparatorForMinHeap>> Heap;
typedef Heap::handle_type handle;


inline Point rotateByThetaAround(const Point &pivot, const Point &p, const long double theta) {
    Transformation rotate(ROTATION, sin(theta), cos(theta));
    Transformation translate1(TRANSLATION, Vector(-pivot.x(),-pivot.y()));
    Transformation translate2(TRANSLATION, Vector(pivot.x(),pivot.y()));

    Point r = translate1(p);
    r = rotate(r);
    r = translate2(r);

    return r;
}

inline void createNewEdge( size_tPairSet &E, const size_t i, const size_t j, const size_t n ) {
    if( std::max(i,j) > n-1)
        cout <<  "Ooops! out-of-range pointID found! -> " << std::max(i,j) << endl;
    E.insert(make_pair(std::min(i,j), std::max(i,j) ));
}

} // namespace lw2004_2

// alpha is set to pi/2
template< typename RandomAccessIterator, typename OutputIterator >
void LW2004_2( RandomAccessIterator pointsBegin, RandomAccessIterator pointsEnd, OutputIterator result, double alpha = PI/2 ) {
    using namespace lw2004_2;
    Timer t(",");
    vector<Point> points( pointsBegin, pointsEnd );
    const size_t n = points.size();
    //cout << "Step 1 starts...\n";
    Delaunay T;

    T.insert( boost::make_transform_iterator(points.begin(),lw2004_2::AutoCount()),
              boost::make_transform_iterator(points.end(),  lw2004_2::AutoCount()));

    //cout << "Step 1 is over...\n";
    // TriangulationPrinter tp(T);
    // tp.draw("del");
    //************* Step 2 ****************//
    vector<Finite_vertices_iterator> pointID2VertexHandle(n);
    for (auto vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
        pointID2VertexHandle[ vit->info() ] = vit;

    Heap H;
    vector<handle> handleToHeap(n);
    vector<size_t> piIndexedByV(n), piIndexedByPiU(n);
    vector<unordered_set<size_t>> currentNeighbours(n);

   // size_t maxDegree = 0;
    // Initialize the vector currentNeighbours with appropriate neighbours for every vertex
    for(size_t it = 0; it < n; it++) {
        Delaunay::Vertex_circulator vc = T.incident_vertices(pointID2VertexHandle[it]), done(vc);
        if (vc != 0) {
            do {
                if(!T.is_infinite(vc)) {
                    //degree++;
                    currentNeighbours[it].insert(vc->info());
                }
            } while(++vc != done);
        }
       // if(degree > maxDegree)
        //    maxDegree = degree;

        size_t degree = currentNeighbours[it].size();
        handleToHeap[it] = H.push(make_pair(degree,it));
    }

    //cout << "Maximum degree in the Delaunay Triangulation: " << maxDegree << endl;
    // Use a heap to walk through G_0 to G_{n-1} and set up the Pi for every vertex
    size_t i = 1;
    while(!H.empty()) {
        size_tPair p = H.top();
        H.pop();

        for(size_t neighbour : currentNeighbours[p.second]) {
            if(!currentNeighbours[neighbour].empty()) {
                currentNeighbours[neighbour].erase(p.second);
            }
            size_tPair q = make_pair((*handleToHeap[neighbour]).first-1,neighbour);
            H.update(handleToHeap[neighbour],q);
            H.update(handleToHeap[neighbour]);
        }
        currentNeighbours[p.second].clear();
        piIndexedByV[p.second] = n-i;
        piIndexedByPiU[n-i] = p.second;
        i++;
    }

    handleToHeap.clear();
    H.clear();
    currentNeighbours.clear();
    //cout << "Step 2 is over...\n";

    //************* Step 3 ****************//
    // In this step we assume alpha = pi/2 in order to minimize the degree
    size_tPairSet ePrime; // without set duplicate edges could be inserted (use the example down below)
    vector<bool> isProcessed(n, false);

    for(size_t u : piIndexedByPiU) {
        //cout << "\n" << u <<": \n===================\n";

        bool sectorBoundaryVertexFound = false;
        Delaunay::Vertex_circulator vc = T.incident_vertices(pointID2VertexHandle[u]), done(vc), hold(vc);

        // Note: vc++ moves counterclockwise and vc-- moves clockwise
        if (vc != 0) {
            // Using this do-while loop, find the first boundary point adjacent to u
            do {
                if(!T.is_infinite(vc) && (piIndexedByV[vc->info()] < piIndexedByV[u]) && isProcessed[vc->info()]) {
                    sectorBoundaryVertexFound = true;
                    hold = vc;
                    break;
                }
            } while(++vc != done);

            // Now using the boundary point make a circular journey around u
            Point coneStart = vc->point();
            Point coneEnd   = rotateByThetaAround( points[u], coneStart, M_PI_2 );
            long double closestDistanceSquared = LDBL_MAX; // stores the distance to the closest vertex in the current cone
            size_t closestToUInTheCurrentCone; // stores the ID of the point that achieves closestDistanceSquared
            size_t prev; // stores the prev point in the cone
            bool currentConeIsEmptySoFar = true;

            if( !sectorBoundaryVertexFound ) { // every neighbour is unprocessed; this is required otherwise gives incorrect graphs sometimes
                // counterexample:
                // 0.182387 2.19867 8.55689 -4.48489 -8.67484 1.56996
                if(T.is_infinite(hold))
                    hold++;

                Delaunay::Vertex_circulator tempVC = hold, done(tempVC);

                if ( tempVC != 0) {
                    coneStart = tempVC->point();
                    createNewEdge(ePrime,tempVC->info(), u, n);
                    coneEnd = rotateByThetaAround( points[u], coneStart, M_PI_2 );
                    closestDistanceSquared = LDBL_MAX;
                    currentConeIsEmptySoFar = true;
                    tempVC++;

                    while(tempVC != done) {
                        if(T.is_infinite(tempVC) ) {
                            tempVC++;
                            continue;
                        }
                        // lies inside the cone
                        if( left_turn( coneEnd, points[u], tempVC->point()  ) && !left_turn( coneStart, points[u], tempVC->point()  )) {
                            long double distanceFromU = squared_distance(points[u],tempVC->point());

                            if( distanceFromU < closestDistanceSquared) {
                                closestDistanceSquared     = distanceFromU;
                                closestToUInTheCurrentCone = tempVC->info();
                            }

                            if( !currentConeIsEmptySoFar )
                                createNewEdge(ePrime,tempVC->info(), prev, n);

                            currentConeIsEmptySoFar = false;

                            prev = tempVC->info();
                            tempVC++;
                        } else if ( collinear(coneEnd, points[u], tempVC->point()) ) {
                            //cout << "Case 2-Special " << endl;
                            createNewEdge(ePrime,tempVC->info(), u, n);
                            currentConeIsEmptySoFar = true;
                            closestDistanceSquared = LDBL_MAX;
                            coneStart = coneEnd;
                            coneEnd = rotateByThetaAround( points[u], coneStart, M_PI_2 );
                            tempVC++;
                            continue;
                        } else { // outside cone
                            if(!currentConeIsEmptySoFar)
                                createNewEdge(ePrime,closestToUInTheCurrentCone,u, n);

                            coneStart = coneEnd;
                            coneEnd   = rotateByThetaAround( points[u], coneStart, M_PI_2 );
                            currentConeIsEmptySoFar = true;
                            closestDistanceSquared = LDBL_MAX;
                            continue;
                        }
                    }
                    if(!currentConeIsEmptySoFar )
                        createNewEdge(ePrime,closestToUInTheCurrentCone,u, n);
                }
                isProcessed[u] = true;
                continue;
            }

            vc++; // required
            while ( vc != hold ) {
                // do nothing if vc is an infinite vertex
                if( T.is_infinite(vc) ) {
                    vc++;
                    continue;
                }
                // vc is a sector boundary vertex
                if( (piIndexedByV[vc->info()] < piIndexedByV[u]) && isProcessed[vc->info()] ) {
                    coneStart = vc->point();
                    coneEnd = rotateByThetaAround( points[u], coneStart, M_PI_2 );
                    if( !currentConeIsEmptySoFar )
                        createNewEdge(ePrime,u,closestToUInTheCurrentCone, n);

                    closestDistanceSquared = LDBL_MAX;
                    currentConeIsEmptySoFar = true;
                }
                // vc lies on the cone boundary (rarely happens)
                else if( collinear( coneEnd, points[u], vc->point() ) ) {
                    coneStart = vc->point();
                    if( !currentConeIsEmptySoFar )
                        createNewEdge(ePrime,u,closestToUInTheCurrentCone, n);

                    coneEnd = rotateByThetaAround( points[u], coneStart, M_PI_2 );
                    currentConeIsEmptySoFar = true;
                    closestDistanceSquared = LDBL_MAX;
                }
                // vc is inside the current cone having angle <= pi/2 (one side open)
                else if( !isProcessed[u] && right_turn( coneStart, points[u], vc->point() )
                         && left_turn ( coneEnd, points[u], vc->point() ) ) {

                    if(!currentConeIsEmptySoFar)
                        createNewEdge(ePrime,prev,vc->info(), n);
                    else
                        currentConeIsEmptySoFar = false;

                    long double distBetweenUandVC = squared_distance( points[u], vc->point());
                    prev = vc->info();

                    if( distBetweenUandVC < closestDistanceSquared ) {
                        closestDistanceSquared = distBetweenUandVC;
                        closestToUInTheCurrentCone = vc->info();
                    }
                }
                // vc lies outside the cone
                else { //if( right_turn( coneEnd, points[u], vc->point() ) && right_turn( coneStart, points[u], vc->point() ) )
                    if(!currentConeIsEmptySoFar)
                        createNewEdge(ePrime,closestToUInTheCurrentCone,u, n);

                    coneStart = coneEnd;
                    coneEnd = rotateByThetaAround( points[u], coneStart, M_PI_2 );
                    currentConeIsEmptySoFar = true;
                    closestDistanceSquared = LDBL_MAX;
                    continue;
                }

                vc++;
            }
            if(!currentConeIsEmptySoFar )
                createNewEdge(ePrime,closestToUInTheCurrentCone,u, n);
            isProcessed[u] = true;
        }
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

#endif // GSNUNF_LW2004_2_H
