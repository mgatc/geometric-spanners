//
// Created by matt on 7/27/21.
//

#ifndef PLANESPANNERS_POINTGENERATORS_H
#define PLANESPANNERS_POINTGENERATORS_H


#include <iostream>
#include <vector>
#include <chrono>
#include <queue>
#include <limits>
#include <cmath>
#include <random>
#include <unordered_set>


#include <CGAL/point_generators_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/ch_jarvis.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>

#include "utilities.h"

namespace planespanners {

    using namespace std;

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::FT FT;

    const number_t PERTURBATION_VALUE = 0.0001;
//typedef K::Segment_2                                         Segment;
//typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
//typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
//typedef CGAL::Triangulation_data_structure_2<Vb,Fb>          Tds;
//typedef CGAL::Alpha_shape_2<CGAL::Delaunay_triangulation_2<K,Tds>> Alpha_shape;
//
//typedef CGAL::Search_traits_2<K> TreeTraits;
//typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
//typedef Neighbor_search::Tree Tree;

//typedef pair<index_t long,index_t long > Edge;
//
//class Point : public K::Point_2 {
//public:
//    index_t long id = -99; // -99 is a dummy value for id
//    string color = "black";
//    long int v = 0, h = 0;
//
//    Point() = default;
//
//    Point(double X, double Y) : K::Point_2(X, Y) { }
//    Point(double X, double Y, index_t long k) : K::Point_2(X, Y) { id = k; }
//
//    friend ostream &operator<<(ostream &strm, const Point &p) {
//        return strm << p.id << ": (" << p.x() << ", " << p.y() << ")";
//    }
//
//    void setID(index_t long i) {
//        this->id = i;
//    }
//};
//
//inline double L2distance(const Point &p1, const Point &p2) {
//    return std::sqrt(squared_distance(p1,p2));
//}
//
//bool turnOnPdfs = false;

    template<class Container>
    void perturb(Container &P, number_t val) {
        perturb_points_2(P.begin(), P.end(), val, val);
    }


    void generatePointsInsideASquare(const index_t n, const double sizeOfSquare, vector<Point> &P) {
        typedef CGAL::Random_points_in_square_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
        Point_generator g(sizeOfSquare / 2);

        std::copy_n(g, n, back_inserter(P));

        perturb(P, PERTURBATION_VALUE);
    }

    void generatePointsInsideADisc(const index_t n, const double radiusOfDisk, vector<Point> &P) {
        typedef CGAL::Random_points_in_disc_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
        Point_generator g(radiusOfDisk);

        std::copy_n(g, n, back_inserter(P));
        perturb(P, PERTURBATION_VALUE);
    }


    void generatePointsInsideASquareNormal(const index_t pointsInACuster,
                                           const index_t numberOfClusters,
                                           vector<Point> &P,
                                           const number_t xStdDev = 2.0,
                                           const number_t yStdDev = 2.0) {

        std::mt19937 rngX(std::random_device{}());
        std::mt19937 rngY(std::random_device{}());
        std::default_random_engine generatorX(rngX()), generatorY(rngY());
        std::normal_distribution<double> distributionX(0.0, xStdDev), distributionY(2.0, yStdDev);

        std::mt19937 rngShift(std::random_device{}());

        std::uniform_int_distribution shiftDistribution(0, INT32_MAX);

        index_t shiftX, shiftY;
        unordered_set<pair<index_t, index_t>, boost::hash<pair<index_t, index_t>>> S;

        for (index_t c = 0; c < numberOfClusters; c++) {
            if (c != 0) {
                shiftX = shiftDistribution(rngShift) % (20 * numberOfClusters);
                shiftY = shiftDistribution(rngShift) % (20 * numberOfClusters);

                while (!(S.find(make_pair(shiftX, shiftY)) == S.end())) {
                    shiftX = shiftDistribution(rngShift) % (20 * numberOfClusters);
                    shiftY = shiftDistribution(rngShift) % (20 * numberOfClusters);
                }
            } else
                shiftX = shiftY = 0;

            S.insert(make_pair(shiftX, shiftY));

            for (index_t i = 0; i < pointsInACuster; i++) {
                double x = distributionX(generatorX) + shiftX;
                double y = distributionY(generatorY) + shiftY;
                P.emplace_back(K::Point_2(x, y));
            }
        }

        perturb(P, PERTURBATION_VALUE);
    }


    void generateContiguousPointsOnAGrid(const index_t n, vector<Point> &P) {
        points_on_square_grid_2(ceil(std::sqrt(n)), n, std::back_inserter(P), CGAL::Creator_uniform_2<number_t, Point>());
        perturb(P, PERTURBATION_VALUE);
    }

    void generateRandomPointsOnAGrid(const index_t n, vector<Point> &P) {

        unordered_set<pair<int, int>, boost::hash<pair<int, int>>> S;
        std::mt19937 rngX(std::random_device{}());
        std::mt19937 rngY(std::random_device{}());
        std::uniform_int_distribution xDistribution(0, (int) ceil(0.7 * n)), yDistribution(0, (int) ceil(0.7 * n));

        index_t count = 0;

        while (count < n) {
            int x = xDistribution(rngX), y = yDistribution(rngY);

            if (S.find(make_pair(x, y)) == S.end()) {
                P.emplace_back(x, y);
                S.insert(make_pair(x, y));
                count++;
            }
        }

        perturb(P, PERTURBATION_VALUE);
    }

    void generateRandomInsideAnnulus(const index_t n, const double r2, const double r1, vector<Point> &P) {

        assert(r2 > r1);
        //assert(n > 1);

        std::default_random_engine generator;
        std::uniform_real_distribution<double> distributionR(r1, r2), distributionT(0, 1);

        for (index_t i = 0; i < n; i++) {
            double t = 2 * M_PI * distributionT(generator);
            double r = distributionR(generator);
            P.emplace_back(r * cos(t), r * sin(t));
        }

        perturb(P, PERTURBATION_VALUE);
        //assert(P.size()==n);
    }

    void generatePointsFromFile(const index_t n, vector<Point> &P) {

    }

} // planespanners
#endif //PLANESPANNERS_POINTGENERATORS_H
