#ifndef SPANNERS_POINTGENERATORS_H
#define SPANNERS_POINTGENERATORS_H

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <unordered_set>

#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>

#include "tools/Utilities.h"

namespace spanners {

    using namespace std;

    enum DistributionType {
        DistributionTypeFirst = 0,
        Synthetic = DistributionTypeFirst,
        Real,
        DistributionTypeLast
    };

    enum SyntheticDistribution {
        SyntheticDistributionFirst=0,
//        UniformInsideSquare = SyntheticDistributionFirst,
//        UniformInsideDisc,
////        UniformOnSquare,
////        UniformOnCircle,
//        NormalInsideSquare,
//        NormalClustersInsideSquare,
//        ContiguousGrid,
//        UniformRandomGrid,
//        UniformInsideAnnulus,
//        Galaxy,
        ConvexHullInDisc = SyntheticDistributionFirst,
        SyntheticDistributionLast
    };

    const vector<string> SYNTHETIC_DISTRIBUTION_NAMES = {
//            "Uniform Inside Square",
//            "Uniform Inside Disc",
////            "Uniform On Square",
////            "Uniform On Circle",
//            "Normal Inside Square",
//            "Normal Inside Square with Clusters",
//            "Contiguous Grid",
//            "Uniform Random Grid",
//            "Uniform Inside Annulus",
//            "Galaxy",
            "Convex Hull In Disc",
    };
    vector<string> REAL_POINTSET_NAMES;

    class PointGenerator_2 {

    public:
        PointGenerator_2() : m_randCgal(std::rand()) {}

        void loadFromFile(const string filename, vector<Point> &P) {
            std::ifstream in(filename);

            if (!in.is_open())
                cout<<"Error opening file!\n";

            number_t x,y;
            while ( in >> x >> y ) {
                P.emplace_back(x,y);
            }
            in.close();

            perturb(P, m_perturbationValue);
        }

        void insideSquare(const index_t n, const double sizeOfSquare, vector<Point> &P) {
            typedef CGAL::Random_points_in_square_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            Point_generator g(sizeOfSquare / 2, m_randCgal);

            unordered_set<Point> P_unique;
            size_t remaining;
            while((remaining = n - P_unique.size()) > 0)
                std::copy_n(g, n, inserter(P_unique, P_unique.end()));

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }
        void onSquare(const index_t n, const double sizeOfSquare, vector<Point> &P) {
            typedef CGAL::Random_points_on_square_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            Point_generator g(sizeOfSquare / 2, m_randCgal);

            unordered_set<Point> P_unique;
            size_t remaining;
            while((remaining = n - P_unique.size()) > 0)
                std::copy_n(g, n, inserter(P_unique, P_unique.end()));

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
        }
        void onCircle(const index_t n, const double sizeOfSquare, vector<Point> &P) {
            typedef CGAL::Random_points_on_circle_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            Point_generator g(sizeOfSquare / 2, m_randCgal);

            unordered_set<Point> P_unique;
            size_t remaining;
            while((remaining = n - P_unique.size()) > 0)
                std::copy_n(g, n, inserter(P_unique, P_unique.end()));

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
        }

        void onSpokes(const index_t n, const unsigned numSpokes, vector<Point> &P) {
            //srand(seed());
            double spokeAngle = 2*PI / numSpokes;

            unordered_set<Point> P_unique;

            while( P_unique.size() < n ) {
                double distance = randFloat();
                double angle = randFloat() * 2 * PI;
                angle = ((unsigned) (angle / spokeAngle) )*spokeAngle;
                P_unique.emplace(cos(angle) * distance,
                                 sin(angle) * distance);
            }

            copy(P_unique.begin(),P_unique.end(),back_inserter(P));
        }
        void inGalaxy(const index_t n, const unsigned numSpokes, vector<Point> &P) {
            // see https://itinerantgames.tumblr.com/post/78592276402/a-2d-procedural-galaxy-with-c
            //srand(seed());
            const double spokeAngle = 2*PI / numSpokes,
                    armOffsetMax = 0.5,
                    rotationFactor = 5,
                    perturbationValue = 0.02;

            unordered_set<Point> P_unique;

            while( P_unique.size() < n ) {
                //for(index_t i=0; i<n; ++i) {
                double distance = randFloat();
                distance = pow(distance,2);

                double angle = randFloat() * 2 * PI;
                double armOffset = randFloat() * armOffsetMax;
                armOffset -= armOffsetMax / 2;
                armOffset *= (1/distance);

                double squaredArmOffset = pow(armOffset,2);
                squaredArmOffset *= -1 * int(armOffset < 0);
                armOffset = squaredArmOffset;

                double rotation = distance * rotationFactor;

                angle = ((unsigned) (angle / spokeAngle) )*spokeAngle;
                angle += armOffset + rotation;

                P_unique.emplace(cos(angle) * distance,
                                 sin(angle) * distance);
            }

            copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, perturbationValue);
        }

        void insideDisc(const index_t n, const double radiusOfDisk, vector<Point> &P) {
            typedef CGAL::Random_points_in_disc_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            Point_generator g(radiusOfDisk, m_randCgal);

            unordered_set<Point> P_unique;
            size_t remaining;
            while((remaining = n - P_unique.size()) > 0)
                std::copy_n(g, n, inserter(P_unique, P_unique.end()));

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void onConvexHullInDisc(const index_t n,
                                      const double radius,
                                      vector<Point> &P) {
            typedef CGAL::Random_points_in_disc_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
            unordered_set<Point> P_unique;
            size_t remaining;

            while((remaining = n - P_unique.size()) > 0)
                random_convex_set_2(remaining,inserter(P_unique, P_unique.end()),Point_generator(radius, m_randCgal));

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }


        void insideSquareNormal(const index_t pointsInACuster,
                                               const index_t numberOfClusters,
                                               vector<Point> &P,
                                               const number_t xStdDev = 2.0,
                                               const number_t yStdDev = 2.0) {
            std::mt19937 rngX(seed());
            std::mt19937 rngY(seed());
            std::default_random_engine generatorX(rngX()), generatorY(rngY());
            std::normal_distribution<double> distributionX(0.0, xStdDev), distributionY(2.0, yStdDev);

            std::mt19937 rngShift(std::random_device{}());

            std::uniform_int_distribution shiftDistribution(0, INT32_MAX);

            index_t shiftX, shiftY;
            unordered_set<pair<index_t, index_t>, boost::hash<pair<index_t, index_t>>> S;

            unordered_set<Point> P_unique;

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
                    P_unique.emplace(x, y);
                }
            }

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void contiguousOnGrid(const index_t n, vector<Point> &P) {
            points_on_square_grid_2(ceil(std::sqrt(n)), n, std::back_inserter(P), CGAL::Creator_uniform_2<number_t, Point>());
            perturb(P, m_perturbationValue);
        }

        void randomOnGrid(const index_t n, vector<Point> &P) {
            unordered_set<pair<int, int>, boost::hash<pair<int, int>>> S;
            unordered_set<Point> P_unique;

            std::mt19937 rngX(seed());
            std::mt19937 rngY(seed());
            std::uniform_int_distribution xDistribution(0, (int) ceil(0.7 * n)), yDistribution(0, (int) ceil(0.7 * n));

            index_t count = 0;

            while (count < n) {
                int x = xDistribution(rngX), y = yDistribution(rngY);

                if (S.find(make_pair(x, y)) == S.end()) {
                    P_unique.emplace(x, y);
                    S.insert(make_pair(x, y));
                    count++;
                }
            }

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }

        void insideAnnulus(const index_t n, const double r2, const double r1, vector<Point> &P) {
            assert(r2 > r1);
            unordered_set<Point> P_unique;

            std::default_random_engine generator(seed());
            std::uniform_real_distribution<double> distributionR(r1, r2), distributionT(0, 1);

            for (index_t i = 0; i < n; i++) {
                double t = 2 * M_PI * distributionT(generator);
                double r = distributionR(generator);
                P_unique.emplace(r * cos(t), r * sin(t));
            }

            std::copy(P_unique.begin(),P_unique.end(),back_inserter(P));
            perturb(P, m_perturbationValue);
        }

    private:
        CGAL::Random m_randCgal;
        inline static number_t m_perturbationValue = 0.0001;

        size_t seed() {
            return std::rand();
        }
        double randFloat() {
            return static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
        }
        template<class Container>
        void perturb(Container &P, number_t val) {
            perturb_points_2(P.begin(), P.end(), val, val,m_randCgal);
        }
    };

} // spanners

#endif //SPANNERS_POINTGENERATORS_H
