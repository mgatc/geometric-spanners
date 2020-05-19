#ifndef GEOMETRICSPANNERPRINTER_H
#define GEOMETRICSPANNERPRINTER_H

#include "CgalComponents.h"

#include<vector>

using namespace std;

class GeometricSpannerPrinter
{
    private:
        vector<Point> pointSet;
        list<Point> centersOfPlacedDisks;
        double r = 0.0;
        string fileName;

    public:
        GeometricSpannerPrinter(vector<Point> &P, list<Point> &centersOfPlacedDisks, double radiusOfDisks, string outputFileName);
        void displayPDF();
};

#endif
