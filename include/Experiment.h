#ifndef GEOMETRIC_SPANNERS_EXPERIMENT_H
#define GEOMETRIC_SPANNERS_EXPERIMENT_H

#include <iostream>
#include <list>
#include <optional>
#include <utility>

#include "LineGraphPrinter.h"
#include "metrics.h"
#include "names.h"
#include "utilities.h"

namespace unf_planespanners {

    GraphPrinter dookie;

    GraphPrinter::OptionsList activeEdgeOptions = { // active edge options
            {"color",      dookie.activeEdgeColor},
            {"line width", to_string(dookie.activeEdgeWidth)}
    },
    highlightEdgeOptions = { // active edge options
            {"densely dashed",        nullopt},
            {"color",      dookie.worstPathEdgeColor},
            {"line width", to_string(dookie.activeEdgeWidth)}
    },
    highlightVertexOptions = {
            {"diamond",    nullopt},
            {"vertex",     make_optional(to_string(dookie.vertexRadius*1.61))}, // vertex width
            {"color",      make_optional(dookie.worstPathEdgeColor)}, // text color
            {"fill",       make_optional(dookie.worstPathEdgeColor)}, // vertex color
            {"line width", make_optional(to_string(0))} // vertex border (same color as text)
    },
    activeVertexOptions = {
            {"circle",     nullopt},
            {"vertex",     make_optional(to_string(dookie.vertexRadius))}, // vertex width
            {"color",      make_optional(dookie.backgroundColor)}, // text color
            {"fill",       make_optional(dookie.activeVertexColor)}, // vertex color
            {"line width", make_optional(to_string(0))} // vertex border (same color as text)
    },
    borderOptions = {
            {"border",     make_optional(to_string(dookie.vertexRadius))}, // choose shape of vertex
            {"color",      dookie.activeEdgeColor}, // additional border color
            {"line width", to_string(dookie.inactiveEdgeWidth)}, // additional border width
    };
    template<class Container>
    void PlaneSpannerExperiment(const Container &points,
                                const Algorithm algorithm ) {
        using namespace std;

        list<pair<size_t, size_t> > spanner;

        auto pointsBegin = points.begin(),
                pointsEnd = points.end();

        Timer tim;

        switch (algorithm) {
            case Algorithm::Bgs2005:
                BGS2005(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Lw2004:
                LW2004(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bsx2009:
                BSX2009(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Kpx2010:
                KPX2010(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Kx2012:
                KX2012(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bcc2012_7:
                BCC2012<7>(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bcc2012_6:
                BCC2012<6>(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bhs2017:
                BHS2017(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bghp2010:
                BGHP2010(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Kpt2017:
                KPT2017(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Bkpx2015:
                BKPX2015(pointsBegin, pointsEnd, back_inserter(spanner));
                break;
            case Algorithm::Last:
                assert(false);
        }

        size_t runtime = tim.stop();
        size_t n = points.size();
        size_t deg = degree( spanner.begin(), spanner.end() );
        number_t degAvg = degreeAvg( spanner.begin(), spanner.end() );
        number_t lightness = getLightness( points.begin(), points.end(), spanner.begin(), spanner.end() );
        number_t t = StretchFactorDijkstraReduction( points.begin(), points.end(), spanner.begin(), spanner.end() );
//        number_t t = StretchFactorUsingHeuristic( points.begin(), points.end(), spanner.begin(), spanner.end() );

        list<Edge> WorstPath;
        SFWorstPath( points.begin(), points.end(), spanner.begin(), spanner.end(),
                     make_optional(inserter(WorstPath,WorstPath.begin())) );

        Result result( algorithm, n, runtime, deg, degAvg, lightness, t );
        lineGraph.registerResult(result);
        cout << result << endl;




        GraphPrinter printer;
        double documentSizeInCm = 10;

        printer.autoscale( points.begin(), points.end(), documentSizeInCm);

        printer.beginFigure();
        printer.setCaption(result);

        // remove worst path edges from spanner edges
        spanner.erase( remove_if( spanner.begin(), spanner.end(), [&]( const Edge& edge ){
            auto spannerEdge = normalize_pair(edge);
            return find_if(WorstPath.begin(), WorstPath.end(), [&] ( const Edge& wpEdge ) {
                return spannerEdge == normalize_pair(wpEdge);
            }) != WorstPath.end();
        }), spanner.end());

        printer.addLatexComment("spanner edges (except for worst path)");
        printer.drawEdges( spanner.begin(), spanner.end(), points, activeEdgeOptions);

        printer.addLatexComment("worst path edges");
        printer.drawEdges(WorstPath.begin(),WorstPath.end(), points, highlightEdgeOptions);

        printer.addLatexComment("vertices");
        printer.drawVertices( points.begin(), points.end(), activeVertexOptions);

        printer.addLatexComment("worst pair vertices");
        vector<Point> worstPair = { points[WorstPath.front().first], points[WorstPath.back().second] };
        printer.drawVertices( worstPair.begin(), worstPair.end(), highlightVertexOptions);

        printer.endFigure();


        string outputname = string("RPIS-")
                + to_string(n)
                + "-"
                + Names.at(algorithm);
        printer.print(outputname);
    }

    bool singleRun( size_t n, double width, string resultFilename, optional<string> filename, bool forcePrint, bool printLog )
    {
        double size = width/2; // cgal's generators produce width 2x given value

        // SET POINT SET
        vector<Point> points;
        optional<string> generatedFile = nullopt;
        filename = make_optional("data-150_61.237244x61.237244.txt");

        if( filename )
            readPointsFromFile( back_inserter( points ), *filename );
        else
            generatedFile = make_optional( generateRandomPoints( n, size, back_inserter(points) ) );

        GraphPrinter printer;
        double documentSizeInCm = 10;
        printer.autoscale(points.begin(), points.end(), documentSizeInCm);
        printer.clearCaptionFile();

        // draw point set
        printer.beginFigure();
        string caption = string("$N=")
                        + to_string(n)
                        + "$";
        printer.setCaption(caption);

        printer.drawVertices(points.begin(), points.end(), activeVertexOptions);

        printer.endFigure();
        string outputname = string("RPIS-") + to_string(n);
        printer.print(outputname);

        for( int alg=Algorithm::First;
             alg!=Algorithm::Last; ++alg )
        {
            PlaneSpannerExperiment( points, static_cast<Algorithm>(alg) );
        }



//    if( deg > 11 || t > 7 || forcePrint ) {
//
//       string resultFileName = ( filename ? *filename : *generatedFile );
//       // strip file extension
//       const std::string ext(".txt");
//       if ( resultFileName != ext &&
//            resultFileName.size() > ext.size() &&
//            resultFileName.substr(resultFileName.size() - ext.size()) == ext )
//       {
//          // if so then strip them off
//          resultFileName = resultFileName.substr(0, resultFileName.size() - ext.size());
//       }
//       resultFileName += "_result-";
//       resultFileName += ( filename ? "redo" : "orig" );
//
//       cout << "DEGREE ERROR!!! DEGREE:" << deg << "\n"<<endl;
//       cout << *generatedFile <<endl;
//
//       if( generatedFile )
//           singleRun( n, width, resultFileName, *generatedFile, true, true );
//
//       return false;
//    }

        cout<<endl;

        return true;
    }

    bool experiment( size_t trials, size_t n_start, size_t n_end, size_t increment ) {
        const double width = 10;

        size_t invalid = 0;
        size_t trialNum = 1;

        for( size_t trial=0; trial<trials; ++trial ) {
            for( size_t n=n_start; n<=n_end; n+=increment ) {
                trialNum++;
                if( !singleRun( n, width*sqrt(n), "output", nullopt, false, false ) ) {
                    ++invalid;
                    return false;
                }
            }
        }
        //cout<<"\nTesting Complete. "<< invalid << " of "<<(trials*(n_end-n_start))<<" invalid results.\n\n";
        return invalid == 0;
    }





} // unf_planespanners

#endif //GEOMETRIC_SPANNERS_EXPERIMENT_H
