#ifndef GRAPHSCENE_H
#define GRAPHSCENE_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // CGAL
#include <CGAL/Delaunay_triangulation_2.h>                      // Triangulations
#include <QAbstractGraphicsShapeItem>
#include <QGraphicsScene>
#include <QPainter>
#include <QPen>
#include <QWidget>

#include "tools/DelaunayGraph.h"#include "GraphAlgorithmEvent.h"



class GraphScene : public QGraphicsScene
{
    Q_OBJECT
public:
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_2 Point;

    typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation_2;
    typedef gsnunf::DelaunayGraph<Delaunay_triangulation_2> DelaunayGraph;
    typedef DelaunayGraph::Vertex_handle Vertex_handle;

    enum GraphAlgorithm {
        SpanningGraph, PolygonSpanner, PlanarSpanner
    };


    GraphScene( QObject *parent = nullptr );
    ~GraphScene();

    void runAlgorithm( GraphAlgorithm alg );
    void setPointset( QStringList fileNames );
    void changeItemSize( double size );

    DelaunayGraph *_DG;

public slots:
    void advance( int position );

signals:
    void send_playRequest( bool play );
    void send_lengthInfo( int length );

protected:

private:

    std::list< Point > _P;
    Delaunay_triangulation_2 _DT;

    std::unordered_map<Vertex_handle, QGraphicsEllipseItem*> _Vitems;
    std::unordered_map<Vertex_handle, std::map<Vertex_handle, QGraphicsLineItem*> > _Eitems;

    DelaunayGraph::VertexHash _Vstate;
    DelaunayGraph::AdjacencyList _Estate;

    std::vector <std::pair< Vertex_handle, QAbstractGraphicsShapeItem* > > _focus;

    std::list<GraphAlgorithmEvent*>::iterator _pos;
    int _posInt;
    int _length;
    double _size;
    int _scale = 800;

    QColor _vertexColor[2] = { Qt::lightGray, Qt::red };
    QColor _edgeColor[2] = { Qt::lightGray, Qt::red };
    QColor _focusColor = Qt::green;

    QBrush _vertexBrush;
    QBrush _focusBrush;
    QPen _defaultPen;
    QPen _focusPen;
    QPen _edgePen;

    void prepareGraphObjects();
};


#endif // GRAPHSCENE_H
