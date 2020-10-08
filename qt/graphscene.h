#ifndef GRAPHSCENE_H
#define GRAPHSCENE_H

#include <vector>

#include <boost/functional/hash.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // CGAL
#include <CGAL/Delaunay_triangulation_2.h>                      // Triangulations
#include <QAbstractGraphicsShapeItem>
#include <QGraphicsScene>
#include <QPainter>
#include <QPen>
#include <QWidget>

#include "../include/utilities.h"



class GraphScene : public QGraphicsScene
{
    Q_OBJECT
public:

    GraphScene( QObject *parent = nullptr );
    ~GraphScene();

    void changeItemSize( double size );


public slots:
    void advance( int position );

signals:
    void send_playRequest( bool play );
    void send_lengthInfo( int length );

protected:

private:

    std::list< Point > _P;
   // Delaunay_triangulation_2 _DT;

    std::unordered_map<Vertex_handle, QGraphicsEllipseItem*> _Vitems;
    std::unordered_map<Vertex_handle, std::map<Vertex_handle, QGraphicsLineItem*> > _Eitems;


    std::vector <std::pair< Vertex_handle, QAbstractGraphicsShapeItem* > > _focus;

    std::list<TV_Event*>::iterator _pos;
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



#endif // GSNUNF_GRAPHSCENE_H



