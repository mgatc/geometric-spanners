#ifndef GRAPHVIEW_H
#define GRAPHVIEW_H

#include <QGraphicsView>
#include <QTimer>
#include <QWidget>

#include "graphscene.h"
#include "../include/DelaunayGraph.h"
#include "../include/GraphAlgoTV.h"



class GraphView : public QGraphicsView
{
    Q_OBJECT

public:
    QTimer _timer;

    GraphView(GraphScene *scene, QWidget *parent = nullptr);
    GraphView(QWidget *parent = nullptr);

    void scaleItemSize( double factor = 1.0 );

public slots:
    void rcv_playRequest( bool play );
    void rcv_stepRequest( bool forward );
    void setPosition( int position ); // connect to rcv_scrubRequest of GraphView
    void rcv_lengthInfo( int length );
//    void setZoom( double level ); // connect to setZoom of GraphView
//    void setFramerate( int fps ); // connect to setFramerate of GraphView

signals:
    void signalAdvance( int step );
    void send_playingInfo( bool playing );
    void send_lengthInfo( int length );
    void send_positionInfo( int position );
    void send_fpsInfo( int fps );

private slots:
    void timeoutAdvance();

private:
    int _posInt = 0;
    int _length;
    bool _playing;

    unsigned int _planeWidth = 600;
    unsigned int _planeHeight = 300;
    unsigned int _vertexSize = 15;
    unsigned int _scale = 800;
    double _size;
    int _frameRate = 10;// 1000.0 / 30;

    QAction *playAct;

};

#endif // GRAPHVIEW_H
