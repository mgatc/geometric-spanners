#include "graphview.h"

#include <iostream>
#include <QWidget>

GraphView::GraphView(GraphScene *scene, QWidget *parent) : QGraphicsView(scene,parent) {
    connect( &this->_timer, &QTimer::timeout, this, &GraphView::timeoutAdvance );
}

GraphView::GraphView(QWidget *parent) : QGraphicsView(parent) {

}

void GraphView::scaleItemSize( double factor ) {
    QRectF sceneBound = sceneRect();

    double width = abs( sceneBound.right() - sceneBound.left() ),
          height = abs( sceneBound.top() - sceneBound.bottom() );

    // Make sure we aren't going to divide by zero (can happen if point set is empty)
    if( width < 0.00001 || height < 0.00001 )
        return;

    double xFactor = (  _planeWidth ) /  width,
           yFactor = ( _planeHeight ) / height;

    factor = std::min( xFactor, yFactor ) * factor;

    // change vertex size
    static_cast<GraphScene*>(scene())->changeItemSize( factor );
}

void GraphView::timeoutAdvance() {
    emit signalAdvance( ++_posInt );
    emit send_positionInfo( _posInt );
}

void GraphView::setPosition( int position ) {
    _posInt = position;
    emit signalAdvance( position );
}

void GraphView::rcv_playRequest( bool play ) {
    if( play ) {
        // don't reset position until play button pushed again after reaching end
        if( _posInt >= _length )
            setPosition(0);
        _timer.start( _frameRate );
    } else {
        _timer.stop();
    }
    this->_playing = play;
    emit send_playingInfo( play );
    std::cout<<"play"<<"\n";
}

void GraphView::rcv_stepRequest( bool forward ) {
    emit signalAdvance( forward ? ++_posInt : --_posInt );
}

void GraphView::rcv_lengthInfo( int length ) {
    _length = length;
    emit send_lengthInfo( length );
    //std::cout<<"view rcv_lengthInfo:"<<length<<"\n";
}
