#include "graphscene.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <QGraphicsItem>
#include <QObject>
#include <QTimer>
#include <QWidget>

#include "GraphAlgorithmEvent.h"
#include "mainwindow.h"
#include "SpanningGraph.h"
#include "TransformPolygon.h"
#include "PolygonSpanner.h"
#include "PlanarSpanner.h"

GraphScene::GraphScene( QObject *parent )
    : QGraphicsScene(parent),
      _vertexBrush( _vertexColor[1], Qt::SolidPattern ),
      _focusBrush( _focusColor, Qt::SolidPattern ),
      _defaultPen( _vertexColor[0], 0, Qt::SolidLine, Qt::SquareCap, Qt::BevelJoin ),
      _focusPen( _focusColor, 0, Qt::SolidLine, Qt::SquareCap, Qt::BevelJoin ),
      _edgePen( _edgeColor[0], 0, Qt::SolidLine, Qt::SquareCap, Qt::BevelJoin ) {}

GraphScene::~GraphScene() {
    //delete _DG;
}

void GraphScene::advance( int position ) {

    //position = std::max( 0, std::min( position, _length ) );
    int step = position - _posInt;
    int direction = int( (0<step)-(step<0) ); // get sign of step
    //bool forward = direction > 0;
    auto eventListEnd = _DG->eventListEnd();

    //std::cout<<"advanceto"<<position<<"\n";

    /* _pos (position iterator) and _posInt (position integer) are the iteration variables.
     * The first condition allows for positive or negative iteration. The others ensure
     * the variables stay in bounds.
     */
    while( 0<(position-_posInt)*direction && 0<=_posInt && _posInt<_length ) {

        switch( (*_pos)->type ) {

        case EventType::Vertex: {
            // cast the event to the proper type, _pos is an itr to a ptr
            auto event = static_cast<VertexStatusEvent<DelaunayGraph>* >(*_pos);
            // change state
            bool state = gsnunf::toggle( _Vstate, event->v );
            // change brush color, find the item, and paint it
            _vertexBrush.setColor( _vertexColor[state] );
            auto item = _Vitems.find( event->v )->second;
            item->setBrush( _vertexBrush );
            item->update();
            //std::cout<<"Event type:VertexStatus point:"<<event->v->point()<<" active:"<<event->active<<"\n";

            break;
        }
        case EventType::Edge: {
            // cast the event to the proper type, _pos is an itr to a ptr
            auto event = static_cast<EdgeStatusEvent<DelaunayGraph>* >(*_pos);
            Vertex_handle v1 = event->e.first,
                          v2 = event->e.second;
            // change state table for v1->v2 and v2->v1
            if( !gsnunf::contains( _Estate, v1 ) )
                _Estate.emplace( v1, DelaunayGraph::VertexSet() );
            DelaunayGraph::VertexSet* incident = &_Estate.find(v1)->second;
            bool state = gsnunf::toggle( *incident, v2 );

            if( !gsnunf::contains( _Estate, v2 ) )
                _Estate.emplace( v2, DelaunayGraph::VertexSet() );
            incident = &_Estate.find(v2)->second;
            gsnunf::toggle( *incident, v1 );

            // find item and paint it
            auto item = _Eitems.find(v1)->second.find(v2)->second;
            _edgePen.setColor( _edgeColor[ state ] );
            item->setPen( _edgePen );
            item->update();
            //std::cout<<"Event type:EdgeStatus point:"<<event->e.first->point()<<" "<<event->e.second->point()<<" active:"<<event->active<<"\n";

            break;
        }
        case EventType::Focus: {
            if( direction > 0 ) { // going forward, focus
                // cast the event to the proper type, _pos is an itr to a ptr
                auto event = static_cast<VertexFocusEvent<DelaunayGraph>* >(*_pos);
                size_t lvl = event->level;
                if( lvl <= _focus.size() ) {
                    // change all from lvl end, according to _Vstate
                    for( size_t i=lvl; i<_focus.size(); i++ ) {
                        _vertexBrush.setColor( _vertexColor[ size_t(gsnunf::contains( _Vstate, _focus.at(i).first )) ] );
                        _focus.at(i).second->setBrush( _vertexBrush );
                        _focus.at(i).second->update();
                    }
                    // remove all from e.nextFocus to end
                    _focus.erase( _focus.begin()+int(lvl), _focus.end() );

                    // change vertex style to focus
                    if( gsnunf::contains( _Vitems, event->v ) ) {
                        auto item = _Vitems.find( event->v )->second;
                        item->setBrush( _focusBrush );
                        item->update();
                        _focus.emplace_back( event->v, item );
                    }
                    // NOTE: if the event->v is not a valid vertex, it will be ignored,
                    // but the preceeding focus at lvl will still be removed
                }
            } else { // going backward, remove all focus

                // change all from lvl end, according to _Vstate
                for( auto item : _focus ) {
                    _vertexBrush.setColor( _vertexColor[ size_t(gsnunf::contains( _Vstate, item.first )) ] );
                    item.second->setBrush( _vertexBrush );
                    item.second->update();
                }
                // remove all from e.nextFocus to end
                _focus.clear();
            }
            //std::cout<<"Event type:VertexFocus point:"<<event->v->point()<<" focus:"<<event->level<<"\n";

            break;
        }
        }

        // Update iteration vars
        if(direction>0) ++_pos;
        else --_pos;
        _posInt += direction;
    }
    if( _pos == eventListEnd || _posInt == _length ) {
        std::cout<<"done\n";
        emit send_playRequest( false );
    }
}

void GraphScene::setPointset( QStringList fileNames ) {
    std::list< Point > points;
    std::string delimiter = " ";

    for( auto fName : fileNames ) {
        std::string line;
        std::ifstream file( fName.toStdString() );
        if( file.is_open() ) {
            double x, y;
            while( file >> x >> y )
                points.emplace_back( Point( x, y ) );

            file.close();
        }
    }

    // Set the list of points
    _P.insert( _P.begin(), points.begin(), points.end() );

    // Set the Delaunay triangulation
    _DT.insert( _P.begin(), _P.end() );

    prepareGraphObjects();
}

void GraphScene::changeItemSize( double size ) {
    QRectF bound;
    QPointF center;
    // set new ellipse properties for vertices
    for( auto v : _Vitems ) {
        bound = v.second->rect();
        center = bound.center();

        bound.setWidth( size );
        bound.setHeight( size );
        bound.moveCenter( center );
        v.second->setRect( bound );
    }
}

void GraphScene::runAlgorithm( GraphAlgorithm alg ) {

    switch(alg) {
    case SpanningGraph:
        gsnunf::SpanningGraph( *_DG );
        break;

    case PolygonSpanner:
        gsnunf::SpanningGraph( *_DG );
        gsnunf::PolygonSpanner( *_DG );
        break;

    case PlanarSpanner:
        gsnunf::PlanarSpanner( *_DG, 2 );
        break;
    }
    _pos = _DG->processEventQueue();
    _length = int(_DG->_eventList.size());

    //std::cout<<"ready, size="<<<<"\n";
    emit send_lengthInfo( _length );
    emit send_playRequest( true ); // send play signal
}

void GraphScene::prepareGraphObjects() {
    // Create new DelaunayGraph
    _DG = new gsnunf::DelaunayGraph<Delaunay_triangulation_2>(_DT);

    // add edges
    for( auto it=_DT.finite_edges_begin(); it!=_DT.finite_edges_end(); ++it ) {
        Vertex_handle v1 = it->first->vertex( (it->second+1)%3 )->handle(),
                      v2 = it->first->vertex( (it->second+2)%3 )->handle();
        Point p1 = v1->point(),
              p2 = v2->point();
        QPointF p1q( p1.x(), p1.y() ),
                p2q( p2.x(), p2.y() );
        QLineF line( p1q, p2q );

        QGraphicsLineItem *e = addLine( line, _edgePen );

        std::unordered_map<Vertex_handle, std::map<Vertex_handle, QGraphicsLineItem*> >::iterator adjacent;

        std::tie(adjacent, std::ignore) = _Eitems.emplace( v1, std::map<Vertex_handle,QGraphicsLineItem*>() );
        adjacent->second.insert( { v2, e } );

        std::tie(adjacent, std::ignore) = _Eitems.emplace( v2, std::map<Vertex_handle,QGraphicsLineItem*>() );
        adjacent->second.insert( { v1, e } );
    }
    // add vertices
    for( auto it=_DT.finite_vertices_begin(); it!=_DT.finite_vertices_end(); ++it ) {
        Vertex_handle v1 = it;
        Point p1 = v1->point();
        QPointF center( p1.x(), p1.y() );
        QPointF offset( _defaultPen.width()/2, _defaultPen.width()/2 );
        QRectF bound( center-offset, center+offset );

        QGraphicsEllipseItem *v = addEllipse( bound, _defaultPen, _vertexBrush );
        _Vitems.insert( {v1, v} );
        _Vstate.insert(v1);
    }
}
