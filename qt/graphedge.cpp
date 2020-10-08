#include "graphedge.h"

#include <QPainter>
#include <QStyleOptionGraphicsItem>

GraphEdge::GraphEdge(QGraphicsItem *parent)
    : GraphItem(parent)
{
    _activeColor = Qt::cyan;
    _inactiveColor = Qt::lightGray;
    _active = false;
}

void GraphEdge::paint( QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget )
{
    this->scene();
    painter->setPen( QPen( Qt::green, 3, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin ) );
}

QRectF GraphEdge::boundingRect() const {

}

QPainterPath GraphEdge::shape() const {

}

void GraphEdge::toggleState() {
    _active = !_active;
}
