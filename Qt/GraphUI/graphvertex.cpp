#include "graphvertex.h"

#include <QPainter>
#include <QStyleOptionGraphicsItem>

GraphVertex::GraphVertex(QGraphicsItem *parent)
    : GraphItem(parent)
{
    _activeColor = Qt::magenta;
    _inactiveColor = Qt::lightGray;
    _active = true;
}

void GraphVertex::paint( QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget )
{

}

void GraphVertex::toggleState() {
    _active = !_active;
}
