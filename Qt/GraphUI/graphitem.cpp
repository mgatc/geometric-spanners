#include "graphitem.h"

GraphItem::GraphItem(QGraphicsItem *parent) : QGraphicsItem(parent)
{

}


void GraphItem::toggleState() {
    _active = !_active;
}
