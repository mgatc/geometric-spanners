#ifndef GRAPHFOCUS_H
#define GRAPHFOCUS_H

#include "graphitem.h"
#include "graphvertex.h"

class GraphFocus : public QGraphicsItem
{
public:
    int _focusLevel = -1;
    GraphVertex *_vertex;

    GraphFocus(QGraphicsItem *parent);

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
               QWidget *widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;
};

#endif // GRAPHFOCUS_H
