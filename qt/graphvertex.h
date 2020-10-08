#ifndef GRAPHVERTEX_H
#define GRAPHVERTEX_H

#include "graphitem.h"

#include "DelaunayGraph.h"

class GraphVertex : public GraphItem
{
    Q_OBJECT
public:
    bool _active;
    QPointF _p;
    QColor _activeColor;
    QColor _inactiveColor;

    GraphVertex(QGraphicsItem *parent = nullptr);

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
               QWidget *widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

public slots:
    void toggleState();
};

#endif // GRAPHVERTEX_H
