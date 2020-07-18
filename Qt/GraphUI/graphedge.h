#ifndef GRAPHEDGE_H
#define GRAPHEDGE_H

#include <QGraphicsLineItem>

class GraphEdge : public QGraphicsLineItem
{
    Q_OBJECT
public:
    bool _active;
    std::pair<QPointF,QPointF> _e;
    QColor _activeColor;
    QColor _inactiveColor;

    GraphEdge(QGraphicsItem *parent = nullptr);

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
               QWidget *widget) override;
    QRectF boundingRect() const override;
    QPainterPath shape() const override;

public slots:
    void toggleState();
private:

};

#endif // GRAPHEDGE_H
