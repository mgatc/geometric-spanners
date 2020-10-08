#ifndef GRAPHITEM_H
#define GRAPHITEM_H

class GraphItem
{
public:
    bool _active;

    GraphItem();


public slots:
    void toggleState();

protected:
    QColor _activeColor;
    QColor _inactiveColor;

};

#endif // GRAPHITEM_H
