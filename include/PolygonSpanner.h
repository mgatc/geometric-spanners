#ifndef GSNUNF_POLYGONSPANNER_H
#define GSNUNF_POLYGONSPANNER_H

#include "CGALComponents.h"
#include "Graph.h"

#include <utility>



namespace gsnunf {

class PolygonSpanner : public Graph {
public:

    PolygonSpanner() {}
    ~PolygonSpanner() {}

    PolygonSpanner( SpanningGraph &SG );

protected:
    void add_cross_delaunay_edges( Vertex_handle a, Vertex_handle b, Vertex_handle c );
    void add_incident_delaunay_edges( Vertex_handle a, Vertex_handle b, Vertex_handle c );
    double get_angle( Vertex_handle a, Vertex_handle b, Vertex_handle c );
    std::optional<std::pair<Vertex_handle,Vertex_handle> >get_next_valid_pair( Vertex_circulator &C );

private:


}; // class PolygonSpanner

} // namespace gsnunf

#endif // GSNUNF_POLYGONSPANNER_H
