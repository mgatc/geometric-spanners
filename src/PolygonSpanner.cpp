#include "PolygonSpanner.h"

#include <list>



namespace gsnunf {

PolygonSpanner::PolygonSpanner( SpanningGraph &SG ) : Graph( SG ) {



}




void PolygonSpanner::add_cross_delaunay_edges( Vertex_handle a, Vertex_handle b, Vertex_handle c ) {

}

void PolygonSpanner::add_incident_delaunay_edges( Vertex_handle a, Vertex_handle b, Vertex_handle c ) {

}

double PolygonSpanner::get_angle( Vertex_handle a, Vertex_handle b, Vertex_handle c ) {

}

std::optional<pair<Vertex_handle,Vertex_handle> > PolygonSpanner::get_next_valid_pair( Vertex_circulator &C ) {

}

} // namespace gsnunf

