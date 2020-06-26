#ifndef GSNUNF_TRANSFORMPOLYGON_H
#define GSNUNF_TRANSFORMPOLYGON_H

#include <unordered_map>

#include "CGALComponents.h"
#include "SpanningGraph.h"


namespace gsnunf {

typedef unordered_map<Vertex_handle, short int> VisitsAllowedTable;

VisitsAllowedTable TransformPolygon( const SpanningGraph &SG );

} // namespace gsnunf

#endif // GSNUNF_TRANSFORMPOLYGON_H

