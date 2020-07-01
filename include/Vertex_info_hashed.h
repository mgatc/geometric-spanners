#ifndef GSNUNF_VERTEX_INFO_H
#define GSNUNF_VERTEX_INFO_H

namespace gsnunf {

typedef std::unordered_set<Vertex_handle> Vertex_set;

template <typename T>
using Vertex_info = unordered_set<Vertex_handle, T>;

bool contains( const Vertex_set& V, const Vertex_handle& v ) {
    return V.find(v) != V.end();
}

} // namespace gsnunf

#endif // GSNUNF_VERTEX_INFO_H
