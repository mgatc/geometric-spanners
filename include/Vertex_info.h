#ifndef VERTEX_INFO_H
#define VERTEX_INFO_H

namespace gsnunf {

/**
 *  The object that will be stored in the Triangulation along with each vertex
 */

struct Vertex_info { // accessed via Vertex_handle_instance.info()
    bool is_removed = false;      // Vertex_handle_instance.info().is_removed
    bool on_outer_face = false;   // Vertex_handle_instance.info().on_outer_face
    int incident_chords = 0;      // etc.
};

//v.info() = {true,true,2};       // how to set entire info() struct
//v.info().is_removed = true;     // how to set a single member only

}

#endif
