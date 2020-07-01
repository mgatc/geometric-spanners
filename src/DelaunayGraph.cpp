#include "DelaunayGraph.h"

#include <list>
#include <unordered_set>
#include <unordered_map>
#include <memory>



namespace gsnunf {

DelaunayGraph::DelaunayGraph( DelaunayGraph& G ) : _DT( G._DT ), _E( G._E ) {

}

DelaunayGraph::DelaunayGraph( DelaunayTriangulation* DT ) : _DT( DT ) {

}

DelaunayGraph::DelaunayGraph( shared_ptr<DelaunayTriangulation> DT ) : _DT( DT ) {

}



void DelaunayGraph::add_edge( Vertex_handle v1, Vertex_handle v2 ) {
    //cout<<"  add_edge( " << v1->point() << " - " << v2->point() << " )\n";
    add_half_edge( v1, v2 );
    add_half_edge( v2, v1 );

}

void DelaunayGraph::add_half_edge( Vertex_handle v1, Vertex_handle v2 ) {

    Adjacency_list::iterator v1_incident_vertices = _E.find(v1);

    if( v1_incident_vertices == _E.end() ) // v1 not found in g
        std::tie( v1_incident_vertices, std::ignore ) = _E.insert( make_pair( v1, Incident_vertices() ) );

    v1_incident_vertices->second.insert(v2);
}

/**
 *  Given a Delaunay Triangulation DT and an output list out, compute the canonical ordering of
 *  the underlying point set.
 */
void DelaunayGraph::canonical_order( std::list<Vertex_handle> &out ) {
    out.clear();

    if( _DT->number_of_vertices() <= 3 ) {
        for( auto v=_DT->finite_vertices_begin(); v!=_DT->finite_vertices_end(); ++v )
            out.push_back( v->handle() );
        return;
    }

    std::unordered_set<Vertex_handle> eligible_vertices, reserved_vertices, empty_set;
    Vertex_handle v_1, v_2, v_k;

    Vertex_circulator v_convex_hull = _DT->incident_vertices( _DT->infinite_vertex() ), // create a circulator of the convex hull
                      done( v_convex_hull );

    do { // cycle through convex hull to set vertex info
        v_convex_hull->info().on_outer_face = true;  // incident_chords (called in next loop) relies on accurate on_outer_face values
    } while( ++v_convex_hull != done );

    // v_convex_hull is at the first vertex, v_1
    // v_1 and v_2 must be reserved as the last two vertices to be added

    v_1 = (v_convex_hull++)->handle(); // capture the handle, then increment
    v_2 = (v_convex_hull++)->handle(); // same

    reserved_vertices.insert( v_1 );   // put v_1 and v_2 into a container so we can pass our
    reserved_vertices.insert( v_2 );   // reserved vertices to functions easily

    do { // cycle through convex hull again to add eligible vertices to list
        update_incident_chords( v_convex_hull, eligible_vertices, reserved_vertices, _DT->number_of_vertices(), false );
        std::unordered_set<Vertex_handle> empty_set;
        update_eligible_vertices( eligible_vertices, v_convex_hull, empty_set );
    } while( ++v_convex_hull != done );

    while( !eligible_vertices.empty() ) {
        v_k = *eligible_vertices.begin();
        eligible_vertices.erase( v_k );

        v_k->info().is_removed = true;
        v_k->info().on_outer_face = false;
        out.push_front( v_k );

        // Update neighbors
        Vertex_circulator v_n = _DT->incident_vertices( v_k ),
                          done( v_n );
        do { // update outer face
            if( !_DT->is_infinite(v_n) && !v_n->info().is_removed )
                v_n->info().on_outer_face = true;
        } while( ++v_n != done );

        do { // test incidence of new outer face vertices
            if( !_DT->is_infinite(v_n) && !v_n->info().is_removed && v_n != v_1 && v_n != v_2 ) {
                update_incident_chords( v_n, eligible_vertices, reserved_vertices, _DT->number_of_vertices()-out.size(), true );
                update_eligible_vertices( eligible_vertices, v_n, reserved_vertices );
            }
        } while( ++v_n != done );
    }
    out.push_front( v_2 );
    out.push_front( v_1 );

    assert( _DT->number_of_vertices() == out.size() );
}

int DelaunayGraph::count_valid_neighbors( Vertex_circulator C ) {
    Vertex_circulator done(C);
    int k = 0; // count N_i path length k

    do {
        if( !C->info().is_removed && !_DT->is_infinite(C) )
            k++;
    } while( ++C != done );

    return k;
}

double DelaunayGraph::get_angle( const Vertex_handle &p, const Vertex_handle &q, const Vertex_handle &r ) {
    assert( !_DT->is_infinite(p) && !_DT->is_infinite(q) && !_DT->is_infinite(r) );

    Vector2D pq( p->point(), q->point() );
    Vector2D rq( r->point(), q->point() );

    double result = atan2(rq.y(), rq.x()) - atan2(pq.y(), pq.x());

    // atan() returns a value between -PI and PI. From zero ("up"), CCW rotation is negative and CW is positive.
    // Our zero is also "up," but we only want positive values between 0 and 2*PI:

    result *= -1; // First, invert the result. This will associate CW rotation with positive values.
    if( result < EPSILON ) // Then, if the result is less than 0 (or epsilon for floats) add 2*PI.
        result += 2*PI;

    return result;
}

void DelaunayGraph::normalize_circulator( Vertex_circulator &C ) {
    Vertex_circulator done = C;
    // Position circulator so that we are guaranteed to be on the first vertex on the path N_i
    // First, loop until the circulator reaches an invalid vertex or completes a full rotation
    while( !C->info().is_removed && !_DT->is_infinite(C) && ++C != done );// cout<<v_n->point()<<"\n";
    done = C;
    // Loop until the circulator reaches a valid vertex
    while( ( C->info().is_removed || _DT->is_infinite(C) ) && ++C != done );// cout<<v_n->point()<<"\n";
}

void DelaunayGraph::normalize_properties() {
    for( auto it = _DT->finite_vertices_begin(); it!=_DT->finite_vertices_end(); ++it ) {
        it->info().is_removed = false;    // always true
        it->info().on_outer_face = false; // needs verification
    }

    // The vertices incident to the infinite vertex are the convex hull
    Vertex_circulator v_convex_hull = _DT->incident_vertices( _DT->infinite_vertex() ),
                done( v_convex_hull );
    do {
        v_convex_hull->info().on_outer_face = true;
    } while( ++v_convex_hull != done );
}

void DelaunayGraph::remove_edge( Vertex_handle v1, Vertex_handle v2 ) {
    remove_half_edge( v1, v2 );
    remove_half_edge( v2, v1 );
}

void DelaunayGraph::remove_half_edge( Vertex_handle v1, Vertex_handle v2 ) {

    Adjacency_list::iterator v1_incident_vertices = _E.find(v1);

    // if v1 in present in _E, and v2 is present in v1, remove v2 from v1
    if( v1_incident_vertices != _E.end() && v1_incident_vertices->second.find(v2) != v1_incident_vertices->second.end() )
        v1_incident_vertices->second.erase(v2);

}

/**
 *  Given an unordered set of eligible vertices, a vertex handle v, and an unordered set of
 *  reserved vertices, perform any necessary insertion or removal of v from eligible
 *  based on its current presence in eligible, its incident chord status, and its reserved status.
 */
void DelaunayGraph::update_eligible_vertices( std::unordered_set<Vertex_handle> &eligible, Vertex_handle v, std::unordered_set<Vertex_handle> &reserved ) {
    auto reserved_iter = reserved.find(v);
    bool is_reserved = reserved_iter != reserved.end();

    auto eligible_iter = eligible.find(v);
    bool in_eligible = eligible_iter != eligible.end();

    if( v->info().incident_chords == 0 && !in_eligible && !is_reserved )
        eligible.insert(v);

    if( ( v->info().incident_chords > 0 || is_reserved ) && in_eligible )
        eligible.erase(v);
}

/**
 *  Given a Delaunay Triangulation DT, a vertex handle v_update, an unordered set of eligible
 *  vertices, an unordered set of reserved vertices, the number of vertices left in DT,
 *  and a boolean value for controlling recursive behavior, update v_update's incident chord
 *  count. Additionally, if update_neighbors is true, recursive call update_incident_chords on
 *  each neighbor.
 */
void DelaunayGraph::update_incident_chords( Vertex_handle v_update, std::unordered_set<Vertex_handle> &eligible_vertices, std::unordered_set<Vertex_handle> &reserved_vertices, int vertices_remaining, bool update_neighbors ) {

    v_update->info().incident_chords = 0;

    // Conditions that will result in 0 incident chords to v_update
    if( _DT->is_infinite( v_update ) ) return;      // Vertex is the infinite vertex
    if( v_update->info().is_removed ) return;     // Vertex is removed already
    if( !v_update->info().on_outer_face ) return; // If a vertex is not on the outer face, it cannot be part of a chord
    if( vertices_remaining <= 3 ) return;         // Graphs must contain at least 4 vertices to contain a chord

    // Otherwise, v_update is guaranteed to be incident to 2 other vertices with on_outer_face==true.
    // Any additional vertex incident to v_update with on_outer_face==true forms a chord with v_update.
    // These vertices must be updated too.

    v_update->info().incident_chords = -2;

    Vertex_circulator v_neighbor = _DT->incident_vertices( v_update ),
                      done( v_neighbor );
    do {
        if( !_DT->is_infinite( v_neighbor ) && !v_neighbor->info().is_removed && v_neighbor->info().on_outer_face ) {
            if( update_neighbors )
                update_incident_chords( v_neighbor, eligible_vertices, reserved_vertices, vertices_remaining, false );

            v_update->info().incident_chords++;

            update_eligible_vertices( eligible_vertices, v_neighbor, reserved_vertices );
        }
    } while( ++v_neighbor != done );

    assert( v_update->info().incident_chords >= 0 ); // Make sure our assumptions are correct
}

}
