#ifndef GSNUNF_SPLITVERTEX_H
#define GSNUNF_SPLITVERTEX_H

#include <memory>

namespace gsnunf {

using namespace std;

struct SplitVertexSet;
struct SplitVertexHasher;

/* Represents a single split of a split vertex where
 * first = the vertex to be split
 * second = the split vertex's s_1
 */
struct SplitVertex {
    using v_type = Vertex_handle;
    using key_type = size_t;

    SplitVertexSet* V;
    v_type v;
    key_type s_1;
    v_type s_1_handle;
    key_type key;

    SplitVertex() {}
    SplitVertex( const v_type& v )
        : v(v) {}
    SplitVertex( const v_type& v, const SplitVertex& s_1 )
        : v(v), s_1( s_1.key ), s_1_handle(s_1.v) {}
    SplitVertex( const v_type& v, const key_type& s_1 )
        : v(v), s_1( s_1 ) {}
    SplitVertex( const SplitVertex& other )
        : V(other.V), v(other.v), s_1(other.s_1), s_1_handle( other.s_1_handle ), key(other.key) {}

    SplitVertex& operator=( const SplitVertex& other ) { //copy assignment
        if( this != &other ) {
            V = other.V;
            v = other.v;
            key = other.key;
            s_1 = other.s_1;
            s_1_handle = other.s_1_handle;
        }
        return *this;
    }
};

ostream& operator<<( ostream& os, const SplitVertex& v ){
    os << v.v->point() << " (" << v.s_1_handle->point() << ")";
    return os;
}

struct CompareSplitVertex {
    bool operator()( const SplitVertex& a, const SplitVertex& b ) {
        return a.v < b.v || ( a.v == b.v && a.s_1 < b.s_1 );
    }
};

// Compare two split vertices for equality
inline bool operator==( const SplitVertex& lhs,
                        const SplitVertex& rhs ) {
    return lhs.v == rhs.v && lhs.s_1_handle == rhs.s_1_handle;
}

// Compare two split vertices for inequality
inline bool operator!=( const SplitVertex& lhs,
                        const SplitVertex& rhs ) {
    return !(lhs == rhs);
}

using key_type = typename SplitVertex::key_type;
using IncidentSplitVertexContainer = std::set< key_type >;
using SplitVertexEdgeMap = std::unordered_map< key_type, IncidentSplitVertexContainer >;

struct SplitVertexSet {
    VertexMap< VertexMap< size_t > > index;

    vector< SplitVertex > V;

    size_t insert( const SplitVertex& v ) {
        size_t vector_key;
        SplitVertex v_insert = v;
        if( contains( index, v_insert.v ) && contains( index.at(v_insert.v), v_insert.s_1_handle ) ) {
            vector_key = index[v_insert.v][v_insert.s_1_handle];
            v_insert.key = vector_key;
            V[vector_key] = v_insert;
        } else {
            vector_key = insert_in_container(v_insert);
            insert_in_map(v_insert);
        }
        return vector_key;
    }
    size_t insert_in_container( SplitVertex& v ) {
        size_t vector_key = V.size();
        v.key = vector_key;
        v.V = this;
        V.push_back(v);
        return vector_key;
    }
    void insert_in_map( const SplitVertex& v ) {
        index[v.v].insert_or_assign( v.s_1_handle, v.key );
    }
    SplitVertex at( const size_t i ) const {
        return V.at(i);
    }
    size_t at( const SplitVertex& v ) const {
        return index.at(v.v).at( V.at(v.s_1).v );
    }
    SplitVertex& at( const Vertex_handle v, const Vertex_handle s_1 ) {
        return V.at( index.at(v).at(s_1) );
    }
};

void add_half_edge( SplitVertexEdgeMap& E, const SplitVertex& a, const SplitVertex& b ) {
    E[a.key].emplace(b.key);
}

size_t add_vertex( SplitVertexSet& V, const Vertex_handle& v, const SplitVertex s_1 ) {
    SplitVertex v_split( v, s_1 );
    //cout<<v_split<<"\n";
    return V.insert(v_split);
}

void add_edge( const DelaunayGraph& SG, SplitVertexEdgeMap& E, const SplitVertex& a, const SplitVertex& b ) {
    assert( SG._DT.is_edge( a.v, b.v ) );

    add_half_edge( E, b, a );
    add_half_edge( E, a, b );
}

// Print the contents of a split vertex edge list
void print( const SplitVertexSet& V, const SplitVertexEdgeMap& E ) {
    for( auto v1 : E ) {
        key_type k = v1.first;
        SplitVertex u = V.at(k);
        cout<< "u: "<<u.v->point()<<" s1: "<<V.at(u.s_1).v->point()<<"\n";
        for( auto v2 : v1.second )
            cout<<"  v: "<<V.at(v2).v->point()<<" s1: "<<V.at(V.at(v2).s_1).v->point()<<"\n";
        cout<<"\n";

    }
}

// Print the contents of a split vertex list
void print( const SplitVertexSet& V ) {
    for( auto v : V.V ) {
        cout<<"v"<<v.key<<": "<<v;
        cout<<" s1_key: "<<v.s_1<<"\n";
    }
}

struct SplitVertexHasher {
    std::size_t operator()( const SplitVertex& k ) const {
        // Compute individual hash values for first, second and third
            // http://stackoverflow.com/a/1646913/126995
        size_t res = 17;
         res = res * 31 + hash< Vertex_handle >()( k.v );
         res = res * 31 + hash< Vertex_handle >()( k.V->at(k.s_1).v );
        return res;
    }
};

} // namespace gsnunf

#endif // GSNUNF_SPLITVERTEX_H
