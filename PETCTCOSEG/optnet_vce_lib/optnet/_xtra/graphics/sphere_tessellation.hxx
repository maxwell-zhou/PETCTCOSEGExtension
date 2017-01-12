/*
 ==========================================================================
 |   
 |   $Id: sphere_tessellation.hxx 186 2005-02-27 02:04:22Z kangli $
 |
 |   Written by Kang Li <kangl@cmu.edu>
 |   Department of Electrical and Computer Engineering
 |   Carnegie Mellon University
 |   
 ==========================================================================
 |   This file is a part of the OptimalNet library.
 ==========================================================================
 | Copyright (c) 2004-2005 Kang Li <kangl@cmu.edu>. All Rights Reserved.
 | 
 | This software is supplied under the terms of a license agreement or
 | nondisclosure agreement  with the author  and may not be copied  or
 | disclosed except in accordance with the terms of that agreement.
 ==========================================================================
 */
 
#ifndef ___SPHERE_TESSELLATION_HXX___
#   define ___SPHERE_TESSELLATION_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/point3.hxx>


/// @namespace optnet
namespace optnet {
    /// @namespace xtra
    namespace xtra {
        /// @namespace graphics
        namespace graphics {

///////////////////////////////////////////////////////////////////////////
///  @class sphere_tessellation
///  @brief Class that provides the sphere tessellation service.
///////////////////////////////////////////////////////////////////////////
template <typename _Real>
class sphere_tessellation
{
public:

    typedef _Real               value_type;

    typedef ptrdiff_t           difference_type;
    typedef size_t              size_type;

    /// A point in 3-D space.
    typedef point3<_Real>       point_type;
    typedef point_type&         point_reference;
    typedef const point_type&   point_const_reference;

    struct edge_type;
    struct triangle_type;

    ///////////////////////////////////////////////////////////////////////
    struct vertex_type      /// Definition of a vertex type.
    {
        point_type*         p;      /// Vertex coordinates.
        edge_type**         e;      /// Incident edges.
        size_type           m;      /// Number of incident edges.
    };

    ///////////////////////////////////////////////////////////////////////
    struct edge_type        /// Definition of an edge type.
    {
        vertex_type*        v[2];   /// Adjacent vertices.
        triangle_type*      t[2];   /// Adjacent triangles.
    };

    ///////////////////////////////////////////////////////////////////////
    struct triangle_type    /// Definition of a triangle type.
    {
        edge_type*          e[3];   /// Triangle edges.
        vertex_type*        v[3];   /// Triangle vertices.
        triangle_type*      a[3];   /// Neighboring triangles.
    };

    typedef edge_type*             edge_pointer;
    typedef vertex_type*           vertex_pointer;
    typedef triangle_type*         triangle_pointer;

    typedef const edge_type*       edge_const_pointer;
    typedef const vertex_type*     vertex_const_pointer;
    typedef const triangle_type*   triangle_const_pointer;

    typedef edge_type*             edge_iterator;
    typedef vertex_type*           vertex_iterator;
    typedef triangle_type*         triangle_iterator;

    typedef const edge_type*       edge_const_iterator;
    typedef const vertex_type*     vertex_const_iterator;
    typedef const triangle_type*   triangle_const_iterator;

    typedef edge_type&             edge_reference;
    typedef vertex_type&           vertex_reference;
    typedef triangle_type&         triangle_reference;

    typedef const edge_type&       edge_const_reference;
    typedef const vertex_type&     vertex_const_reference;
    typedef const triangle_type&   triangle_const_reference;

    ///////////////////////////////////////////////////////////////////////
    struct polyhedron_type  /// Definition of a polyhedron.
    {
        point_type          c;

        size_type           m;
        edge_pointer        e;

        size_type           n;
        vertex_pointer      v;

        size_type           r;
        triangle_pointer    t;
    };


    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    sphere_tessellation();

    ///////////////////////////////////////////////////////////////////////
    /// Create initial polyhedron.
    ///
    /// @param shape Type of initial polyhedron.
    ///
    ///////////////////////////////////////////////////////////////////////
    void initialize(int shape = 8);

    ///////////////////////////////////////////////////////////////////////
    /// Perform sphere tessellation.
    ///
    /// @param level Level of subdivisions.
    ///
    ///////////////////////////////////////////////////////////////////////
    void tessellate(int level);

    ///////////////////////////////////////////////////////////////////////
    /// Clears the tessellation.
    ///////////////////////////////////////////////////////////////////////
    void clear();

    ///////////////////////////////////////////////////////////////////////
    ///  Computes the offset of the edge in the edge vector.
    inline difference_type edge_offset(edge_const_pointer e) const
            { return edge_offset(m_poly, e); }

    ///  Computes the offset of the vertex in the vertex vector.
    inline difference_type vertex_offset(vertex_const_pointer v) const
            { return vertex_offset(m_poly, v); }

    ///  Computes the offset of the triangle in the triangle vector.
    inline difference_type triangle_offset(triangle_const_pointer t) const
            { return triangle_offset(m_poly, t); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the number of edges.
    inline size_type num_edges()     const  { return m_poly.m; }

    ///  Returns the number of vertices.
    inline size_type num_vertices()  const  { return m_poly.n; }

    ///  Returns the number of triangles.
    inline size_type num_triangles() const  { return m_poly.r; }

    ///  Returns an edge_iterator addressing the first element in
    ///  the edge container.
    inline edge_iterator            edge_begin()
                                    { return m_poly.e;  }
    ///  Returns an edge_const_iterator addressing the first element in
    ///  the edge container.
    inline edge_const_iterator      edge_begin() const
                                    { return m_poly.e;  }
    ///  Returns an edge_iterator addressing the location succeeding
    ///  the last element in the edge container.
    inline edge_iterator            edge_end()
                                    { return m_poly.e + m_poly.m;   }
    ///  Returns an edge_const_iterator addressing the location succeeding
    ///  the last element in the edge container.
    inline edge_const_iterator      edge_end() const
                                    { return m_poly.e + m_poly.m;   }

    ///  Returns a vertex_iterator addressing the first element in
    ///  the vertex container.
    inline vertex_iterator          vertex_begin()
                                    { return m_poly.v;  }
    ///  Returns a vertex_const_iterator addressing the first element in
    ///  the vertex container.
    inline vertex_const_iterator    vertex_begin() const
                                    { return m_poly.v;  }
    ///  Returns a vertex_iterator addressing the location succeeding
    ///  the last element in the vertex container.
    inline vertex_iterator          vertex_end()
                                    { return m_poly.v + m_poly.n;   }
    ///  Returns a vertex_const_iterator addressing the location succeeding
    ///  the last element in the vertex container.
    inline vertex_const_iterator    vertex_end() const
                                    { return m_poly.v + m_poly.n;   }

    ///  Returns a triangle_iterator addressing the first element in
    ///  the triangle container.
    inline triangle_iterator        triangle_begin()
                                    { return m_poly.t;  }
    ///  Returns a triangle_const_iterator addressing the first element in
    ///  the triangle container.
    inline triangle_const_iterator  triangle_begin() const
                                    { return m_poly.t;  }
    ///  Returns a triangle_iterator addressing the location succeeding
    ///  the last element in the triangle container.
    inline triangle_iterator        triangle_end()
                                    { return m_poly.t + m_poly.r;   }
    ///  Returns a triangle_const_iterator addressing the location succeeding
    ///  the last element in the triangle container.
    inline triangle_const_iterator  triangle_end() const
                                    { return m_poly.t + m_poly.r;   }

    ///  Returns the centoid of the current polyhedron.
    inline point_reference centroid()               { return m_poly.c; }

    ///  Returns the centoid of the current polyhedron.
    inline point_const_reference centroid() const   { return m_poly.c; }
    

protected:
    
    ///////////////////////////////////////////////////////////////////////
    inline difference_type edge_offset(const polyhedron_type& poly,
                                       edge_const_pointer     e
                                       ) const
    {
        return (difference_type)(e - poly.e);
    }

    ///////////////////////////////////////////////////////////////////////
    inline difference_type vertex_offset(const polyhedron_type& poly,
                                         vertex_const_pointer   v
                                         ) const
    {
        return (difference_type)(v - poly.v);
    }

    ///////////////////////////////////////////////////////////////////////
    inline difference_type triangle_offset(const polyhedron_type& poly,
                                           triangle_const_pointer t
                                           ) const
    {
        return (difference_type)(t - poly.t);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Stores the tessellation result.
    ///////////////////////////////////////////////////////////////////////
    polyhedron_type m_poly;

    ///////////////////////////////////////////////////////////////////////
    int  adjacent_orient(triangle_const_pointer t,
                         triangle_const_pointer ta);
    void update_centroid(size_type n);
    void expand(int level);

};

        } // namespace
    } // namespace
} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_xtra/graphics/sphere_tessellation.cxx>
#   endif

#endif // ___SPHERE_TESSELLATION_CXX___
