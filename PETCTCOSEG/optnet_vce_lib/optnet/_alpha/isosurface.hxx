/*
 ==========================================================================
 |   
 |   $Id: isosurface.hxx 186 2005-02-27 02:04:22Z kangli $
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

#ifndef ___ISOSURFACE_HXX___
#   define ___ISOSURFACE_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_base/point3.hxx>
#   include <fstream>
#   include <vector>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class isosurface
///  @brief Isosurfacing using marching cubes.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty,
          typename _Real = double,
          typename _Tg = net_f_xy>
class isosurface
{
public:

    typedef _Ty                                 value_type;
    typedef value_type&                         reference;
    typedef const value_type&                   const_reference;

    typedef _Real                               real_value_type;
    typedef real_value_type&                    real_reference;
    typedef const real_value_type&              real_const_reference;

    typedef array_base<value_type, _Tg>         array_base_type;
    typedef array_ref<value_type, _Tg>          array_ref_type;
    typedef array<value_type, _Tg>              array_type;

    typedef size_t                              size_type;

    struct edge_type
    {
        size_type   v[2];
    };

    typedef std::vector<edge_type>              edge_vector_type;
    typedef typename
        edge_vector_type::iterator              edge_iterator;
    typedef typename
        edge_vector_type::const_iterator        edge_const_iterator;
    typedef typename
        edge_vector_type::reference             edge_reference;
    typedef typename
        edge_vector_type::const_reference       edge_const_reference;

    struct triangle_type
    {
        size_type   v[3];
    };

    typedef std::vector<triangle_type>          triangle_vector_type;
    typedef typename
        triangle_vector_type::iterator          triangle_iterator;
    typedef typename
        triangle_vector_type::const_iterator    triangle_const_iterator;
    typedef typename
        triangle_vector_type::reference         triangle_reference;
    typedef typename
        triangle_vector_type::const_reference   triangle_const_reference;

    typedef point3<real_value_type>             point_type;

    typedef std::vector<point_type>             point_vector_type;
    typedef typename
        point_vector_type::iterator             point_iterator;
    typedef typename
        point_vector_type::const_iterator       point_const_iterator;
    typedef typename
        point_vector_type::reference            point_reference;
    typedef typename
        point_vector_type::const_reference      point_const_reference;


    ///////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    ///////////////////////////////////////////////////////////////////////
    isosurface() :
        m_isovalue(value_type())
    {
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Extract iso-surface using the marching cubes algorithm.
    ///
    ///  @param volume    The 3-D volume in which the isosurface is to be
    ///                   extracted.
    ///  @param isovalue  The value of the isosurface.
    ///
    ///////////////////////////////////////////////////////////////////////
    void find(array_base_type&  volume,
              const value_type& isovalue = value_type()
              )
    {
        if (volume.size_0() <= 0 ||
            volume.size_1() <= 0 ||
            volume.size_2() <= 0)
            return;

        clear();

        m_isovalue = isovalue;

        for (size_type i2 = 1; i2 + 1 < volume.size_2(); ++i2) {
            for (size_type i1 = 1; i1 + 1 < volume.size_1(); ++i1) {
                for (size_type i0 = 1; i0 + 1 < volume.size_0(); ++i0) {

                    // march over a single cube
                    march_cube(volume, i0, i1, i2);

                } // i0
            } // i1
        } // i2
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Save the detected isosurface into an XML file.
    ///
    ///  @param filename  The file name.
    ///////////////////////////////////////////////////////////////////////
    void save(const char* filename,
              double      scale_x = 1.0,
              double      scale_y = 1.0,
              double      scale_z = 1.0
              );

    ///////////////////////////////////////////////////////////////////////
    ///  Save the detected isosurface into an output stream.
    ///
    ///  @param os  The output stream.
    ///////////////////////////////////////////////////////////////////////
    void save(std::ostream& os,
              double        scale_x = 1.0,
              double        scale_y = 1.0,
              double        scale_z = 1.0
              );

    ///////////////////////////////////////////////////////////////////////
    ///  Clear all stored isosurfaces.
    ///////////////////////////////////////////////////////////////////////
    void clear()
    {
        m_triangles.clear();
        m_vertices.clear();
        m_normals.clear();
        m_edges.clear();
    }

    ///////////////////////////////////////////////////////////////////////
    inline triangle_const_iterator triangles_begin() const
                { return m_triangles.begin(); }
    inline triangle_const_iterator triangles_end() const
                { return m_triangles.end();   }

    ///////////////////////////////////////////////////////////////////////
    inline point_const_reference vertex(size_type i) const
                { return m_vertices[i]; }
    inline point_const_reference normal(size_type i) const
                { return m_normals[i];  }
    inline edge_const_reference  edge  (size_type i) const
                { return m_edges  [i];  }

    ///////////////////////////////////////////////////////////////////////
    inline size_type num_triangles() const  { return m_triangles.size(); }
    inline size_type num_vertices() const   { return m_vertices.size();  }
    inline size_type num_normals() const    { return m_normals.size();   }
    inline size_type num_edges() const      { return m_edges.size();     }


private:

    value_type              m_isovalue;
    point_vector_type       m_vertices, m_normals;
    triangle_vector_type    m_triangles;
    edge_vector_type        m_edges;

    ///////////////////////////////////////////////////////////////////////
    static const int             CUBE_EDGE_TBL[256];
    static const int             TRIANGLES_TBL[256][16];

    static const real_value_type CUBE_EDGE_DIR[12][3];
    static const size_type       CUBE_EDGE_LNK[12][2];
    static const size_type       CUBE_VERT_OFS[ 8][3];
    static const size_type       TETR_EDGE_LNK[ 6][2];
    static const size_type       TETR_CUBE_MAP[ 6][4];

    ///////////////////////////////////////////////////////////////////////
    inline real_value_type get_isosurface_offset(const real_value_type& v1,
                                                 const real_value_type& v2
                                                 )
    {
        real_value_type delta = v2 - v1;

        return (delta == 0) ?
            (real_value_type)(0.5) : (m_isovalue - v1) / delta;
    }

    ///////////////////////////////////////////////////////////////////////
    inline void march_cube(array_base_type& volume,
                           size_type        i0,
                           size_type        i1,
                           size_type        i2
                           )
    {
        edge_type       edge;
        point_type      vertex;
        triangle_type   triangle;
        real_value_type cube_value[8];
        size_type       vertex_id_map[12];
        bool            edge_map[12][12];
        int             i, flag_index, edge_flag;

        flag_index = 0;
        for (i = 0; i < 8; ++i) {
            cube_value[i] = volume(i0 + CUBE_VERT_OFS[i][0],
                                   i1 + CUBE_VERT_OFS[i][1],
                                   i2 + CUBE_VERT_OFS[i][2]
                                   );
            if (cube_value[i] <= m_isovalue) {
                flag_index |= 1 << i;
            } // if
        } // i

        edge_flag = CUBE_EDGE_TBL[flag_index];

        // If the cube is entirely inside or outside of the surface,
        // there will be no intersections.
        if (edge_flag == 0) return;

        memset(edge_map, 0, sizeof(edge_map));
        
        for (i = 0; i < 12; ++i) {
            if (edge_flag & (1 << i)) {
                
                real_value_type offset = get_isosurface_offset(
                        cube_value[CUBE_EDGE_LNK[i][0]],
                        cube_value[CUBE_EDGE_LNK[i][1]]
                    );
                
                vertex.v[0] = i0 + (CUBE_VERT_OFS[CUBE_EDGE_LNK[i][0]][0]
                                       + offset * CUBE_EDGE_DIR[i][0]);
                vertex.v[1] = i1 + (CUBE_VERT_OFS[CUBE_EDGE_LNK[i][0]][1]
                                       + offset * CUBE_EDGE_DIR[i][1]);
                vertex.v[2] = i2 + (CUBE_VERT_OFS[CUBE_EDGE_LNK[i][0]][2]
                                       + offset * CUBE_EDGE_DIR[i][2]);

                m_vertices.push_back(vertex);

                vertex_id_map[i] = m_vertices.size() - 1;

                //
                // FIXME: Get normals -- there may be a bug here.
                //
                point_type normal;
                normal[0] = volume(i0 + 1, i1, i2) - volume(i0 - 1, i1, i2);
                normal[1] = volume(i0, i1 + 1, i2) - volume(i0, i1 - 1, i2);
                normal[2] = volume(i0, i1, i2 + 1) - volume(i0, i1, i2 - 1);
                normal.normalize();

                m_normals.push_back(normal);
            } // if
        } // for i

        for (i = 0; i < 15; i += 3) {

            if (TRIANGLES_TBL[flag_index][i] < 0) break;

            int v0 = TRIANGLES_TBL[flag_index][i    ];
            int v1 = TRIANGLES_TBL[flag_index][i + 1];
            int v2 = TRIANGLES_TBL[flag_index][i + 2];
            int u0 = (int)vertex_id_map[v0];
            int u1 = (int)vertex_id_map[v1];
            int u2 = (int)vertex_id_map[v2];

            if (!edge_map[v0][v1]) {
                edge.v[0] = u0; edge.v[1] = u1;
                m_edges.push_back(edge);
                edge_map[v0][v1] =
                    edge_map[v1][v0] = true;
            }

            if (!edge_map[v0][v2]) {
                edge.v[0] = u0; edge.v[1] = u2;
                m_edges.push_back(edge);
                edge_map[v0][v2] = 
                    edge_map[v2][v0] = true;
            }

            if (!edge_map[v1][v2]) {
                edge.v[0] = u1; edge.v[1] = u2;
                m_edges.push_back(edge);
                edge_map[v1][v2] = 
                    edge_map[v2][v1] = true;
            }

            triangle.v[0] = u0;
            triangle.v[1] = u1;
            triangle.v[2] = u2;

            m_triangles.push_back(triangle);

        } // for i
    }

};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_alpha/isosurface.cxx>
#   endif

#endif // ___ISOSURFACE_HXX___
