/*
 ==========================================================================
 |   
 |   $Id: sphere_tessellation.cxx 21 2005-01-14 15:52:31Z kangli $
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
 
#ifndef ___SPHERE_TESSELLATION_CXX___
#   define ___SPHERE_TESSELLATION_CXX___

#   include <optnet/_base/memory.hxx>
#   include <optnet/_xtra/graphics/sphere_tessellation.hxx>
#   include <optnet/_xtra/graphics/intersect.hxx>
#   include <cstdlib>

/// @namespace optnet
namespace optnet {
    /// @namespace xtra
    namespace xtra {
        /// @namespace graphics
        namespace graphics {

///////////////////////////////////////////////////////////////////////////
template <typename _Real>
sphere_tessellation<_Real>::sphere_tessellation()
{
    memset(&m_poly, 0, sizeof(polyhedron_type));
}

///////////////////////////////////////////////////////////////////////////
template <typename _Real>
void sphere_tessellation<_Real>::initialize(int shape)
{
    clear();
    
    switch(shape) {
    case 20:
        {
            // ============================================================
            // Create an icosahedron.
            //
            
            const value_type ICOSAHEDRON_X = 
                sqrt(sqrt(5.0) + 1.0) / sqrt(2.0 * sqrt(5.0));
            const value_type ICOSAHEDRON_Y = 
                sqrt(2.0) / sqrt(5.0 + sqrt(5.0));
            const value_type ICOSAHEDRON_Z = 0.0;

            m_poly.e = new edge_type[30];       // 30 edges
            m_poly.m = 30;
            m_poly.v = new vertex_type[12];     // 12 vertices
            m_poly.n = 12;
            m_poly.t = new triangle_type[20];   // 20 triangles
            m_poly.r = 20;

            // ------------------------------------------------------------
            // vertices
            m_poly.v[0].p        = new point_type;
            m_poly.v[0].p->v[0]  = +ICOSAHEDRON_Z;
            m_poly.v[0].p->v[1]  = +ICOSAHEDRON_X;
            m_poly.v[0].p->v[2]  = -ICOSAHEDRON_Y;
            m_poly.v[0].e        = new edge_type*[5];
            m_poly.v[0].e[0]     = &(m_poly.e[0]);
            m_poly.v[0].e[1]     = &(m_poly.e[2]);
            m_poly.v[0].e[2]     = &(m_poly.e[16]);
            m_poly.v[0].e[3]     = &(m_poly.e[17]);
            m_poly.v[0].e[4]     = &(m_poly.e[29]);
            m_poly.v[0].m        = 5;
            
            m_poly.v[1].p        = new point_type;
            m_poly.v[1].p->v[0]  = +ICOSAHEDRON_X;
            m_poly.v[1].p->v[1]  = +ICOSAHEDRON_Y;
            m_poly.v[1].p->v[2]  = +ICOSAHEDRON_Z;
            m_poly.v[1].e        = new edge_type*[5];
            m_poly.v[1].e[0]     = &(m_poly.e[0]);
            m_poly.v[1].e[1]     = &(m_poly.e[1]);
            m_poly.v[1].e[2]     = &(m_poly.e[4]);
            m_poly.v[1].e[3]     = &(m_poly.e[5]);
            m_poly.v[1].e[4]     = &(m_poly.e[20]);
            m_poly.v[1].m        = 5;

            m_poly.v[2].p        = new point_type;
            m_poly.v[2].p->v[0]  = +ICOSAHEDRON_Y;
            m_poly.v[2].p->v[1]  = +ICOSAHEDRON_Z;
            m_poly.v[2].p->v[2]  = -ICOSAHEDRON_X;
            m_poly.v[2].e        = new edge_type*[5];
            m_poly.v[2].e[0]     = &(m_poly.e[1]);
            m_poly.v[2].e[1]     = &(m_poly.e[2]);
            m_poly.v[2].e[2]     = &(m_poly.e[9]);
            m_poly.v[2].e[3]     = &(m_poly.e[10]);
            m_poly.v[2].e[4]     = &(m_poly.e[15]);
            m_poly.v[2].m        = 5;

            m_poly.v[3].p        = new point_type;
            m_poly.v[3].p->v[0]  = +ICOSAHEDRON_Y;
            m_poly.v[3].p->v[1]  = +ICOSAHEDRON_Z;
            m_poly.v[3].p->v[2]  = +ICOSAHEDRON_X;
            m_poly.v[3].e        = new edge_type*[5];
            m_poly.v[3].e[0]     = &(m_poly.e[3]);
            m_poly.v[3].e[1]     = &(m_poly.e[5]);
            m_poly.v[3].e[2]     = &(m_poly.e[7]);
            m_poly.v[3].e[3]     = &(m_poly.e[8]);
            m_poly.v[3].e[4]     = &(m_poly.e[25]);
            m_poly.v[3].m        = 5;

            m_poly.v[4].p        = new point_type;
            m_poly.v[4].p->v[0]  = +ICOSAHEDRON_X;
            m_poly.v[4].p->v[1]  = -ICOSAHEDRON_Y;
            m_poly.v[4].p->v[2]  = +ICOSAHEDRON_Z;
            m_poly.v[4].e        = new edge_type*[5];
            m_poly.v[4].e[0]     = &(m_poly.e[3]);
            m_poly.v[4].e[1]     = &(m_poly.e[4]);
            m_poly.v[4].e[2]     = &(m_poly.e[10]);
            m_poly.v[4].e[3]     = &(m_poly.e[11]);
            m_poly.v[4].e[4]     = &(m_poly.e[28]);
            m_poly.v[4].m        = 5;

            m_poly.v[5].p        = new point_type;
            m_poly.v[5].p->v[0]  = +ICOSAHEDRON_Z;
            m_poly.v[5].p->v[1]  = +ICOSAHEDRON_X;
            m_poly.v[5].p->v[2]  = +ICOSAHEDRON_Y;
            m_poly.v[5].e        = new edge_type*[5];
            m_poly.v[5].e[0]     = &(m_poly.e[6]);
            m_poly.v[5].e[1]     = &(m_poly.e[8]);
            m_poly.v[5].e[2]     = &(m_poly.e[12]);
            m_poly.v[5].e[3]     = &(m_poly.e[16]);
            m_poly.v[5].e[4]     = &(m_poly.e[20]);
            m_poly.v[5].m        = 5;

            m_poly.v[6].p        = new point_type;
            m_poly.v[6].p->v[0]  = -ICOSAHEDRON_Y;
            m_poly.v[6].p->v[1]  = +ICOSAHEDRON_Z;
            m_poly.v[6].p->v[2]  = +ICOSAHEDRON_X;
            m_poly.v[6].e        = new edge_type*[5];
            m_poly.v[6].e[0]     = &(m_poly.e[6]);
            m_poly.v[6].e[1]     = &(m_poly.e[7]);
            m_poly.v[6].e[2]     = &(m_poly.e[13]);
            m_poly.v[6].e[3]     = &(m_poly.e[23]);
            m_poly.v[6].e[4]     = &(m_poly.e[24]);
            m_poly.v[6].m        = 5;

            m_poly.v[7].p        = new point_type;
            m_poly.v[7].p->v[0]  = +ICOSAHEDRON_Z;
            m_poly.v[7].p->v[1]  = -ICOSAHEDRON_X;
            m_poly.v[7].p->v[2]  = -ICOSAHEDRON_Y;
            m_poly.v[7].e        = new edge_type*[5];
            m_poly.v[7].e[0]     = &(m_poly.e[9]);
            m_poly.v[7].e[1]     = &(m_poly.e[11]);
            m_poly.v[7].e[2]     = &(m_poly.e[14]);
            m_poly.v[7].e[3]     = &(m_poly.e[18]);
            m_poly.v[7].e[4]     = &(m_poly.e[22]);
            m_poly.v[7].m        = 5;

            m_poly.v[8].p        = new point_type;
            m_poly.v[8].p->v[0]  = -ICOSAHEDRON_X;
            m_poly.v[8].p->v[1]  = +ICOSAHEDRON_Y;
            m_poly.v[8].p->v[2]  = +ICOSAHEDRON_Z;
            m_poly.v[8].e        = new edge_type*[5];
            m_poly.v[8].e[0]     = &(m_poly.e[12]);
            m_poly.v[8].e[1]     = &(m_poly.e[13]);
            m_poly.v[8].e[2]     = &(m_poly.e[17]);
            m_poly.v[8].e[3]     = &(m_poly.e[26]);
            m_poly.v[8].e[4]     = &(m_poly.e[27]);
            m_poly.v[8].m        = 5;

            m_poly.v[9].p        = new point_type;
            m_poly.v[9].p->v[0]  = -ICOSAHEDRON_Y;
            m_poly.v[9].p->v[1]  = +ICOSAHEDRON_Z;
            m_poly.v[9].p->v[2]  = -ICOSAHEDRON_X;
            m_poly.v[9].e        = new edge_type*[5];
            m_poly.v[9].e[0]     = &(m_poly.e[14]);
            m_poly.v[9].e[1]     = &(m_poly.e[15]);
            m_poly.v[9].e[2]     = &(m_poly.e[19]);
            m_poly.v[9].e[3]     = &(m_poly.e[27]);
            m_poly.v[9].e[4]     = &(m_poly.e[29]);
            m_poly.v[9].m        = 5;

            m_poly.v[10].p       = new point_type;
            m_poly.v[10].p->v[0] = -ICOSAHEDRON_X;
            m_poly.v[10].p->v[1] = -ICOSAHEDRON_Y;
            m_poly.v[10].p->v[2] = +ICOSAHEDRON_Z;
            m_poly.v[10].e       = new edge_type*[5];
            m_poly.v[10].e[0]    = &(m_poly.e[18]);
            m_poly.v[10].e[1]    = &(m_poly.e[19]);
            m_poly.v[10].e[2]    = &(m_poly.e[21]);
            m_poly.v[10].e[3]    = &(m_poly.e[24]);
            m_poly.v[10].e[4]    = &(m_poly.e[26]);
            m_poly.v[10].m       = 5;

            m_poly.v[11].p       = new point_type;
            m_poly.v[11].p->v[0] = +ICOSAHEDRON_Z;
            m_poly.v[11].p->v[1] = -ICOSAHEDRON_X;
            m_poly.v[11].p->v[2] = +ICOSAHEDRON_Y;
            m_poly.v[11].e       = new edge_type*[5];
            m_poly.v[11].e[0]    = &(m_poly.e[21]);
            m_poly.v[11].e[1]    = &(m_poly.e[22]);
            m_poly.v[11].e[2]    = &(m_poly.e[23]);
            m_poly.v[11].e[3]    = &(m_poly.e[25]);
            m_poly.v[11].e[4]    = &(m_poly.e[28]);
            m_poly.v[11].m       = 5;

            // ------------------------------------------------------------
            // edges
            m_poly.e[0].v[0]     = &(m_poly.v[0]);
            m_poly.e[0].v[1]     = &(m_poly.v[1]);
            m_poly.e[0].t[0]     = &(m_poly.t[0]);
            m_poly.e[0].t[1]     = &(m_poly.t[8]);

            m_poly.e[1].v[0]     = &(m_poly.v[2]);
            m_poly.e[1].v[1]     = &(m_poly.v[1]);
            m_poly.e[1].t[0]     = &(m_poly.t[0]);
            m_poly.e[1].t[1]     = &(m_poly.t[17]);

            m_poly.e[2].v[0]     = &(m_poly.v[0]);
            m_poly.e[2].v[1]     = &(m_poly.v[2]);
            m_poly.e[2].t[0]     = &(m_poly.t[0]);
            m_poly.e[2].t[1]     = &(m_poly.t[19]);

            m_poly.e[3].v[0]     = &(m_poly.v[3]);
            m_poly.e[3].v[1]     = &(m_poly.v[4]);
            m_poly.e[3].t[0]     = &(m_poly.t[1]);
            m_poly.e[3].t[1]     = &(m_poly.t[14]);

            m_poly.e[4].v[0]     = &(m_poly.v[1]);
            m_poly.e[4].v[1]     = &(m_poly.v[4]);
            m_poly.e[4].t[0]     = &(m_poly.t[1]);
            m_poly.e[4].t[1]     = &(m_poly.t[17]);

            m_poly.e[5].v[0]     = &(m_poly.v[3]);
            m_poly.e[5].v[1]     = &(m_poly.v[1]);
            m_poly.e[5].t[0]     = &(m_poly.t[1]);
            m_poly.e[5].t[1]     = &(m_poly.t[10]);

            m_poly.e[6].v[0]     = &(m_poly.v[5]);
            m_poly.e[6].v[1]     = &(m_poly.v[6]);
            m_poly.e[6].t[0]     = &(m_poly.t[2]);
            m_poly.e[6].t[1]     = &(m_poly.t[4]);

            m_poly.e[7].v[0]     = &(m_poly.v[3]);
            m_poly.e[7].v[1]     = &(m_poly.v[6]);
            m_poly.e[7].t[0]     = &(m_poly.t[2]);
            m_poly.e[7].t[1]     = &(m_poly.t[12]);

            m_poly.e[8].v[0]     = &(m_poly.v[5]);
            m_poly.e[8].v[1]     = &(m_poly.v[3]);
            m_poly.e[8].t[0]     = &(m_poly.t[2]);
            m_poly.e[8].t[1]     = &(m_poly.t[10]);

            m_poly.e[9].v[0]     = &(m_poly.v[7]);
            m_poly.e[9].v[1]     = &(m_poly.v[2]);
            m_poly.e[9].t[0]     = &(m_poly.t[3]);
            m_poly.e[9].t[1]     = &(m_poly.t[5]);

            m_poly.e[10].v[0]    = &(m_poly.v[2]);
            m_poly.e[10].v[1]    = &(m_poly.v[4]);
            m_poly.e[10].t[0]    = &(m_poly.t[3]);
            m_poly.e[10].t[1]    = &(m_poly.t[17]);

            m_poly.e[11].v[0]    = &(m_poly.v[7]);
            m_poly.e[11].v[1]    = &(m_poly.v[4]);
            m_poly.e[11].t[0]    = &(m_poly.t[3]);
            m_poly.e[11].t[1]    = &(m_poly.t[16]);

            m_poly.e[12].v[0]    = &(m_poly.v[5]);
            m_poly.e[12].v[1]    = &(m_poly.v[8]);
            m_poly.e[12].t[0]    = &(m_poly.t[4]);
            m_poly.e[12].t[1]    = &(m_poly.t[6]);

            m_poly.e[13].v[0]    = &(m_poly.v[6]);
            m_poly.e[13].v[1]    = &(m_poly.v[8]);
            m_poly.e[13].t[0]    = &(m_poly.t[4]);
            m_poly.e[13].t[1]    = &(m_poly.t[15]);

            m_poly.e[14].v[0]    = &(m_poly.v[7]);
            m_poly.e[14].v[1]    = &(m_poly.v[9]);
            m_poly.e[14].t[0]    = &(m_poly.t[5]);
            m_poly.e[14].t[1]    = &(m_poly.t[7]);

            m_poly.e[15].v[0]    = &(m_poly.v[2]);
            m_poly.e[15].v[1]    = &(m_poly.v[9]);
            m_poly.e[15].t[0]    = &(m_poly.t[5]);
            m_poly.e[15].t[1]    = &(m_poly.t[19]);

            m_poly.e[16].v[0]    = &(m_poly.v[5]);
            m_poly.e[16].v[1]    = &(m_poly.v[0]);
            m_poly.e[16].t[0]    = &(m_poly.t[6]);
            m_poly.e[16].t[1]    = &(m_poly.t[8]);

            m_poly.e[17].v[0]    = &(m_poly.v[0]);
            m_poly.e[17].v[1]    = &(m_poly.v[8]);
            m_poly.e[17].t[0]    = &(m_poly.t[6]);
            m_poly.e[17].t[1]    = &(m_poly.t[18]);

            m_poly.e[18].v[0]    = &(m_poly.v[7]);
            m_poly.e[18].v[1]    = &(m_poly.v[10]);
            m_poly.e[18].t[0]    = &(m_poly.t[7]);
            m_poly.e[18].t[1]    = &(m_poly.t[9]);

            m_poly.e[19].v[0]    = &(m_poly.v[9]);
            m_poly.e[19].v[1]    = &(m_poly.v[10]);
            m_poly.e[19].t[0]    = &(m_poly.t[7]);
            m_poly.e[19].t[1]    = &(m_poly.t[13]);

            m_poly.e[20].v[0]    = &(m_poly.v[5]);
            m_poly.e[20].v[1]    = &(m_poly.v[1]);
            m_poly.e[20].t[0]    = &(m_poly.t[8]);
            m_poly.e[20].t[1]    = &(m_poly.t[10]);

            m_poly.e[21].v[0]    = &(m_poly.v[11]);
            m_poly.e[21].v[1]    = &(m_poly.v[10]);
            m_poly.e[21].t[0]    = &(m_poly.t[9]);
            m_poly.e[21].t[1]    = &(m_poly.t[11]);

            m_poly.e[22].v[0]    = &(m_poly.v[11]);
            m_poly.e[22].v[1]    = &(m_poly.v[7]);
            m_poly.e[22].t[0]    = &(m_poly.t[9]);
            m_poly.e[22].t[1]    = &(m_poly.t[16]);

            m_poly.e[23].v[0]    = &(m_poly.v[11]);
            m_poly.e[23].v[1]    = &(m_poly.v[6]);
            m_poly.e[23].t[0]    = &(m_poly.t[11]);
            m_poly.e[23].t[1]    = &(m_poly.t[12]);

            m_poly.e[24].v[0]    = &(m_poly.v[6]);
            m_poly.e[24].v[1]    = &(m_poly.v[10]);
            m_poly.e[24].t[0]    = &(m_poly.t[11]);
            m_poly.e[24].t[1]    = &(m_poly.t[15]);

            m_poly.e[25].v[0]    = &(m_poly.v[11]);
            m_poly.e[25].v[1]    = &(m_poly.v[3]);
            m_poly.e[25].t[0]    = &(m_poly.t[12]);
            m_poly.e[25].t[1]    = &(m_poly.t[14]);

            m_poly.e[26].v[0]    = &(m_poly.v[8]);
            m_poly.e[26].v[1]    = &(m_poly.v[10]);
            m_poly.e[26].t[0]    = &(m_poly.t[13]);
            m_poly.e[26].t[1]    = &(m_poly.t[15]);

            m_poly.e[27].v[0]    = &(m_poly.v[9]);
            m_poly.e[27].v[1]    = &(m_poly.v[8]);
            m_poly.e[27].t[0]    = &(m_poly.t[13]);
            m_poly.e[27].t[1]    = &(m_poly.t[18]);

            m_poly.e[28].v[0]    = &(m_poly.v[11]);
            m_poly.e[28].v[1]    = &(m_poly.v[4]);
            m_poly.e[28].t[0]    = &(m_poly.t[14]);
            m_poly.e[28].t[1]    = &(m_poly.t[16]);

            m_poly.e[29].v[0]    = &(m_poly.v[0]);
            m_poly.e[29].v[1]    = &(m_poly.v[9]);
            m_poly.e[29].t[0]    = &(m_poly.t[18]);
            m_poly.e[29].t[1]    = &(m_poly.t[19]);

            // ------------------------------------------------------------
            // triangles
            m_poly.t[0].e[0]     = &(m_poly.e[0]);
            m_poly.t[0].e[1]     = &(m_poly.e[1]);
            m_poly.t[0].e[2]     = &(m_poly.e[2]);
            m_poly.t[0].v[0]     = &(m_poly.v[0]);
            m_poly.t[0].v[1]     = &(m_poly.v[1]);
            m_poly.t[0].v[2]     = &(m_poly.v[2]);
            m_poly.t[0].a[0]     = &(m_poly.t[8]);
            m_poly.t[0].a[1]     = &(m_poly.t[17]);
            m_poly.t[0].a[2]     = &(m_poly.t[19]);

            m_poly.t[1].e[0]     = &(m_poly.e[3]);
            m_poly.t[1].e[1]     = &(m_poly.e[4]);
            m_poly.t[1].e[2]     = &(m_poly.e[5]);
            m_poly.t[1].v[0]     = &(m_poly.v[3]);
            m_poly.t[1].v[1]     = &(m_poly.v[4]);
            m_poly.t[1].v[2]     = &(m_poly.v[1]);
            m_poly.t[1].a[0]     = &(m_poly.t[14]);
            m_poly.t[1].a[1]     = &(m_poly.t[17]);
            m_poly.t[1].a[2]     = &(m_poly.t[10]);

            m_poly.t[2].e[0]     = &(m_poly.e[6]);
            m_poly.t[2].e[1]     = &(m_poly.e[7]);
            m_poly.t[2].e[2]     = &(m_poly.e[8]);
            m_poly.t[2].v[0]     = &(m_poly.v[5]);
            m_poly.t[2].v[1]     = &(m_poly.v[6]);
            m_poly.t[2].v[2]     = &(m_poly.v[3]);
            m_poly.t[2].a[0]     = &(m_poly.t[4]);
            m_poly.t[2].a[1]     = &(m_poly.t[12]);
            m_poly.t[2].a[2]     = &(m_poly.t[10]);

            m_poly.t[3].e[0]     = &(m_poly.e[9]);
            m_poly.t[3].e[1]     = &(m_poly.e[10]);
            m_poly.t[3].e[2]     = &(m_poly.e[11]);
            m_poly.t[3].v[0]     = &(m_poly.v[7]);
            m_poly.t[3].v[1]     = &(m_poly.v[2]);
            m_poly.t[3].v[2]     = &(m_poly.v[4]);
            m_poly.t[3].a[0]     = &(m_poly.t[5]);
            m_poly.t[3].a[1]     = &(m_poly.t[17]);
            m_poly.t[3].a[2]     = &(m_poly.t[16]);

            m_poly.t[4].e[0]     = &(m_poly.e[12]);
            m_poly.t[4].e[1]     = &(m_poly.e[13]);
            m_poly.t[4].e[2]     = &(m_poly.e[6]);
            m_poly.t[4].v[0]     = &(m_poly.v[5]);
            m_poly.t[4].v[1]     = &(m_poly.v[8]);
            m_poly.t[4].v[2]     = &(m_poly.v[6]);
            m_poly.t[4].a[0]     = &(m_poly.t[6]);
            m_poly.t[4].a[1]     = &(m_poly.t[15]);
            m_poly.t[4].a[2]     = &(m_poly.t[2]);

            m_poly.t[5].e[0]     = &(m_poly.e[14]);
            m_poly.t[5].e[1]     = &(m_poly.e[15]);
            m_poly.t[5].e[2]     = &(m_poly.e[9]);
            m_poly.t[5].v[0]     = &(m_poly.v[7]);
            m_poly.t[5].v[1]     = &(m_poly.v[9]);
            m_poly.t[5].v[2]     = &(m_poly.v[2]);
            m_poly.t[5].a[0]     = &(m_poly.t[7]);
            m_poly.t[5].a[1]     = &(m_poly.t[19]);
            m_poly.t[5].a[2]     = &(m_poly.t[3]);

            m_poly.t[6].e[0]     = &(m_poly.e[16]);
            m_poly.t[6].e[1]     = &(m_poly.e[17]);
            m_poly.t[6].e[2]     = &(m_poly.e[12]);
            m_poly.t[6].v[0]     = &(m_poly.v[5]);
            m_poly.t[6].v[1]     = &(m_poly.v[0]);
            m_poly.t[6].v[2]     = &(m_poly.v[8]);
            m_poly.t[6].a[0]     = &(m_poly.t[8]);
            m_poly.t[6].a[1]     = &(m_poly.t[18]);
            m_poly.t[6].a[2]     = &(m_poly.t[4]);

            m_poly.t[7].e[0]     = &(m_poly.e[18]);
            m_poly.t[7].e[1]     = &(m_poly.e[19]);
            m_poly.t[7].e[2]     = &(m_poly.e[14]);
            m_poly.t[7].v[0]     = &(m_poly.v[7]);
            m_poly.t[7].v[1]     = &(m_poly.v[10]);
            m_poly.t[7].v[2]     = &(m_poly.v[9]);
            m_poly.t[7].a[0]     = &(m_poly.t[9]);
            m_poly.t[7].a[1]     = &(m_poly.t[13]);
            m_poly.t[7].a[2]     = &(m_poly.t[5]);

            m_poly.t[8].e[0]     = &(m_poly.e[20]);
            m_poly.t[8].e[1]     = &(m_poly.e[0]);
            m_poly.t[8].e[2]     = &(m_poly.e[16]);
            m_poly.t[8].v[0]     = &(m_poly.v[5]);
            m_poly.t[8].v[1]     = &(m_poly.v[1]);
            m_poly.t[8].v[2]     = &(m_poly.v[0]);
            m_poly.t[8].a[0]     = &(m_poly.t[10]);
            m_poly.t[8].a[1]     = &(m_poly.t[0]);
            m_poly.t[8].a[2]     = &(m_poly.t[6]);

            m_poly.t[9].e[0]     = &(m_poly.e[21]);
            m_poly.t[9].e[1]     = &(m_poly.e[18]);
            m_poly.t[9].e[2]     = &(m_poly.e[22]);
            m_poly.t[9].v[0]     = &(m_poly.v[11]);
            m_poly.t[9].v[1]     = &(m_poly.v[10]);
            m_poly.t[9].v[2]     = &(m_poly.v[7]);
            m_poly.t[9].a[0]     = &(m_poly.t[11]);
            m_poly.t[9].a[1]     = &(m_poly.t[7]);
            m_poly.t[9].a[2]     = &(m_poly.t[16]);

            m_poly.t[10].e[0]    = &(m_poly.e[8]);
            m_poly.t[10].e[1]    = &(m_poly.e[5]);
            m_poly.t[10].e[2]    = &(m_poly.e[20]);
            m_poly.t[10].v[0]    = &(m_poly.v[5]);
            m_poly.t[10].v[1]    = &(m_poly.v[3]);
            m_poly.t[10].v[2]    = &(m_poly.v[1]);
            m_poly.t[10].a[0]    = &(m_poly.t[2]);
            m_poly.t[10].a[1]    = &(m_poly.t[1]);
            m_poly.t[10].a[2]    = &(m_poly.t[8]);

            m_poly.t[11].e[0]    = &(m_poly.e[23]);
            m_poly.t[11].e[1]    = &(m_poly.e[24]);
            m_poly.t[11].e[2]    = &(m_poly.e[21]);
            m_poly.t[11].v[0]    = &(m_poly.v[11]);
            m_poly.t[11].v[1]    = &(m_poly.v[6]);
            m_poly.t[11].v[2]    = &(m_poly.v[10]);
            m_poly.t[11].a[0]    = &(m_poly.t[12]);
            m_poly.t[11].a[1]    = &(m_poly.t[15]);
            m_poly.t[11].a[2]    = &(m_poly.t[9]);

            m_poly.t[12].e[0]    = &(m_poly.e[25]);
            m_poly.t[12].e[1]    = &(m_poly.e[7]);
            m_poly.t[12].e[2]    = &(m_poly.e[23]);
            m_poly.t[12].v[0]    = &(m_poly.v[11]);
            m_poly.t[12].v[1]    = &(m_poly.v[3]);
            m_poly.t[12].v[2]    = &(m_poly.v[6]);
            m_poly.t[12].a[0]    = &(m_poly.t[14]);
            m_poly.t[12].a[1]    = &(m_poly.t[2]);
            m_poly.t[12].a[2]    = &(m_poly.t[11]);

            m_poly.t[13].e[0]    = &(m_poly.e[19]);
            m_poly.t[13].e[1]    = &(m_poly.e[26]);
            m_poly.t[13].e[2]    = &(m_poly.e[27]);
            m_poly.t[13].v[0]    = &(m_poly.v[9]);
            m_poly.t[13].v[1]    = &(m_poly.v[10]);
            m_poly.t[13].v[2]    = &(m_poly.v[8]);
            m_poly.t[13].a[0]    = &(m_poly.t[7]);
            m_poly.t[13].a[1]    = &(m_poly.t[15]);
            m_poly.t[13].a[2]    = &(m_poly.t[18]);

            m_poly.t[14].e[0]    = &(m_poly.e[28]);
            m_poly.t[14].e[1]    = &(m_poly.e[3]);
            m_poly.t[14].e[2]    = &(m_poly.e[25]);
            m_poly.t[14].v[0]    = &(m_poly.v[11]);
            m_poly.t[14].v[1]    = &(m_poly.v[4]);
            m_poly.t[14].v[2]    = &(m_poly.v[3]);
            m_poly.t[14].a[0]    = &(m_poly.t[16]);
            m_poly.t[14].a[1]    = &(m_poly.t[1]);
            m_poly.t[14].a[2]    = &(m_poly.t[12]);

            m_poly.t[15].e[0]    = &(m_poly.e[13]);
            m_poly.t[15].e[1]    = &(m_poly.e[26]);
            m_poly.t[15].e[2]    = &(m_poly.e[24]);
            m_poly.t[15].v[0]    = &(m_poly.v[6]);
            m_poly.t[15].v[1]    = &(m_poly.v[8]);
            m_poly.t[15].v[2]    = &(m_poly.v[10]);
            m_poly.t[15].a[0]    = &(m_poly.t[4]);
            m_poly.t[15].a[1]    = &(m_poly.t[13]);
            m_poly.t[15].a[2]    = &(m_poly.t[11]);

            m_poly.t[16].e[0]    = &(m_poly.e[22]);
            m_poly.t[16].e[1]    = &(m_poly.e[11]);
            m_poly.t[16].e[2]    = &(m_poly.e[28]);
            m_poly.t[16].v[0]    = &(m_poly.v[11]);
            m_poly.t[16].v[1]    = &(m_poly.v[7]);
            m_poly.t[16].v[2]    = &(m_poly.v[4]);
            m_poly.t[16].a[0]    = &(m_poly.t[9]);
            m_poly.t[16].a[1]    = &(m_poly.t[3]);
            m_poly.t[16].a[2]    = &(m_poly.t[14]);

            m_poly.t[17].e[0]    = &(m_poly.e[1]);
            m_poly.t[17].e[1]    = &(m_poly.e[4]);
            m_poly.t[17].e[2]    = &(m_poly.e[10]);
            m_poly.t[17].v[0]    = &(m_poly.v[2]);
            m_poly.t[17].v[1]    = &(m_poly.v[1]);
            m_poly.t[17].v[2]    = &(m_poly.v[4]);
            m_poly.t[17].a[0]    = &(m_poly.t[0]);
            m_poly.t[17].a[1]    = &(m_poly.t[1]);
            m_poly.t[17].a[2]    = &(m_poly.t[3]);

            m_poly.t[18].e[0]    = &(m_poly.e[29]);
            m_poly.t[18].e[1]    = &(m_poly.e[27]);
            m_poly.t[18].e[2]    = &(m_poly.e[17]);
            m_poly.t[18].v[0]    = &(m_poly.v[0]);
            m_poly.t[18].v[1]    = &(m_poly.v[9]);
            m_poly.t[18].v[2]    = &(m_poly.v[8]);
            m_poly.t[18].a[0]    = &(m_poly.t[19]);
            m_poly.t[18].a[1]    = &(m_poly.t[13]);
            m_poly.t[18].a[2]    = &(m_poly.t[6]);

            m_poly.t[19].e[0]    = &(m_poly.e[2]);
            m_poly.t[19].e[1]    = &(m_poly.e[15]);
            m_poly.t[19].e[2]    = &(m_poly.e[29]);
            m_poly.t[19].v[0]    = &(m_poly.v[0]);
            m_poly.t[19].v[1]    = &(m_poly.v[2]);
            m_poly.t[19].v[2]    = &(m_poly.v[9]);
            m_poly.t[19].a[0]    = &(m_poly.t[0]);
            m_poly.t[19].a[1]    = &(m_poly.t[5]);
            m_poly.t[19].a[2]    = &(m_poly.t[18]);
        }
        break;
    case 4:
        {
            // ============================================================
            // Create a tetrahedron whose vertices are located at:
            //   (1,0,0), (0,1,0), (0,0,1), (sqrt(1/3),sqrt(1/3),sqrt(1/3))
            //

            m_poly.e = new edge_type[6];        // 6 edges
            m_poly.m = 6;
            m_poly.v = new vertex_type[4];      // 4 vertices
            m_poly.n = 4;
            m_poly.t = new triangle_type[4];    // 4 triangles
            m_poly.r = 4;
            
            // ------------------------------------------------------------
            // vertices
            m_poly.v[0].p        = new point_type;
            m_poly.v[0].p->v[0]  = value_type(1.0);
            m_poly.v[0].p->v[1]  = value_type(0.0);
            m_poly.v[0].p->v[2]  = value_type(0.0);
            m_poly.v[0].e        = new edge_type*[3];
            m_poly.v[0].e[0]     = &(m_poly.e[0]);
            m_poly.v[0].e[1]     = &(m_poly.e[2]);
            m_poly.v[0].e[2]     = &(m_poly.e[3]);
            m_poly.v[0].m        = 3;

            m_poly.v[1].p        = new point_type;
            m_poly.v[1].p->v[0]  = value_type(0.0);
            m_poly.v[1].p->v[1]  = value_type(1.0);
            m_poly.v[1].p->v[2]  = value_type(0.0);
            m_poly.v[1].e        = new edge_type*[3];
            m_poly.v[1].e[0]     = &(m_poly.e[0]);
            m_poly.v[1].e[1]     = &(m_poly.e[1]);
            m_poly.v[1].e[2]     = &(m_poly.e[4]);
            m_poly.v[1].m        = 3;

            m_poly.v[2].p        = new point_type;
            m_poly.v[2].p->v[0]  = value_type(0.0);
            m_poly.v[2].p->v[1]  = value_type(0.0);
            m_poly.v[2].p->v[2]  = value_type(1.0);
            m_poly.v[2].e        = new edge_type*[3];
            m_poly.v[2].e[0]     = &(m_poly.e[1]);
            m_poly.v[2].e[1]     = &(m_poly.e[2]);
            m_poly.v[2].e[2]     = &(m_poly.e[5]);
            m_poly.v[2].m        = 3;

            m_poly.v[3].p        = new point_type;
            m_poly.v[3].p->v[0]  = value_type(sqrt(1.0/3.0));
            m_poly.v[3].p->v[1]  = value_type(sqrt(1.0/3.0));
            m_poly.v[3].p->v[2]  = value_type(sqrt(1.0/3.0));
            m_poly.v[3].e        = new edge_type*[3];
            m_poly.v[3].e[0]     = &(m_poly.e[3]);
            m_poly.v[3].e[1]     = &(m_poly.e[4]);
            m_poly.v[3].e[2]     = &(m_poly.e[5]);
            m_poly.v[3].m        = 3;
            
            // ------------------------------------------------------------
            // edges
            m_poly.e[0].v[0]     = &(m_poly.v[0]);
            m_poly.e[0].v[1]     = &(m_poly.v[1]);
            m_poly.e[0].t[0]     = &(m_poly.t[0]);
            m_poly.e[0].t[1]     = &(m_poly.t[1]);

            m_poly.e[1].v[0]     = &(m_poly.v[1]);
            m_poly.e[1].v[1]     = &(m_poly.v[2]);
            m_poly.e[1].t[0]     = &(m_poly.t[0]);
            m_poly.e[1].t[1]     = &(m_poly.t[2]);

            m_poly.e[2].v[0]     = &(m_poly.v[2]);
            m_poly.e[2].v[1]     = &(m_poly.v[0]);
            m_poly.e[2].t[0]     = &(m_poly.t[0]);
            m_poly.e[2].t[1]     = &(m_poly.t[3]);

            m_poly.e[3].v[0]     = &(m_poly.v[0]);
            m_poly.e[3].v[1]     = &(m_poly.v[3]);
            m_poly.e[3].t[0]     = &(m_poly.t[1]);
            m_poly.e[3].t[1]     = &(m_poly.t[3]);

            m_poly.e[4].v[0]     = &(m_poly.v[1]);
            m_poly.e[4].v[1]     = &(m_poly.v[3]);
            m_poly.e[4].t[0]     = &(m_poly.t[1]);
            m_poly.e[4].t[1]     = &(m_poly.t[2]);

            m_poly.e[5].v[0]     = &(m_poly.v[2]);
            m_poly.e[5].v[1]     = &(m_poly.v[3]);
            m_poly.e[5].t[0]     = &(m_poly.t[2]);
            m_poly.e[5].t[1]     = &(m_poly.t[3]);
            
            // ------------------------------------------------------------
            // triangles
            m_poly.t[0].e[0]     = &(m_poly.e[2]);
            m_poly.t[0].e[1]     = &(m_poly.e[1]);
            m_poly.t[0].e[2]     = &(m_poly.e[0]);
            m_poly.t[0].v[0]     = &(m_poly.v[0]);
            m_poly.t[0].v[1]     = &(m_poly.v[2]);
            m_poly.t[0].v[2]     = &(m_poly.v[1]);
            m_poly.t[0].a[0]     = &(m_poly.t[3]);
            m_poly.t[0].a[1]     = &(m_poly.t[2]);
            m_poly.t[0].a[2]     = &(m_poly.t[1]);

            m_poly.t[1].e[0]     = &(m_poly.e[0]);
            m_poly.t[1].e[1]     = &(m_poly.e[4]);
            m_poly.t[1].e[2]     = &(m_poly.e[3]);
            m_poly.t[1].v[0]     = &(m_poly.v[0]);
            m_poly.t[1].v[1]     = &(m_poly.v[1]);
            m_poly.t[1].v[2]     = &(m_poly.v[3]);
            m_poly.t[1].a[0]     = &(m_poly.t[0]);
            m_poly.t[1].a[1]     = &(m_poly.t[2]);
            m_poly.t[1].a[2]     = &(m_poly.t[3]);

            m_poly.t[2].e[0]     = &(m_poly.e[1]);
            m_poly.t[2].e[1]     = &(m_poly.e[5]);
            m_poly.t[2].e[2]     = &(m_poly.e[4]);
            m_poly.t[2].v[0]     = &(m_poly.v[1]);
            m_poly.t[2].v[1]     = &(m_poly.v[2]);
            m_poly.t[2].v[2]     = &(m_poly.v[3]);
            m_poly.t[2].a[0]     = &(m_poly.t[0]);
            m_poly.t[2].a[1]     = &(m_poly.t[3]);
            m_poly.t[2].a[2]     = &(m_poly.t[1]);

            m_poly.t[3].e[0]     = &(m_poly.e[3]);
            m_poly.t[3].e[1]     = &(m_poly.e[5]);
            m_poly.t[3].e[2]     = &(m_poly.e[2]);
            m_poly.t[3].v[0]     = &(m_poly.v[0]);
            m_poly.t[3].v[1]     = &(m_poly.v[3]);
            m_poly.t[3].v[2]     = &(m_poly.v[2]);
            m_poly.t[3].a[0]     = &(m_poly.t[1]);
            m_poly.t[3].a[1]     = &(m_poly.t[2]);
            m_poly.t[3].a[2]     = &(m_poly.t[0]);

        }
        break;
    case 8:
    default:
        {
            // ============================================================
            // Create an octahedron whose vertices are located at:
            //   (0,0,1) (1,0,0) (0,1,0) (-1,0,0) (0,-1,0) (0,0,-1)
            //
            
            m_poly.e = new edge_type[12];       // 12 edges
            m_poly.m = 12;
            m_poly.v = new vertex_type[6];      // 6  vertices
            m_poly.n = 6;
            m_poly.t = new triangle_type[8];    // 8  triangles
            m_poly.r = 8;

            // ------------------------------------------------------------
            // vertices
            m_poly.v[0].p        = new point_type;
            m_poly.v[0].p->v[0]  = value_type(0.0);
            m_poly.v[0].p->v[1]  = value_type(0.0);
            m_poly.v[0].p->v[2]  = value_type(1.0);
            m_poly.v[0].e        = new edge_type*[4];
            m_poly.v[0].e[0]     = &(m_poly.e[0]);
            m_poly.v[0].e[1]     = &(m_poly.e[1]);
            m_poly.v[0].e[2]     = &(m_poly.e[2]);
            m_poly.v[0].e[3]     = &(m_poly.e[3]);
            m_poly.v[0].m        = 4;

            m_poly.v[1].p        = new point_type;
            m_poly.v[1].p->v[0]  = value_type(1.0);
            m_poly.v[1].p->v[1]  = value_type(0.0);
            m_poly.v[1].p->v[2]  = value_type(0.0);
            m_poly.v[1].e        = new edge_type*[4];
            m_poly.v[1].e[0]     = &(m_poly.e[0]);
            m_poly.v[1].e[1]     = &(m_poly.e[4]);
            m_poly.v[1].e[2]     = &(m_poly.e[7]);
            m_poly.v[1].e[3]     = &(m_poly.e[8]);
            m_poly.v[1].m        = 4;

            m_poly.v[2].p        = new point_type;
            m_poly.v[2].p->v[0]  = value_type(0.0);
            m_poly.v[2].p->v[1]  = value_type(1.0);
            m_poly.v[2].p->v[2]  = value_type(0.0);
            m_poly.v[2].e        = new edge_type*[4];
            m_poly.v[2].e[0]     = &(m_poly.e[1]);
            m_poly.v[2].e[1]     = &(m_poly.e[4]);
            m_poly.v[2].e[2]     = &(m_poly.e[5]);
            m_poly.v[2].e[3]     = &(m_poly.e[9]);
            m_poly.v[2].m        = 4;

            m_poly.v[3].p        = new point_type;
            m_poly.v[3].p->v[0]  = value_type(-1.0);
            m_poly.v[3].p->v[1]  = value_type(0.0);
            m_poly.v[3].p->v[2]  = value_type(0.0);
            m_poly.v[3].e        = new edge_type*[4];
            m_poly.v[3].e[0]     = &(m_poly.e[2]);
            m_poly.v[3].e[1]     = &(m_poly.e[5]);
            m_poly.v[3].e[2]     = &(m_poly.e[6]);
            m_poly.v[3].e[3]     = &(m_poly.e[10]);
            m_poly.v[3].m        = 4;
            
            m_poly.v[4].p        = new point_type;
            m_poly.v[4].p->v[0]  = value_type(0.0);
            m_poly.v[4].p->v[1]  = value_type(-1.0);
            m_poly.v[4].p->v[2]  = value_type(0.0);
            m_poly.v[4].e        = new edge_type*[4];
            m_poly.v[4].e[0]     = &(m_poly.e[3]);
            m_poly.v[4].e[1]     = &(m_poly.e[6]);
            m_poly.v[4].e[2]     = &(m_poly.e[7]);
            m_poly.v[4].e[3]     = &(m_poly.e[11]);
            m_poly.v[4].m        = 4;

            m_poly.v[5].p        = new point_type;
            m_poly.v[5].p->v[0]  = value_type(0.0);
            m_poly.v[5].p->v[1]  = value_type(0.0);
            m_poly.v[5].p->v[2]  = value_type(-1.0);
            m_poly.v[5].e        = new edge_type*[4];
            m_poly.v[5].e[0]     = &(m_poly.e[8]);
            m_poly.v[5].e[1]     = &(m_poly.e[9]);
            m_poly.v[5].e[2]     = &(m_poly.e[10]);
            m_poly.v[5].e[3]     = &(m_poly.e[11]);
            m_poly.v[5].m        = 4;

            // ------------------------------------------------------------
            // edges
            m_poly.e[0].v[0]     = &(m_poly.v[0]);
            m_poly.e[0].v[1]     = &(m_poly.v[1]);
            m_poly.e[0].t[0]     = &(m_poly.t[3]);
            m_poly.e[0].t[1]     = &(m_poly.t[0]);

            m_poly.e[1].v[0]     = &(m_poly.v[0]);
            m_poly.e[1].v[1]     = &(m_poly.v[2]);
            m_poly.e[1].t[0]     = &(m_poly.t[0]);
            m_poly.e[1].t[1]     = &(m_poly.t[1]);

            m_poly.e[2].v[0]     = &(m_poly.v[0]);
            m_poly.e[2].v[1]     = &(m_poly.v[3]);
            m_poly.e[2].t[0]     = &(m_poly.t[1]);
            m_poly.e[2].t[1]     = &(m_poly.t[2]);

            m_poly.e[3].v[0]     = &(m_poly.v[0]);
            m_poly.e[3].v[1]     = &(m_poly.v[4]);
            m_poly.e[3].t[0]     = &(m_poly.t[2]);
            m_poly.e[3].t[1]     = &(m_poly.t[3]);

            m_poly.e[4].v[0]     = &(m_poly.v[1]);
            m_poly.e[4].v[1]     = &(m_poly.v[2]);
            m_poly.e[4].t[0]     = &(m_poly.t[0]);
            m_poly.e[4].t[1]     = &(m_poly.t[4]);

            m_poly.e[5].v[0]     = &(m_poly.v[2]);
            m_poly.e[5].v[1]     = &(m_poly.v[3]);
            m_poly.e[5].t[0]     = &(m_poly.t[1]);
            m_poly.e[5].t[1]     = &(m_poly.t[5]);

            m_poly.e[6].v[0]     = &(m_poly.v[3]);
            m_poly.e[6].v[1]     = &(m_poly.v[4]);
            m_poly.e[6].t[0]     = &(m_poly.t[2]);
            m_poly.e[6].t[1]     = &(m_poly.t[6]);

            m_poly.e[7].v[0]     = &(m_poly.v[4]);
            m_poly.e[7].v[1]     = &(m_poly.v[1]);
            m_poly.e[7].t[0]     = &(m_poly.t[3]);
            m_poly.e[7].t[1]     = &(m_poly.t[7]);

            m_poly.e[8].v[0]     = &(m_poly.v[1]);
            m_poly.e[8].v[1]     = &(m_poly.v[5]);
            m_poly.e[8].t[0]     = &(m_poly.t[7]);
            m_poly.e[8].t[1]     = &(m_poly.t[4]);

            m_poly.e[9].v[0]     = &(m_poly.v[2]);
            m_poly.e[9].v[1]     = &(m_poly.v[5]);
            m_poly.e[9].t[0]     = &(m_poly.t[4]);
            m_poly.e[9].t[1]     = &(m_poly.t[5]);

            m_poly.e[10].v[0]    = &(m_poly.v[3]);
            m_poly.e[10].v[1]    = &(m_poly.v[5]);
            m_poly.e[10].t[0]    = &(m_poly.t[5]);
            m_poly.e[10].t[1]    = &(m_poly.t[6]);

            m_poly.e[11].v[0]    = &(m_poly.v[4]);
            m_poly.e[11].v[1]    = &(m_poly.v[5]);
            m_poly.e[11].t[0]    = &(m_poly.t[6]);
            m_poly.e[11].t[1]    = &(m_poly.t[7]);

            // ------------------------------------------------------------
            // triangles
            m_poly.t[0].v[0]     = &(m_poly.v[0]);
            m_poly.t[0].v[1]     = &(m_poly.v[1]);
            m_poly.t[0].v[2]     = &(m_poly.v[2]);
            m_poly.t[0].e[0]     = &(m_poly.e[0]);
            m_poly.t[0].e[1]     = &(m_poly.e[4]);
            m_poly.t[0].e[2]     = &(m_poly.e[1]);
            m_poly.t[0].a[0]     = &(m_poly.t[3]);
            m_poly.t[0].a[1]     = &(m_poly.t[4]);
            m_poly.t[0].a[2]     = &(m_poly.t[1]);

            m_poly.t[1].v[0]     = &(m_poly.v[0]);
            m_poly.t[1].v[1]     = &(m_poly.v[2]);
            m_poly.t[1].v[2]     = &(m_poly.v[3]);
            m_poly.t[1].e[0]     = &(m_poly.e[1]);
            m_poly.t[1].e[1]     = &(m_poly.e[5]);
            m_poly.t[1].e[2]     = &(m_poly.e[2]);
            m_poly.t[1].a[0]     = &(m_poly.t[0]);
            m_poly.t[1].a[1]     = &(m_poly.t[5]);
            m_poly.t[1].a[2]     = &(m_poly.t[2]);

            m_poly.t[2].v[0]     = &(m_poly.v[0]);
            m_poly.t[2].v[1]     = &(m_poly.v[3]);
            m_poly.t[2].v[2]     = &(m_poly.v[4]);
            m_poly.t[2].e[0]     = &(m_poly.e[2]);
            m_poly.t[2].e[1]     = &(m_poly.e[6]);
            m_poly.t[2].e[2]     = &(m_poly.e[3]);
            m_poly.t[2].a[0]     = &(m_poly.t[1]);
            m_poly.t[2].a[1]     = &(m_poly.t[6]);
            m_poly.t[2].a[2]     = &(m_poly.t[3]);

            m_poly.t[3].v[0]     = &(m_poly.v[0]);
            m_poly.t[3].v[1]     = &(m_poly.v[4]);
            m_poly.t[3].v[2]     = &(m_poly.v[1]);
            m_poly.t[3].e[0]     = &(m_poly.e[3]);
            m_poly.t[3].e[1]     = &(m_poly.e[7]);
            m_poly.t[3].e[2]     = &(m_poly.e[0]);
            m_poly.t[3].a[0]     = &(m_poly.t[2]);
            m_poly.t[3].a[1]     = &(m_poly.t[7]);
            m_poly.t[3].a[2]     = &(m_poly.t[0]);

            m_poly.t[4].v[0]     = &(m_poly.v[5]);
            m_poly.t[4].v[1]     = &(m_poly.v[2]);
            m_poly.t[4].v[2]     = &(m_poly.v[1]);
            m_poly.t[4].e[0]     = &(m_poly.e[9]);
            m_poly.t[4].e[1]     = &(m_poly.e[4]);
            m_poly.t[4].e[2]     = &(m_poly.e[8]);
            m_poly.t[4].a[0]     = &(m_poly.t[5]);
            m_poly.t[4].a[1]     = &(m_poly.t[0]);
            m_poly.t[4].a[2]     = &(m_poly.t[7]);

            m_poly.t[5].v[0]     = &(m_poly.v[5]);
            m_poly.t[5].v[1]     = &(m_poly.v[3]);
            m_poly.t[5].v[2]     = &(m_poly.v[2]);
            m_poly.t[5].e[0]     = &(m_poly.e[10]);
            m_poly.t[5].e[1]     = &(m_poly.e[5]);
            m_poly.t[5].e[2]     = &(m_poly.e[9]);
            m_poly.t[5].a[0]     = &(m_poly.t[6]);
            m_poly.t[5].a[1]     = &(m_poly.t[1]);
            m_poly.t[5].a[2]     = &(m_poly.t[4]);

            m_poly.t[6].v[0]     = &(m_poly.v[5]);
            m_poly.t[6].v[1]     = &(m_poly.v[4]);
            m_poly.t[6].v[2]     = &(m_poly.v[3]);
            m_poly.t[6].e[0]     = &(m_poly.e[11]);
            m_poly.t[6].e[1]     = &(m_poly.e[6]);
            m_poly.t[6].e[2]     = &(m_poly.e[10]);
            m_poly.t[6].a[0]     = &(m_poly.t[7]);
            m_poly.t[6].a[1]     = &(m_poly.t[2]);
            m_poly.t[6].a[2]     = &(m_poly.t[5]);

            m_poly.t[7].v[0]     = &(m_poly.v[5]);
            m_poly.t[7].v[1]     = &(m_poly.v[1]);
            m_poly.t[7].v[2]     = &(m_poly.v[4]);
            m_poly.t[7].e[0]     = &(m_poly.e[8]);
            m_poly.t[7].e[1]     = &(m_poly.e[7]);
            m_poly.t[7].e[2]     = &(m_poly.e[11]);
            m_poly.t[7].a[0]     = &(m_poly.t[4]);
            m_poly.t[7].a[1]     = &(m_poly.t[3]);
            m_poly.t[7].a[2]     = &(m_poly.t[6]);
        }
        break;
    }
}

///////////////////////////////////////////////////////////////////////////
template <typename _Real>
void sphere_tessellation<_Real>::tessellate(int level)
{
    size_type   n, // number of vertices
                m, // number of edges
                r; // number of triangles
    int l;

    if (level <= 0) return;

    n = m_poly.n;
    m = m_poly.m;
    r = m_poly.r;

    expand(level);

    // subdivide polyhedron
    for (l = 1; l <= level; ++l) {
        
        size_t          i0, i1, i2, i3, im, in, ir, ir1, ir2, ir3, orient;
        triangle_type   *t0, *t1, *t2, *t3, *ta0, *ta1, *ta2;
        vertex_type     *v0, *v1, *v2, *vm0, *vm1, *vm2;
        edge_type       *e0, *e1, *e2, *en0, *en1, *en2;

        for (im = 0; im < m; ++im) {
            e0  = &(m_poly.e[im]);
            e1  = &(m_poly.e[im + m]);

            point_type& p0  = *(e0->v[0]->p);
            point_type& p1  = *(e0->v[1]->p);
            point_type  pm  = (p0 + p1) * value_type(0.5);
            
            vm0 = &(m_poly.v[im + n]);
            vm0->p          = new point_type;
            intersect_sphere(vm0->p, m_poly.c, pm, value_type(1));

            vm0->m          = 0;
            vm0->e          = new edge_type*[6];
            
            // e0 = <v0,v1>, vm0 = (v0+v1)/2
            //    ==> e0 = <v0,vm0>, e1 = <v1,vm0>
            e1->v[0]        = e0->v[1];
            e0->v[1]        = vm0;
            e1->v[1]        = vm0;
            
            e0->t[0]        = NULL;
            e1->t[0]        = NULL;
        } // for im

        update_centroid(m + n);
        
        ir1 = r;
        ir2 = r * 2;
        ir3 = r * 3;
        im  = 2 * m;

        for (ir = 0; ir < r; ++ir, ++ir1, ++ir2, ++ir3) {

            t0  = &(m_poly.t[ir]);
            
            e0          = t0->e[0]; //
            e1          = t0->e[1]; // edges
            e2          = t0->e[2]; //

            v0          = t0->v[0]; //
            v1          = t0->v[1]; // vertices
            v2          = t0->v[2]; //

            ta0         = t0->a[0]; //
            ta1         = t0->a[1]; // adjacent triangles
            ta2         = t0->a[2]; //

            // midpoints of triangle edges
            vm0         = &(m_poly.v[n + edge_offset(m_poly, e0)]);
            vm1         = &(m_poly.v[n + edge_offset(m_poly, e1)]);
            vm2         = &(m_poly.v[n + edge_offset(m_poly, e2)]);

            // new edges
            en0         = &(m_poly.e[im++]);
            en0->v[0]   = vm2;
            en0->v[1]   = vm0;
            en1         = &(m_poly.e[im++]);
            en1->v[0]   = vm0;
            en1->v[1]   = vm1;
            en2         = &(m_poly.e[im++]);
            en2->v[0]   = vm1;
            en2->v[1]   = vm2;
            
            // construct triangle t1
            t1          = &(m_poly.t[ir1]);
            t1->v[0]    = v0;
            t1->v[1]    = vm0;
            t1->v[2]    = vm2;
            t1->e[0]    = (e0->v[0] == v0 ? e0 : e0 + m);
            t1->e[1]    = en0;
            t1->e[2]    = (e2->v[0] == v0 ? e2 : e2 + m);
            t1->a[1]    = t0;

            // construct triangle t2
            t2          = &(m_poly.t[ir2]);
            t2->v[0]    = vm0;
            t2->v[1]    = v1;
            t2->v[2]    = vm1;
            t2->e[0]    = (e0->v[0] == v1 ? e0 : e0 + m);
            t2->e[1]    = (e1->v[0] == v1 ? e1 : e1 + m);
            t2->e[2]    = en1;
            t2->a[2]    = t0;

            // construct triangle t3
            t3          = &(m_poly.t[ir3]);
            t3->v[0]    = vm2;
            t3->v[1]    = vm1;
            t3->v[2]    = v2;
            t3->e[0]    = en2;
            t3->e[1]    = (e1->v[0] == v2 ? e1 : e1 + m);
            t3->e[2]    = (e2->v[0] == v2 ? e2 : e2 + m);
            t3->a[0]    = t0;

            // edge triangle pointers
            en0->t[0]   = t0;
            en0->t[1]   = t1;
            en1->t[0]   = t0;
            en1->t[1]   = t2;
            en2->t[0]   = t0;
            en2->t[1]   = t3;

            // get the indices for the subdivided triangles of ta0
            i0 = triangle_offset(m_poly, ta0);
            i1 = i0 + r;
            i2 = i1 + r;
            i3 = i2 + r;

            // determine the orientation of ta0 relative to t0 and set the
            // adjacency links.
            orient = adjacent_orient(t0, ta0);
            if (orient == 0) {
                t1->a[0]            = &(m_poly.t[i2]);
                t2->a[0]            = &(m_poly.t[i1]);
                m_poly.t[i2].a[0]   = t1;
                m_poly.t[i1].a[0]   = t2;
            }
            else if (orient == 1) {
                t1->a[0]            = &(m_poly.t[i3]);
                t2->a[0]            = &(m_poly.t[i2]);
                m_poly.t[i3].a[1]   = t1;
                m_poly.t[i2].a[1]   = t2;
            }
            else {
                t1->a[0]            = &(m_poly.t[i1]);
                t2->a[0]            = &(m_poly.t[i3]);
                m_poly.t[i1].a[2]   = t1;
                m_poly.t[i3].a[2]   = t2;
            }

            // get the indices for the subdivided triangles of ta1
            i0 = triangle_offset(m_poly, ta1);
            i1 = i0 + r;
            i2 = i1 + r;
            i3 = i2 + r;

            // determine the orientation of ta1 relative to t1 and set the
            // adjacency links.
            orient = adjacent_orient(t0, ta1);
            if (orient == 0) {
                t2->a[1]            = &(m_poly.t[i2]);
                t3->a[1]            = &(m_poly.t[i1]);
                m_poly.t[i2].a[0]   = t2;
                m_poly.t[i1].a[0]   = t3;
            }
            else if (orient == 1) {
                t2->a[1]            = &(m_poly.t[i3]);
                t3->a[1]            = &(m_poly.t[i2]);
                m_poly.t[i3].a[1]   = t2;
                m_poly.t[i2].a[1]   = t3;
            }
            else {
                t2->a[1]            = &(m_poly.t[i1]);
                t3->a[1]            = &(m_poly.t[i3]);
                m_poly.t[i1].a[2]   = t2;
                m_poly.t[i3].a[2]   = t3;
            }

            // get the indices for the subdivided triangles of ta2
            i0 = triangle_offset(m_poly, ta2);
            i1 = i0 + r;
            i2 = i1 + r;
            i3 = i2 + r;

            // determine the orientation of ta2 relative to t2 and set the
            // adjacency links.
            orient = adjacent_orient(t0, ta2);
            if (orient == 0) {
                t3->a[2]            = &(m_poly.t[i2]);
                t1->a[2]            = &(m_poly.t[i1]);
                m_poly.t[i2].a[0]   = t3;
                m_poly.t[i1].a[0]   = t1;
            }
            else if (orient == 1) {
                t3->a[2]            = &(m_poly.t[i3]);
                t1->a[2]            = &(m_poly.t[i2]);
                m_poly.t[i3].a[1]   = t3;
                m_poly.t[i2].a[1]   = t1;
            }
            else {
                t3->a[2]            = &(m_poly.t[i1]);
                t1->a[2]            = &(m_poly.t[i3]);
                m_poly.t[i1].a[2]   = t3;
                m_poly.t[i3].a[2]   = t1;
            }

            // add edge links to midpoint vertices
            if (vm0->m == 0) {
                // add 4 edges
                vm0->m      = 4;
                vm0->e[0]   = t1->e[0];
                vm0->e[1]   = t1->e[1];
                vm0->e[2]   = t2->e[2];
                vm0->e[3]   = t2->e[0];
            }
            else if (vm0->m == 4) {
                // add 2 edges
                vm0->m      = 6;
                vm0->e[4]   = t1->e[1];
                vm0->e[5]   = t2->e[2];
            }

            if (vm1->m == 0) {
                // add 4 edges
                vm1->m      = 4;
                vm1->e[0]   = t2->e[1];
                vm1->e[1]   = t2->e[2];
                vm1->e[2]   = t3->e[0];
                vm1->e[3]   = t3->e[1];
            }
            else if (vm1->m == 4) {
                // add 2 edges
                vm1->m      = 6;
                vm1->e[4]   = t2->e[2];
                vm1->e[5]   = t3->e[0];
            }

            if (vm2->m == 0) {
                // add 4 edges
                vm2->m      = 4;
                vm2->e[0]   = t1->e[2];
                vm2->e[1]   = t1->e[1];
                vm2->e[2]   = t3->e[0];
                vm2->e[3]   = t3->e[2];
            }
            else if (vm2->m == 4) {
                // add 2 edges
                vm2->m      = 6;
                vm2->e[4]   = t1->e[1];
                vm2->e[5]   = t3->e[0];
            }

        } // for ir

        ir1 = r;
        ir2 = r * 2;
        ir3 = r * 3;
        for (ir = 0; ir < r; ++ir, ++ir1, ++ir2, ++ir3) {
            t0  = &(m_poly.t[ir ]); //
            t1  = &(m_poly.t[ir1]); // subdivided triangles
            t2  = &(m_poly.t[ir2]); //
            t3  = &(m_poly.t[ir3]); //
            
            vm0 = &m_poly.v[n + edge_offset(m_poly, t0->e[0])];
            vm1 = &m_poly.v[n + edge_offset(m_poly, t0->e[1])];
            vm2 = &m_poly.v[n + edge_offset(m_poly, t0->e[2])];

            t0->e[0] = t1->e[1];    //
            t0->e[1] = t2->e[2];    // edges for middle triangle
            t0->e[2] = t3->e[0];    //

            t0->v[0] = vm2;         //
            t0->v[1] = vm0;         // vertices for middle triangle
            t0->v[2] = vm1;         //

            t0->a[0] = t1;          //
            t0->a[1] = t2;          // adjacencies for middle triangle
            t0->a[2] = t3;          //

            // edge triangle pointers
            if (t1->e[0]->t[0] == 0) {
                t1->e[0]->t[0] = t1;
                t1->e[0]->t[1] = t1->a[0];
            }
            if (t1->e[2]->t[0] == 0) {
                t1->e[2]->t[0] = t1;
                t1->e[2]->t[1] = t1->a[2];
            }

            if (t2->e[0]->t[0] == 0) {
                t2->e[0]->t[0] = t2;
                t2->e[0]->t[1] = t2->a[0];
            }
            if (t2->e[1]->t[0] == 0) {
                t2->e[1]->t[0] = t2;
                t2->e[1]->t[1] = t2->a[1];
            }

            if (t3->e[1]->t[0] == 0) {
                t3->e[1]->t[0] = t3;
                t3->e[1]->t[1] = t3->a[1];
            }
            if (t3->e[2]->t[0] == 0) {
                t3->e[2]->t[0] = t3;
                t3->e[2]->t[1] = t3->a[2];
            }

        } // for ir

        // adjust edge pointers of original vertices to account for the
        // edge splitting.
        for (in = 0; in < n; ++in) {
            v0 = &m_poly.v[in];
            for (im = 0; im < v0->m; ++im) {
                e0 = v0->e[im];
                if (e0->v[0] != v0 && e0->v[1] != v0)
                    v0->e[im] = &m_poly.e[m + edge_offset(m_poly, e0)];
            }
        }

        // update numbers
        n = n + m;
        m = 2 * m + 3 * r;
        r = 4 * r;

    } // for l
}

///////////////////////////////////////////////////////////////////////////
template <typename _Real>
void sphere_tessellation<_Real>::clear()
{
    for (size_type i = 0; i < m_poly.n; ++i) {
        delete m_poly.v[i].p;
        delete [] m_poly.v[i].e;
    }

    m_poly.c.clear();
    SAFE_DELETE_ARRAY(m_poly.e);
    SAFE_DELETE_ARRAY(m_poly.v);
    SAFE_DELETE_ARRAY(m_poly.t);

    m_poly.n = 0;
    m_poly.m = 0;
    m_poly.r = 0;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Real>
int sphere_tessellation<_Real>::adjacent_orient(triangle_const_pointer t,
                                                triangle_const_pointer ta
                                                )
{
    if (ta->a[0] == t) return 0;
    else if (ta->a[1] == t) return 1;
    return 2;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Real>
void sphere_tessellation<_Real>::update_centroid(size_type n)
{
    size_type   i;
    value_type  ninv;

    m_poly.c.clear(); // set centroid to (0,0,0).
    for (i = 0; i < n; ++i) m_poly.c += *(m_poly.v[i].p);
    ninv = value_type(1.0 / n);
    m_poly.c *= ninv;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Real>
void sphere_tessellation<_Real>::expand(int level)
{
    size_type l, in, im, ir;

    update_centroid(m_poly.n);

    // make a shallow copy of the current polyhedron.
    polyhedron_type seed = m_poly;

    // compute number of vertices/edges/triangles in the final polyhedron.
    for (l = 1; l <= (size_type)level; ++l) {
        m_poly.n = m_poly.n + m_poly.m;
        m_poly.m = 2 * m_poly.m + 3 * m_poly.r;
        m_poly.r = 4 * m_poly.r;
    } // for l

    m_poly.e = new edge_type[m_poly.m];
    m_poly.v = new vertex_type[m_poly.n];
    m_poly.t = new triangle_type[m_poly.r];

    for (in = 0; in < seed.n; ++in) {

        // transfer points from seed to m_poly
        vertex_type& vs = seed.v[in];
        vertex_type& vp = m_poly.v[in];
        vp.p = vs.p;

        // duplicate edge information
        vp.m = vs.m;
        l = vs.m * sizeof(edge_type*);
        vp.e = new edge_type*[l];

        for (im = 0; im < vs.m; ++im)
            vp.e[im] = &m_poly.e[edge_offset(seed, vs.e[im])];

        delete [] vs.e;
    }

    delete [] seed.v;

    // duplicate edges
    for (im = 0; im < seed.m; ++im) {
        edge_type& es = seed.e[im];
        edge_type& ep = m_poly.e[im];

        for (l = 0; l < 2; ++l) {
            ep.v[l] = &m_poly.v[vertex_offset(seed, es.v[l])];
            ep.t[l] = &m_poly.t[triangle_offset(seed, es.t[im])];
        }
    }

    delete [] seed.e;

    // duplicate triangles
    for (ir = 0; ir < seed.r; ++ir) {
        triangle_type& ts = seed.t[ir];
        triangle_type& tp = m_poly.t[ir];

        for (l = 0; l < 3; ++l) {
            tp.e[l] = &m_poly.e[edge_offset(seed, ts.e[l])];
            tp.v[l] = &m_poly.v[vertex_offset(seed, ts.v[l])];
            tp.a[l] = &m_poly.t[triangle_offset(seed, ts.a[l])];
        }
    }

    delete [] seed.t;
}

        } // namespace
    } // namespace
} // namespace

#endif // ___SPHERE_TESSELLATION_CXX___
