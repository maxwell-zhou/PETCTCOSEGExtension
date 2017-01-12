/*
 ==========================================================================
 |   
 |   $Id: intersect.hxx 58 2005-01-21 01:44:10Z kangli $
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
 
#ifndef ___INTERSECT_HXX___
#   define ___INTERSECT_HXX___

#   include <optnet/_base/point3.hxx>
#   include <cassert>
#   include <cmath>

/// @namespace optnet
namespace optnet {
    /// @namespace xtra
    namespace xtra {
        /// @namespace graphics
        namespace graphics {

///////////////////////////////////////////////////////////////////////////
///  Compute the intersecting point of a sphere and a ray emanated
///  from the center of the sphere.
///
///  @param pt      The ray-sphere intersecting point.
///  @param o       The center of the sphere.
///  @param ray     The point through which the ray passes.
///  @param radius  The radius of the sphere.
///
///////////////////////////////////////////////////////////////////////////
template <class _T>
inline void intersect_sphere(point3<_T>*        pt,
                             const point3<_T>&  o,
                             const point3<_T>&  ray,
                             _T                 radius = 1
                             )
{
    assert(pt != 0);

    point3<_T> dir = ray - o;

    _T a = dir.length2();
    _T b = 2 * dir.dot(o);
    _T c = o.length2() - radius*radius;
    
    double s = b * b - 4 * a * c;
    if (s < 0) s = -s;
    s = sqrt(s);
    
    double t = (-b + s) / (2.0 * a);

    //Algorithm:
    //  pt(t) = o + t * (ray - o)
    //  where t satisfies: radius^2 = pt(t) * pt(t)
    *pt = o + dir * t;
}

        } // namespace
    } // namespace
} // namespace

#endif // ___SPHERE_TESSELLATION_CXX___
