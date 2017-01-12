/*
 ==========================================================================
 |   
 |   $Id: levelset3_initializer.hxx 150 2005-02-08 18:25:47Z kangli $
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

#ifndef ___LEVELSET3_INITIALIZER_HXX___
#   define ___LEVELSET3_INITIALIZER_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_base/point3.hxx>
#   include <algorithm>
#   include <vector>
#   include <cmath>

namespace optnet {

///////////////////////////////////////////////////////////////////////
///  Initialize a levelset function representing a sphere.
///
///  @param levelset  The level set function to be initialized.
///  @param center    The center of the sphere.
///  @param radius    The radius of the sphere.
///
///////////////////////////////////////////////////////////////////////
template <typename _Ty>
static void levelset3_init_sphere(array_base<_Ty>&      levelset,
                                  const point3<double>& center,
                                  double                radius
                                  )
{
    size_t i0, i1, i2;

    for (i2 = 0; i2 < levelset.size_2(); ++i2) {
        for (i1 = 0; i1 < levelset.size_1(); ++i1) {
            for (i0 = 0; i0 < levelset.size_0(); ++i0) {

                double d = sqrt((i0 - center.v[0]) * (i0 - center.v[0]) +
                                (i1 - center.v[1]) * (i1 - center.v[1]) +
                                (i2 - center.v[2]) * (i2 - center.v[2]));

                levelset(i0, i1, i2) = static_cast<_Ty>(radius - d);

            } // i0
        } // i1
    } // i2
}

///////////////////////////////////////////////////////////////////////
///  Initialize a levelset function representing two spheres.
///
///  @param levelset  The level set function to be initialized.
///  @param center1   The center of the first  sphere.
///  @param radius1   The radius of the first  sphere.
///  @param center2   The center of the second sphere.
///  @param radius2   The radius of the second sphere.
///
///////////////////////////////////////////////////////////////////////
template <typename _Ty>
static void levelset3_init_two_spheres(array_base<_Ty>&      levelset,
                                       const point3<double>& center1,
                                       double                radius1,
                                       const point3<double>& center2,
                                       double                radius2
                                       )
{
    size_t i0, i1, i2;

    for (i2 = 0; i2 < levelset.size_2(); ++i2) {
        for (i1 = 0; i1 < levelset.size_1(); ++i1) {
            for (i0 = 0; i0 < levelset.size_0(); ++i0) {

                double d1 = radius1 - sqrt((i0 - center1.v[0]) * (i0 - center1.v[0]) +
                                           (i1 - center1.v[1]) * (i1 - center1.v[1]) +
                                           (i2 - center1.v[2]) * (i2 - center1.v[2]));
                double d2 = radius2 - sqrt((i0 - center2.v[0]) * (i0 - center2.v[0]) +
                                           (i1 - center2.v[1]) * (i1 - center2.v[1]) +
                                           (i2 - center2.v[2]) * (i2 - center2.v[2]));

                double d = d1 > d2 ? d1 : d2;

                levelset(i0, i1, i2) = static_cast<_Ty>(d);

            } // i0
        } // i1
    } // i2
}

///////////////////////////////////////////////////////////////////////
///  Initialize a levelset function representing three spheres.
///
///  @param levelset  The level set function to be initialized.
///  @param center1   The center of the first  sphere.
///  @param radius1   The radius of the first  sphere.
///  @param center2   The center of the second sphere.
///  @param radius2   The radius of the second sphere.
///  @param center3   The center of the third  sphere.
///  @param radius3   The radius of the third  sphere.
///
///////////////////////////////////////////////////////////////////////
template <typename _Ty>
static void levelset3_init_three_spheres(array_base<_Ty>&      levelset,
                                         const point3<double>& center1,
                                         double                radius1,
                                         const point3<double>& center2,
                                         double                radius2,
                                         const point3<double>& center3,
                                         double                radius3
                                         )
{
    size_t i0, i1, i2;

    for (i2 = 0; i2 < levelset.size_2(); ++i2) {
        for (i1 = 0; i1 < levelset.size_1(); ++i1) {
            for (i0 = 0; i0 < levelset.size_0(); ++i0) {

                double d1 = radius1 - sqrt((i0 - center1.v[0]) * (i0 - center1.v[0]) +
                                           (i1 - center1.v[1]) * (i1 - center1.v[1]) +
                                           (i2 - center1.v[2]) * (i2 - center1.v[2]));
                double d2 = radius2 - sqrt((i0 - center2.v[0]) * (i0 - center2.v[0]) +
                                           (i1 - center2.v[1]) * (i1 - center2.v[1]) +
                                           (i2 - center2.v[2]) * (i2 - center2.v[2]));
                double d3 = radius3 - sqrt((i0 - center3.v[0]) * (i0 - center3.v[0]) +
                                           (i1 - center3.v[1]) * (i1 - center3.v[1]) +
                                           (i2 - center3.v[2]) * (i2 - center3.v[2]));

                double d = d1 > d2 ? d1 : d2;
                if (d3 > d) d = d3;

                levelset(i0, i1, i2) = static_cast<_Ty>(d);

            } // i0
        } // i1
    } // i2
}

struct sphere_info
{
    point3<double>  center;
    double          radius;
};

///////////////////////////////////////////////////////////////////////
///  Initialize a levelset function representing spheres.
///
///  @param levelset     The level set function to be initialized.
///  @param sphere_infos A vector containing the infomation for each
///                      sphere.
///
///////////////////////////////////////////////////////////////////////
template <typename _Ty>
static void levelset3_init_spheres(
    array_base<_Ty>&                levelset,
    const std::vector<sphere_info>& sphere_infos
    )
{
    size_t              i0, i1, i2, i;
    std::vector<double> distances(sphere_infos.size());

    for (i2 = 0; i2 < levelset.size_2(); ++i2) {
        for (i1 = 0; i1 < levelset.size_1(); ++i1) {
            for (i0 = 0; i0 < levelset.size_0(); ++i0) {
                
                for (i = 0; i < sphere_infos.size(); ++i) {
                    const sphere_info& si = sphere_infos[i];

                    distances[i] = si.radius - 
                        sqrt((i0 - si.center.v[0]) * (i0 - si.center.v[0]) +
                             (i1 - si.center.v[1]) * (i1 - si.center.v[1]) +
                             (i2 - si.center.v[2]) * (i2 - si.center.v[2]));
                }

                double d = *std::max_element(distances.begin(), distances.end());
                
                levelset(i0, i1, i2) = static_cast<_Ty>(d);
            } // i0
        } // i1
    } // i2
}

} // namespace

#endif // ___LEVELSET3_INITIALIZER_HXX___
