/*
 ==========================================================================
 |   
 |   $Id: levelset2_initializer.hxx 164 2005-02-10 09:31:46Z kangli $
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

#ifndef ___LEVELSET2_INITIALIZER_HXX___
#   define ___LEVELSET2_INITIALIZER_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/array2.hxx>
#   include <optnet/_base/array2_ref.hxx>
#   include <optnet/_base/point2.hxx>
#   include <cmath>

namespace optnet {

///////////////////////////////////////////////////////////////////////
///  Initialize a levelset function representing a circle.
///
///  @param levelset  The level set function to be initialized.
///  @param center    The center of the circle.
///  @param radius    The radius of the circle.
///
///////////////////////////////////////////////////////////////////////
template <typename _Ty>
static void levelset2_init_circle(array2_base<_Ty>&     levelset,
                                  const point2<double>& center,
                                  double                radius
                                  )
{
    size_t i0, i1;

    for (i1 = 0; i1 < levelset.size_1(); ++i1) {
        for (i0 = 0; i0 < levelset.size_0(); ++i0) {

            double d = sqrt((i0 - center.v[0]) * (i0 - center.v[0]) +
                            (i1 - center.v[1]) * (i1 - center.v[1]));
            levelset(i0, i1) = static_cast<_Ty>(radius - d);

        } // i0
    } // i1
}

///////////////////////////////////////////////////////////////////////
///  Initialize a levelset function representing two circles.
///
///  @param levelset  The level set function to be initialized.
///  @param center1   The center of the first  circle.
///  @param radius1   The radius of the first  circle.
///  @param center2   The center of the second circle.
///  @param radius2   The radius of the second circle.
///
///////////////////////////////////////////////////////////////////////
template <typename _Ty>
static void levelset2_init_two_circles(array2_base<_Ty>&     levelset,
                                       const point2<double>& center1,
                                       double                radius1,
                                       const point2<double>& center2,
                                       double                radius2
                                       )
{
    size_t i0, i1;

    for (i1 = 0; i1 < levelset.size_1(); ++i1) {
        for (i0 = 0; i0 < levelset.size_0(); ++i0) {

            double d1 = radius1 - 
                sqrt((i0 - center1.v[0]) * (i0 - center1.v[0]) +
                     (i1 - center1.v[1]) * (i1 - center1.v[1]));
            double d2 = radius2 - 
                sqrt((i0 - center2.v[0]) * (i0 - center2.v[0]) +
                     (i1 - center2.v[1]) * (i1 - center2.v[1]));

            double d = d1 > d2 ? d1 : d2;

            levelset(i0, i1) = static_cast<_Ty>(d);

        } // i0
    } // i1
}

struct circle_info
{
    point2<double>  center;
    double          radius;
};

///////////////////////////////////////////////////////////////////////
///  Initialize a levelset function representing circles.
///
///  @param levelset     The level set function to be initialized.
///  @param circle_infos A vector containing the infomation for each
///                      circle.
///
///////////////////////////////////////////////////////////////////////
template <typename _Ty>
static void levelset2_init_circles(
    array2_base<_Ty>&               levelset,
    const std::vector<circle_info>& circle_infos
    )
{
    size_t              i0, i1, i;
    std::vector<double> distances(circle_infos.size());

    for (i1 = 0; i1 < levelset.size_1(); ++i1) {
        for (i0 = 0; i0 < levelset.size_0(); ++i0) {
            
            for (i = 0; i < circle_infos.size(); ++i) {
                const circle_info& ci = circle_infos[i];

                distances[i] = ci.radius - 
                    sqrt((i0 - ci.center.v[0]) * (i0 - ci.center.v[0]) +
                         (i1 - ci.center.v[1]) * (i1 - ci.center.v[1]));
            }

            double d = *std::max_element(distances.begin(), distances.end());
            
            levelset(i0, i1) = static_cast<_Ty>(d);
        } // i0
    } // i1
}

} // namespace

#endif // ___LEVELSET2_INITIALIZER_HXX___
