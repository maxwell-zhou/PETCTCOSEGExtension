/* 
 ==========================================================================
 |   
 |   $Id: regiongrow.hxx 21 2005-01-14 15:52:31Z kangli $
 |
 |   Written by Kang Li <kangl@cmu.edu>
 |   College of Engineering Imaging Group
 |   The University of Iowa
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
 
#ifndef ___REGIONGROW_HXX___
#   define ___REGIONGROW_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#      pragma once
#      pragma warning(disable : 4786)
#   endif

#   include <queue>
#   include <optnet/config.h>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_base/triple.hxx>

/// @namespace optnet
namespace optnet {
    /// @namespace optnet::utils
    namespace utils {

///////////////////////////////////////////////////////////////////////////
///  @class regiongrow_3d_ongrow
///  @brief A functor called during the process of region-grow.
///////////////////////////////////////////////////////////////////////////
template <typename _Tp, typename _Tg = net_f_xy>
struct regiongrow_3d_ongrow
{
    virtual bool operator()(optnet::array_base<_Tp, _Tg>&,
                            size_t, size_t, size_t
                            ) = 0;
};

///////////////////////////////////////////////////////////////////////////
///  3-D region-grow algorithm framework (6-neighbor configuration).
///////////////////////////////////////////////////////////////////////////
template <typename _Tp, typename _Tg>
void
regiongrow_3d_6neighbor(optnet::array_base<_Tp, _Tg>&   a,
                        size_t                          seed_x,
                        size_t                          seed_y,
                        size_t                          seed_z,
                        regiongrow_3d_ongrow<_Tp, _Tg>& ongrow,
                        size_t                          min_x = -1,
                        size_t                          max_x = -1,
                        size_t                          min_y = -1,
                        size_t                          max_y = -1,
                        size_t                          min_z = -1,
                        size_t                          max_z = -1
                        )
{
    typedef optnet::triple<size_t, size_t, size_t>   point_type;

    assert(a.size_0() > 0 && a.size_1() > 0 && a.size_2() > 0);
    assert(seed_x < a.size_0());
    assert(seed_y < a.size_1());
    assert(seed_z < a.size_2());

    if (min_x > (a.size_0() - 1)) min_x = 0;
    if (min_y > (a.size_1() - 1)) min_y = 0;
    if (min_z > (a.size_2() - 1)) min_z = 0;
    
    if (max_x > (a.size_0() - 1)) max_x = a.size_0() - 1;
    if (max_y > (a.size_1() - 1)) max_y = a.size_1() - 1;
    if (max_z > (a.size_2() - 1)) max_z = a.size_2() - 1;
    
    
    std::queue<point_type>   Q;
    point_type               c, n;
    
    // Check if the seed point satisfies the grow condition.
    if (!ongrow(a, seed_x, seed_y, seed_z)) return;
    
    // Enqueue the seed point.
    Q.push(point_type(seed_x, seed_y, seed_z));
    
    while (!Q.empty())
    {
        c = Q.front();
        Q.pop();

        if (c.first  > min_x)  // left
        {
            n = point_type(c.first - 1, c.second, c.third);
            if (ongrow(a, n.first, n.second, n.third)) {
                Q.push(n);
            }
        }

        if (c.first  < max_x)  // right
        {
            n = point_type(c.first + 1, c.second, c.third);
            if (ongrow(a, n.first, n.second, n.third)) {
                Q.push(n);
            }
        }

        if (c.second > min_y)  // above
        {
            n = point_type(c.first, c.second - 1, c.third);
            if (ongrow(a, n.first, n.second, n.third)) {
                Q.push(n);
            }
        }

        if (c.second < max_y)  // below
        {
            n = point_type(c.first, c.second + 1, c.third);
            if (ongrow(a, n.first, n.second, n.third)) {
                Q.push(n);
            }
        }

        if (c.third  > min_z)  // front
        {
            n = point_type(c.first, c.second, c.third - 1);
            if (ongrow(a, n.first, n.second, n.third)) {
                Q.push(n);
            }
        }

        if (c.third  < max_z)  // back
        {
            n = point_type(c.first, c.second, c.third + 1);
            if (ongrow(a, n.first, n.second, n.third)) {
                Q.push(n);
            }
        }

    } // while
}
    
    
    } // namespace
} // namespace

#endif // ___REGIONGROW_HXX___

