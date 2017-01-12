/*
 ==========================================================================
 |   
 |   $Id: arith.hxx 126 2005-02-05 18:06:44Z kangli $
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

#ifndef ___ARITH_HXX___
#   define ___ARITH_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif

#   include <cassert>
#   include <cmath>

/// @namespace optnet
namespace optnet { 
    /// @namespace optnet::utils
    namespace utils {

///////////////////////////////////////////////////////////////////////////
/// Determines if a number is power-of-two.
///
/// @return Returns true of the number is power-of-two, false otherwise.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
bool is_power_of_two(_T x)
{
    return (x != 0) && ((x & -x) == x);
}

///////////////////////////////////////////////////////////////////////////
/// Returns the nearest power-of-two number that is larget or equal to
/// the given number.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
_T ceil_power_of_two(_T x)
{
    if (x < 1) return 1;
    if (is_power_of_two(x)) return x;
    double n = ceil(log(double(x)) / log(2.0));
    return static_cast<_T>(pow(2.0, n));
}

    } // namespace
} // namespace

#endif
