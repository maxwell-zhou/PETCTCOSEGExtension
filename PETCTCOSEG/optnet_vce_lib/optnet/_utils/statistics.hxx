/*
 ==========================================================================
 |   
 |   $Id: statistics.hxx 127 2005-02-06 00:46:11Z kangli $
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

#ifndef ___STATISTICS_HXX___
#   define ___STATISTICS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif

#   include <cmath>
#   include <cassert>
#   include <optnet/define.h>

/// @namespace optnet
namespace optnet { 
    /// @namespace optnet::utils
    namespace utils {

///////////////////////////////////////////////////////////////////////////
///  Normal probability density function (pdf).
///////////////////////////////////////////////////////////////////////////
template <typename _T>
inline double normpdf(const _T& x, double mu, double sigma)
{
    return exp(-(x - mu) * (x - mu) * 0.5 / (sigma * sigma)) /
        (sigma * M_SQRTTWOPI);
}

///////////////////////////////////////////////////////////////////////////
///  Rayleigh probability density function.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
inline double raylpdf(const _T& x, double b)
{
    return x * exp(-x * x * 0.5 / (b * b)) / (b * b);
}

///////////////////////////////////////////////////////////////////////////
///  Rayleigh cumulative density function.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
inline double raylcdf(const _T& x, double b)
{
    return 1 - exp(-0.5 * sqrt(x / b));
}

    } // namespace
} // namespace

#endif
