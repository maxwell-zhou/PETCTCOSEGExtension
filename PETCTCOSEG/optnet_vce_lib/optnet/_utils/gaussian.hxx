/* 
 ==========================================================================
 |   
 |   $Id: gaussian.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___GAUSSIAN_HXX___
#   define ___GAUSSIAN_HXX___

#   include <cmath>

/// @namespace optnet
namespace optnet {
    /// @namespace optnet::utils
    namespace utils {

///////////////////////////////////////////////////////////////////////////
///  Computes the (zero-mean) gaussian function.
///////////////////////////////////////////////////////////////////////////
template<typename _T>
_T
gaussian(_T x, _T sigma) {
    static _T   _2pi = static_cast<_T>(8.0 * atan(1.0));
    return
        (
            static_cast<_T>
                (exp (-(x*x) / (2.0*sigma*sigma)) / (sqrt (_2pi)*sigma))
        );
}


///////////////////////////////////////////////////////////////////////////
///  Computes the first derivative of the (zero-mean) gaussian function.
///////////////////////////////////////////////////////////////////////////
template<typename _T>
_T
dgaussian(_T x, _T sigma) {
    static _T   _2pi = static_cast<_T>(8.0 * atan(1.0));
    return
        (
            static_cast<_T>
                (
                    -x*exp (-(x*x) / (2.0*sigma*sigma)) /
                        (sqrt (_2pi)*sigma*sigma*sigma)
                )
        );
}


///////////////////////////////////////////////////////////////////////////
///  Computes the second derivative of the (zero-mean) gaussian function.
///////////////////////////////////////////////////////////////////////////
template<typename _T>
_T
d2gaussian(_T x, _T sigma) {
    static _T   _2pi = static_cast<_T>(8.0 * atan(1.0));
    return
        (
            static_cast<_T>
                (
                    ((x*x) - (sigma*sigma))*exp (
      -(x*x) / (2.0*sigma*sigma)) /
                    (sqrt (_2pi)*sigma*sigma*sigma*sigma*sigma)
                )
        );
}

    }   // namespace
}   // namespace

#endif // ___GAUSSIAN_HXX__
