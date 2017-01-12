/*
 ==========================================================================
 |   
 |   $Id: cubic_spline.hxx 186 2005-02-27 02:04:22Z kangli $
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

#ifndef ___CUBIC_SPLINE_HXX___
#   define ___CUBIC_SPLINE_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/consts.hxx>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class cubic_spline
///  @brief Cubic spline evaluator.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
struct cubic_spline
{
    typedef _T                  value_type;
    typedef value_type&         reference;
    typedef const value_type&   const_reference;
    typedef size_t              size_type;

    ///////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    ///////////////////////////////////////////////////////////////////////
    cubic_spline()
    {
        b[0] = b[1] = b[2] = b[3] = 0.0;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Constructor.
    ///  @param a  The value where the cubic-spline basis is evaluated.
    ///////////////////////////////////////////////////////////////////////
    cubic_spline(const_reference a)             { evaluate(a); }

    ///////////////////////////////////////////////////////////////////////
    ///  @param a  The value where the cubic-spline basis is evaluated.
    ///////////////////////////////////////////////////////////////////////
    inline void operator()(const_reference a)   { evaluate(a); }

    ///////////////////////////////////////////////////////////////////////
    ///  @param a  The value where the cubic-spline basis is evaluated.
    ///////////////////////////////////////////////////////////////////////
    inline void evaluate(const_reference a)
    {
        _T __a2  = a *    a;
        _T __a3  = a * __a2;
        _T __3a  = 3 *    a;
        _T __3a2 = 3 * __a2;
        _T __3a3 = 3 * __a3;

        b[0] = (1.0 - __3a  + __3a2 - __a3 ) * INV6;
        b[1] = (4.0 + __3a3 -   2.0 * __3a2) * INV6;
        b[2] = (1.0 + __3a  + __3a2 - __3a3) * INV6;
        b[3] = __a3 * INV6;
    }   

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the value of the i-th basis function.
    ///
    ///  @param i  The index of the basis function.
    ///////////////////////////////////////////////////////////////////////
    inline const_reference operator[](size_type i) const
    {
        return b[i];
    }

    _T b[4]; /// The evaluated cubic-spline basis function.
};

} // namespace

#endif // ___CUBIC_SPLINE_HXX___
