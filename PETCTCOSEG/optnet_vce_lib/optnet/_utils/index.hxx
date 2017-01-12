/* 
 ==========================================================================
 |   
 |   $Id: index.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___INDEX_HXX___
#   define ___INDEX_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>
#   include <cassert>

/// @namespace optnet
namespace optnet {
    /// @namespace optnet::utils
    namespace utils {

///////////////////////////////////////////////////////////////////////////
///
///  @return  The converted index under the circularity criterion.
///           The range of the returned index is in [0..imax-1].
///////////////////////////////////////////////////////////////////////////
template <typename _Tp>
inline _Tp
index_circular   (const _Tp& i, const _Tp& imax)
{
    _Tp j = i;
    if (j >= imax) j -= imax;
    else if (j < 0) j += imax;
    return j;
}

///////////////////////////////////////////////////////////////////////////
///
///  @return  The converted index under the circularity criterion.
///           The range of the returned index is in [0..imax-1].
///////////////////////////////////////////////////////////////////////////
template <typename _Tp>
inline _Tp
index_circular_ex(const _Tp& i, const _Tp& imax)
{
    _Tp j = i % imax;
    if (j < 0) j = imax + j;
    return j;
}

///////////////////////////////////////////////////////////////////////////
///
///  @return  The converted index under the symmetric (mirror) criterion.
///           The range of the returned index is in [0..imax-1].
///////////////////////////////////////////////////////////////////////////
template <typename _Tp>
inline _Tp
index_symmetric   (const _Tp& i, const _Tp& imax)
{
    _Tp j = i;
    if (j >= imax) j = imax - (j - imax) - 2;
    else if (j < 0) j = -j;
    return j;
}

///////////////////////////////////////////////////////////////////////////
///
///  @return  The converted index under the symmetric (mirror) criterion.
///           The range of the returned index is in [0..imax-1].
///////////////////////////////////////////////////////////////////////////
template <typename _Tp>
inline _Tp
index_symmetric_ex(const _Tp& i, const _Tp& imax)
{
    assert(false); // Not implemented yet.
    return i;
}

    } // namespace
} // namespace

#endif // ___INDEX_HXX___
