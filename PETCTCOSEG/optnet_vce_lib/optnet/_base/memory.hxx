/*
 ==========================================================================
 |   
 |   $Id: memory.hxx 116 2005-02-04 06:17:04Z kangli $
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

#ifndef ___MEMORY_HXX___
#   define ___MEMORY_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <cstddef>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
/// Allocate an m-by-n (row major) matrix.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
_T** new_matrix(size_t m, size_t n)
{
    _T** pp = new _T*[m];
    *pp     = new _T [m * n];
    for (_T** pp1 = pp; pp1 < pp + m - 1; ++pp1)
        *(pp1 + 1) = *pp1 + n;
    return pp;
}

///////////////////////////////////////////////////////////////////////////
/// Allocate an m-by-n (row major) matrix with initial value.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
_T** new_matrix(size_t m, size_t n, const _T& v)
{
    size_t i;
    _T** pp = new _T*[m];
    *pp     = new _T [m * n];
    for (i = 0; i < m * n; ++i) (*pp)[i] = v;
    for (i = 0; i + 1 < m; ++i) pp[i + 1] = pp[i] + n;
    return pp;
}

///////////////////////////////////////////////////////////////////////////
/// Free a 2-D matrix.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
inline void delete_matrix(_T** pp)
{
    delete[] *pp;
    delete[]  pp;
}

} // namespace

///////////////////////////////////////////////////////////////////////////
//  macros

#   ifndef SAFE_DELETE
/// Delete a variable pointed by p, and set p to NULL.
#       define SAFE_DELETE(p)       { if (p) { delete   (p); (p) = 0; } }
#   endif // SAFE_DELETE

#   ifndef SAFE_DELETE_ARRAY
/// Delete an array   pointed by p, and set p to NULL.
#       define SAFE_DELETE_ARRAY(p) { if (p) { delete[] (p); (p) = 0; } }
#   endif // SAFE_DELETE_ARRAY

#   ifndef SAFE_DELETE_MATRIX
/// Delete a matrix pointed by pp, and set pp to NULL.
#       define SAFE_DELETE_MATRIX(pp) \
                 { if (pp) { optnet::delete_matrix(pp); (pp) = 0; } }
#   endif // SAFE_DELETE_MATRIX

#endif 
