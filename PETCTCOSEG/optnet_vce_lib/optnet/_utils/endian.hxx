/*
 ==========================================================================
 |   
 |   $Id: endian.hxx 126 2005-02-05 18:06:44Z kangli $
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

#ifndef ___ENDIAN_HXX___
#   define ___ENDIAN_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <algorithm>

#   include <cassert>


/// @namespace optnet
namespace optnet { 
    /// @namespace optnet::utils
    namespace utils {

///////////////////////////////////////////////////////////////////////////
///  Adjust the endian of a 32-bit value.
///
///  @param[in,out] p  A pointer to the 32-bit value.
///
///////////////////////////////////////////////////////////////////////////
template<typename _T>
inline void
swap_endian32(_T* p) {
    assert(sizeof(_T) == 4);
    std::swap(*((char*) p), *((char*) p + 3));
    std::swap(*((char*) p + 1), *((char*) p + 2));
}

///////////////////////////////////////////////////////////////////////////
///  Adjust the endian of a 16-bit value.
///
///  @param[in,out] p  A pointer to the 16-bit value.
///
///////////////////////////////////////////////////////////////////////////
template<typename _T>
inline void
swap_endian16(_T* p) {
    assert(sizeof(_T) == 2);
    std::swap(*((char*) p), *((char*) p + 1));
}

    } // namespace
} // namespace

#endif
