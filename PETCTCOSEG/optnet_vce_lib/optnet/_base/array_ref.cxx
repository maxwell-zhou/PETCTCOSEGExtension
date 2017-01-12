/*
 ==========================================================================
 |   
 |   $Id: array_ref.cxx 58 2005-01-21 01:44:10Z kangli $
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

#ifndef ___ARRAY_REF_CXX___
#   define ___ARRAY_REF_CXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4244)
#   endif

#   include <optnet/_base/array_ref.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array_ref<_Ty, _Tg>::array_ref()
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array_ref<_Ty, _Tg>::array_ref(_Ty* p, 
                               size_type s0, 
                               size_type s1, 
                               size_type s2, 
                               size_type s3, 
                               size_type s4
                               )
{
    assert(p != 0);
    assert(s0 > 0);
    assert(s1 > 0);
    assert(s2 > 0);
    assert(s3 > 0);
    assert(s4 > 0);

    _Base::m_p      = p;
    _Base::m_sz     = s0 * s1 * s2 * s3 * s4;
    _Base::m_asz[0] = s0;
    _Base::m_asz[1] = s1;
    _Base::m_asz[2] = s2;
    _Base::m_asz[3] = s3;
    _Base::m_asz[4] = s4;

    _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array_ref<_Ty, _Tg>::array_ref(const array_ref<_Ty, _Tg>& robj)
{
    assert(&robj != this);

    _Base::m_p  = robj.m_p;
    _Base::m_sz = robj.m_sz;
    ::memcpy(_Base::m_asz, 
             robj.m_asz, 
             OPTNET_ARRAY_DIMS * sizeof(size_type));

    // Initialize the look-up-table.
    _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array_ref<_Ty, _Tg>::~array_ref()
{
    _Base::free_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array_ref<_Ty, _Tg>&
array_ref<_Ty, _Tg>::operator=(const array_ref<_Ty, _Tg>& robj)
{
    if (&robj != this) {

        _Base::m_p  = robj.m_p;
        _Base::m_sz = robj.m_sz;
        ::memcpy(_Base::m_asz, 
                 robj.m_asz, 
                 OPTNET_ARRAY_DIMS * sizeof(size_type));

        // Reinitialize the look-up-table.
        _Base::free_lut();
        _Base::init_lut();
    }
    return *this;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array_ref<_Ty, _Tg>::assign (const _Ty* p,
                             size_type s0, 
                             size_type s1, 
                             size_type s2, 
                             size_type s3, 
                             size_type s4
                             )
{
    if (0 == p) return false;

    _Base::m_p = p;

    if (_Base::m_asz[0] != s0 ||
        _Base::m_asz[1] != s1 || 
        _Base::m_asz[2] != s2 ||
        _Base::m_asz[3] != s3 || 
        _Base::m_asz[4] != s4) {

        _Base::free_lut();

        _Base::m_sz     = s0 * s1 * s2 * s3 * s4;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
        _Base::m_asz[2] = s2;
        _Base::m_asz[3] = s3;
        _Base::m_asz[4] = s4;
        
        return _Base::init_lut();
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array_ref<_Ty, _Tg>::reshape(size_type s0, 
                             size_type s1, 
                             size_type s2, 
                             size_type s3, 
                             size_type s4
                             )
{
    assert(s0 > 0);
    assert(s1 > 0);
    assert(s2 > 0);
    assert(s3 > 0);
    assert(s4 > 0);

    if (0 == _Base::m_p) return false;

    size_type sz = s0 * s1 * s2 * s3 * s4;

    // The array sizes must be equal.
    if (_Base::m_sz == sz) {

        _Base::free_lut();

        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
        _Base::m_asz[2] = s2;
        _Base::m_asz[3] = s3;
        _Base::m_asz[4] = s4;

        return _Base::init_lut();
    }

    return false;
}

} // namespace

#endif
