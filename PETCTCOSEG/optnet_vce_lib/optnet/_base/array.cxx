/*
 ==========================================================================
 |   
 |   $Id: array.cxx 146 2005-02-08 05:26:28Z kangli $
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

#ifndef ___ARRAY_CXX___
#   define ___ARRAY_CXX___

#   include <optnet/_base/array.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4244)
#       if (_MSC_VER <= 1200)
#           pragma warning(disable: 4018)
#           pragma warning(disable: 4146)
#       endif
#   endif

#   include <algorithm>


namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array() :
    array_base<_Ty, _Tg>()
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array(const size_type* asz)
{
    assert(asz[0] > 0);
    assert(asz[1] > 0);
    assert(asz[2] > 0);
    assert(asz[3] > 0);
    assert(asz[4] > 0);

    size_type sz = asz[0] * asz[1] * asz[2] * asz[3] * asz[4];

    _Base::m_p = new value_type[sz];

    // Note: see above.
    assert(0 != _Base::m_p);

    std::fill(_Base::m_p, _Base::m_p + sz, value_type());

    _Base::m_sz     = sz;
    _Base::m_asz[0] = asz[0];
    _Base::m_asz[1] = asz[1];
    _Base::m_asz[2] = asz[2];
    _Base::m_asz[3] = asz[3];
    _Base::m_asz[4] = asz[4];

    _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array(size_type s0, 
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

    size_type sz = s0 * s1 * s2 * s3 * s4;

    _Base::m_p = new value_type[sz];

    // Note: see above.
    assert(0 != _Base::m_p);

    std::fill(_Base::m_p, _Base::m_p + sz, value_type());

    _Base::m_sz     = sz;
    _Base::m_asz[0] = s0;
    _Base::m_asz[1] = s1;
    _Base::m_asz[2] = s2;
    _Base::m_asz[3] = s3;
    _Base::m_asz[4] = s4;

    _Base::init_lut();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array(const array& robj)
{
    assert(&robj != this);

    if (0 != robj.m_p && 0 != robj.m_sz) {

        _Base::m_p = new value_type[robj.m_sz];
        
        // This is necessary for non-standard compliant compilers.
        assert(0 != _Base::m_p);

        // Perform deep copy. Use std::copy since the elements may not
        // be simple type.
        std::copy(robj.begin(), robj.end(), this->begin());

        _Base::m_sz = robj.size();
        ::memcpy(_Base::m_asz, 
                 robj.m_asz, 
                 OPTNET_ARRAY_DIMS * sizeof(size_type));

        // Initialize the look-up-table.
        _Base::init_lut();

    } // if
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array(const array_base<_Ty, _Tg>& robj)
{
    assert(&robj != this);

    if (0 != robj.data() && 0 != robj.size()) {

        _Base::m_p = new value_type[robj.size()];
        
        // This is necessary for non-standard compliant compilers.
        assert(0 != _Base::m_p);

        // Perform deep copy. Use std::copy since the elements may not
        // be simple type.
        std::copy(robj.begin(), robj.end(), this->begin());

        _Base::m_sz = robj.size();
        ::memcpy(_Base::m_asz, 
                 robj.sizes(), 
                 OPTNET_ARRAY_DIMS * sizeof(size_type));

        // Initialize the look-up-table.
        _Base::init_lut();

    } // if
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::array(const array2_base<_Ty>& robj)
{
    if (0 != robj.data() && 0 != robj.size()) {

        _Base::m_p = new value_type[robj.size()];
        
        // This is necessary for non-standard compliant compilers.
        assert(0 != _Base::m_p);

        // Perform deep copy. Use std::copy since the elements may not
        // be simple type.
        std::copy(robj.begin(), robj.end(), this->begin());

        _Base::m_sz = robj.size();
        ::memcpy(_Base::m_asz, robj.sizes(), 2 * sizeof(size_type));

        _Base::m_asz[2] = 1;
        _Base::m_asz[3] = 1;
        _Base::m_asz[4] = 1;

        // Initialize the look-up-table.
        _Base::init_lut();

    } // if
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>::~array()
{
    release();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>&
array<_Ty, _Tg>::operator=(const array& robj)
{
    if (&robj != this && 0 != robj.data() && 0 != robj.size()) {

        if (_Base::m_sz != robj.size()) {

            SAFE_DELETE_ARRAY(_Base::m_p);

            _Base::m_p = new value_type[robj.size()];
            assert(0 != _Base::m_p);

            _Base::m_sz = robj.size();
            memcpy(_Base::m_asz, 
                   robj.sizes(), 
                   OPTNET_ARRAY_DIMS * sizeof(size_type));

            // Reinitialize the look-up-table.
            _Base::free_lut();
            _Base::init_lut();
        }

        // Perform deep copy. Use std::copy since the elements may not
        // be simple type.
        std::copy(robj.begin(), robj.end(), this->begin());

    } // if

    return *this;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>&
array<_Ty, _Tg>::operator=(const array_base<_Ty, _Tg>& robj)
{
    if (&robj != this && 0 != robj.data() && 0 != robj.size()) {

        if (_Base::m_sz != robj.size()) {

            SAFE_DELETE_ARRAY(_Base::m_p);
            
            _Base::m_p = new value_type[robj.size()];
            assert(0 != _Base::m_p);

            _Base::m_sz = robj.size();
            memcpy(_Base::m_asz, 
                   robj.sizes(), 
                   OPTNET_ARRAY_DIMS * sizeof(size_type));

            // Reinitialize the look-up-table.
            _Base::free_lut();
            _Base::init_lut();
        }

        // Perform deep copy. Use std::copy since the elements may not
        // be simple type.
        std::copy(robj.begin(), robj.end(), this->begin());

    } // if

    return *this;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
array<_Ty, _Tg>&
array<_Ty, _Tg>::operator=(const array2_base<_Ty>& robj)
{
    if (0 != robj.data() && 0 != robj.size()) {

        if (_Base::m_sz != robj.size()) {

            SAFE_DELETE_ARRAY(_Base::m_p);
            
            _Base::m_p = new value_type[robj.size()];
            assert(0 != _Base::m_p);

            memcpy(_Base::m_asz, robj.sizes(), 2 * sizeof(size_type));

            _Base::m_sz = robj.size();
            _Base::m_asz[2] = 1;
            _Base::m_asz[3] = 1;
            _Base::m_asz[4] = 1;

            // Reinitialize the look-up-table.
            _Base::free_lut();
            _Base::init_lut();
        }

        // Perform deep copy. Use std::copy since the elements may not
        // be simple type.
        std::copy(robj.begin(), robj.end(), this->begin());

    } // if

    return *this;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array<_Ty, _Tg>::create (const size_type* asz)
{
    size_type sz = asz[0] * asz[1] * asz[2] * asz[3] * asz[4];

    if (0 == sz) {
        assert(false); return false;
    }    

    // Relocate memory only if size is changed.
    if (_Base::m_sz != sz) {

        // Free allocated resources first.
        SAFE_DELETE_ARRAY(_Base::m_p);

        _Base::m_p = new value_type[sz];
        if (0 == _Base::m_p) {
            assert(false); return false;
        }

        _Base::m_sz     = sz;
        _Base::m_asz[0] = asz[0];
        _Base::m_asz[1] = asz[1];
        _Base::m_asz[2] = asz[2];
        _Base::m_asz[3] = asz[3];
        _Base::m_asz[4] = asz[4];

        // Reinitialize the look-up-table.
        _Base::free_lut();
        _Base::init_lut();
    }
    
    std::fill(_Base::m_p, _Base::m_p + sz, value_type());

    return true;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array<_Ty, _Tg>::create (size_type s0, 
                         size_type s1, 
                         size_type s2, 
                         size_type s3,
                         size_type s4
                         )
{
    size_type sz = s0 * s1 * s2 * s3 * s4;

    if (0 == sz) {
        assert(false); return false;
    }    

    // Relocate memory only if size is changed.
    if (_Base::m_sz != sz) {

        // Free allocated resources first.
        SAFE_DELETE_ARRAY(_Base::m_p);

        _Base::m_p = new value_type[sz];
        if (0 == _Base::m_p) {
            assert(false); return false;
        }

        _Base::m_sz     = sz;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
        _Base::m_asz[2] = s2;
        _Base::m_asz[3] = s3;
        _Base::m_asz[4] = s4;

        // Reinitialize the look-up-table.
        _Base::free_lut();
        _Base::init_lut();
    }
    
    std::fill(_Base::m_p, _Base::m_p + sz, value_type());

    return true;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array<_Ty, _Tg>::create_and_fill(const value_type& value,
                                 size_type         s0, 
                                 size_type         s1, 
                                 size_type         s2, 
                                 size_type         s3,
                                 size_type         s4
                                 )
{
    size_type sz = s0 * s1 * s2 * s3 * s4;

    if (0 == sz) {
        assert(false); return false;
    }

    // Relocate memory only if size is changed.
    if (_Base::m_sz != sz) {

        // Free allocated resources first.
        SAFE_DELETE_ARRAY(_Base::m_p);

        _Base::m_p = new value_type[sz];
        if (0 == _Base::m_p) {
            assert(false); return false;
        }
        
        _Base::m_sz     = sz;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
        _Base::m_asz[2] = s2;
        _Base::m_asz[3] = s3;
        _Base::m_asz[4] = s4;

        // Reinitialize the look-up-table.
        _Base::free_lut();
        _Base::init_lut();
    }
    
    std::fill(_Base::m_p, _Base::m_p + sz, value);

    return true;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
void
array<_Ty, _Tg>::release()
{
    // Free array data buffer.
    SAFE_DELETE_ARRAY(_Base::m_p);

    // Free look-up-table.
    _Base::free_lut();
    
    // Reset sizes.
    _Base::m_sz     = 0;
    _Base::m_asz[0] = 0;
    _Base::m_asz[1] = 0;
    _Base::m_asz[2] = 0;
    _Base::m_asz[3] = 0;
    _Base::m_asz[4] = 0;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg>
bool
array<_Ty, _Tg>::reshape(size_type s0, 
                         size_type s1, 
                         size_type s2, 
                         size_type s3,
                         size_type s4
                         )
{
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
