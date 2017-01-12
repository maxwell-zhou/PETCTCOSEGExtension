/*
 ==========================================================================
 |   
 |   $Id: array2.hxx 41 2005-01-18 21:23:26Z kangli $
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

#ifndef ___ARRAY2_HXX___
#   define ___ARRAY2_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/_base/array2_base.hxx>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class array2 array2.hxx "_base/array2.hxx"
///  @brief Array template class.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty>
class array2: public array2_base<_Ty>
{
    typedef array2_base<_Ty>    _Base;

public:

    typedef typename _Base::value_type      value_type;
    typedef typename _Base::reference       reference;
    typedef typename _Base::const_reference const_reference;
    typedef typename _Base::iterator        iterator;
    typedef typename _Base::const_iterator  const_iterator;
    typedef typename _Base::pointer         pointer;
    typedef typename _Base::const_pointer   const_pointer;
    typedef typename _Base::difference_type difference_type;
    typedef typename _Base::size_type       size_type;

    // constructor/destructors
    array2() :
        array2_base<_Ty>()
    {
    }
    
    array2(size_type s0, size_type s1)
    {
        assert(s0 > 0);
        assert(s1 > 0);
        
        size_type sz = s0 * s1;

        _Base::m_p = new value_type[sz];

        // This is necessary for non-standard compliant compilers.
        assert(0 != _Base::m_p);

        std::fill(_Base::m_p, _Base::m_p + sz, value_type());

        _Base::m_sz     = sz;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
   }
    
    array2(const array2& robj)
    {
        assert(&robj != this);

        if (0 != robj._Base::m_p && 0 != robj._Base::m_sz) {

            _Base::m_p  = new value_type[robj._Base::m_sz];
            
            // This is necessary for non-standard compliant compilers.
            assert(0 != _Base::m_p);

            // Perform deep copy. Use std::copy since the elements may not
            // be simple type.
            std::copy(robj.begin(), robj.end(), this->begin());

            _Base::m_asz[0] = robj._Base::m_asz[0];
            _Base::m_asz[1] = robj._Base::m_asz[1];
            _Base::m_sz     = robj._Base::m_sz;

        } // if
    }
    
    template <typename _Ta>
    array2(const array2_base<_Ta>& robj)
    {
        assert(&robj != this);

        if (0 != robj.data() && 0 != robj.size()) {

            // Free allocated resources first.
            SAFE_DELETE_ARRAY(_Base::m_p);

            _Base::m_p = new value_type[robj.size()];
            
            // This is necessary for non-standard compliant compilers.
            assert(0 != _Base::m_p);

            // Perform deep copy. Use std::copy since the elements may not
            // be simple type.
            std::copy(robj.begin(), robj.end(), this->begin());

            _Base::m_asz[0] = robj.sizes()[0];
            _Base::m_asz[1] = robj.sizes()[1];
            _Base::m_sz     = robj.size();

        } // if
    }
    
    ~array2()
    {
        release();
    }
    
    // modifiers
    array2& operator=(const array2& robj)
    {
        if (&robj != this && 0 != robj._Base::m_p && 0 != robj._Base::m_sz) {

            if ((_Base::m_sz != robj._Base::m_sz) || (0 == _Base::m_p)) {
                // Reallocate data.
                SAFE_DELETE_ARRAY(_Base::m_p);
                _Base::m_p = new value_type[robj._Base::m_sz];
            }

            // This is necessary for non-standard compliant compilers.
            assert(0 != _Base::m_p);

            // Perform deep copy. Use std::copy since the elements may not
            // be simple type.
            std::copy(robj.begin(), robj.end(), this->begin());
            
            _Base::m_asz[0] = robj._Base::m_asz[0];
            _Base::m_asz[1] = robj._Base::m_asz[1];
            _Base::m_sz     = robj._Base::m_sz;

        } // if

        return *this;
    }

    template <typename _Ta>
    array2& operator=(const array2_base<_Ta>& robj)
    {
        if (0 != robj.data() && 0 != robj.size()) {

            // Free allocated resources first.
            SAFE_DELETE_ARRAY(_Base::m_p);

            _Base::m_p = new value_type[robj.size()];
            
            // This is necessary for non-standard compliant compilers.
            assert(0 != _Base::m_p);

            // Perform deep copy. Use std::copy since the elements may not
            // be simple type.
            std::copy(robj.begin(), robj.end(), this->begin());

            _Base::m_asz[0] = robj.sizes()[0];
            _Base::m_asz[1] = robj.sizes()[1];
            _Base::m_sz     = robj.size();

        } // if

        return *this;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Creates an array of a specified size.
    ///
    ///  @param s0 The size of the first  dimension of the array.
    ///  @param s1 The size of the second dimension of the array.
    ///
    ///  @return Returns true if the array creation is successful.
    ///
    ///////////////////////////////////////////////////////////////////////
    bool create(size_type s0, size_type s1)
    {
        size_type sz = s0 * s1;

        if (0 == sz) return false;

        // Free allocated resources first.
        SAFE_DELETE_ARRAY(_Base::m_p);

        _Base::m_p = new value_type[sz];
        if (0 == _Base::m_p) 
            return false;

        std::fill(_Base::m_p, _Base::m_p + sz, value_type());

        _Base::m_sz     = sz;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;

        return true;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Creates an array of a specified size and fill its elements with
    ///  a specified value.
    ///
    ///  @param value The initial value of the array elements.
    ///  @param s0    The size of the first  dimension of the array.
    ///  @param s1    The size of the second dimension of the array.
    ///
    ///  @return Returns true if the array creation is successful.
    ///
    ///////////////////////////////////////////////////////////////////////
    bool create_and_fill(const value_type& value,
                         size_type         s0,
                         size_type         s1
                         )
    {
        size_type sz = s0 * s1;

        if (0 == sz) return false;

        // Free allocated resources first.
        SAFE_DELETE_ARRAY(_Base::m_p);

        _Base::m_p = new value_type[sz];
        if (0 == _Base::m_p) 
            return false;

        std::fill(_Base::m_p, _Base::m_p + sz, value);

        _Base::m_sz     = sz;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;

        return true;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Empty the array and free all allocated memory.
    ///////////////////////////////////////////////////////////////////////
    void release()
    {
        // Free array data buffer.
        SAFE_DELETE_ARRAY(_Base::m_p);
        
        // Clear sizes.
        _Base::m_sz     = 0;
        _Base::m_asz[0] = 0;
        _Base::m_asz[1] = 0;
    }
};

} // namespace

#endif
