/*
 ==========================================================================
 |   
 |   $Id: array2_ref.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___ARRAY2_REF_HXX___
#   define ___ARRAY2_REF_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/_base/array2_base.hxx>

///  @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class array2_ref array2_ref.hxx "_base/array2_ref.hxx"
///  @brief Array template class.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty>
class array2_ref: public array2_base<_Ty>
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
    array2_ref() :
        array2_base<_Ty>()
    {
    }
    
    explicit array2_ref(_Ty* p, size_type s0, size_type s1)
    {
        assert(p != 0);
        assert(s0 > 0);
        assert(s1 > 0);

        _Base::m_p      = p;
        _Base::m_sz     = s0 * s1;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
    }

    array2_ref(const array2_base<_Ty>& robj)
    {
        assert(&robj != this);
        
        _Base::m_p = const_cast<_Ty*>(robj.data());

        _Base::m_asz[0] = robj.sizes()[0];
        _Base::m_asz[1] = robj.sizes()[1];
        _Base::m_sz     = robj.size();
    }
    
    array2_ref(const array2_ref& robj)
    {
        assert(&robj != this);

        _Base::m_p      = const_cast<_Ty*>(robj._Base::m_p);
        _Base::m_asz[0] = robj._Base::m_asz[0];
        _Base::m_asz[1] = robj._Base::m_asz[1];
        _Base::m_sz     = robj._Base::m_sz;
    }
    
    ~array2_ref()
    {
    }

    // modifiers
    array2_ref& operator=(const array2_ref& robj)
    {
        if (&robj != this) {
            _Base::m_p      = const_cast<_Ty*>(robj._Base::m_p);
            _Base::m_asz[0] = robj._Base::m_asz[0];
            _Base::m_asz[1] = robj._Base::m_asz[1];
            _Base::m_sz     = robj._Base::m_sz;
        }
        return *this;
    }

    array2_ref& operator=(const array2_base<_Ty>& robj)
    {
        _Base::m_p  = const_cast<_Ty*>(robj.data());

        _Base::m_asz[0] = robj.sizes()[0];
        _Base::m_asz[1] = robj.sizes()[1];
        _Base::m_sz     = robj.size();

        return *this;
    }

    bool assign(const _Ty* p, size_type s0, size_type s1)
    {
        assert(s0 > 0);
        assert(s1 > 0);

        if (0 == p) return false;

        _Base::m_p      = p;
        _Base::m_sz     = s0 * s1;
        _Base::m_asz[0] = s0;
        _Base::m_asz[1] = s1;
        
        return true;
    }    
};

} // namespace

#endif
