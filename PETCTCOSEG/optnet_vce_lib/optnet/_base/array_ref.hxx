/*
 ==========================================================================
 |   
 |   $Id: array_ref.hxx 58 2005-01-21 01:44:10Z kangli $
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

#ifndef ___ARRAY_REF_HXX___
#   define ___ARRAY_REF_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4244)
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/array_base.hxx>

///  @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class array_ref array_ref.hxx "_base/array_ref.hxx"
///  @brief Array template class.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg = net_f_xy>
class array_ref: public array_base<_Ty, _Tg>
{
    typedef array_base<_Ty, _Tg>    _Base;

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
    array_ref();
    explicit array_ref(_Ty* p, 
              size_type s0, 
              size_type s1, 
              size_type s2, 
              size_type s3 = 1, 
              size_type s4 = 1
              );

    // I have to put this function inline because of the buggy
    // Visual C++ 6.0.
    template <typename _Ta>
    array_ref(const array_base<_Ty, _Ta>& robj)
    {
        _Base::m_sz = robj.size();
        ::memcpy(_Base::m_asz, 
                 robj.sizes(), 
                 OPTNET_ARRAY_DIMS * sizeof(size_type));

        _Base::m_p = const_cast<_Ty*>(robj.data());

        // Initialize the look-up-table.
        _Base::init_lut();
    }
    
    array_ref(const array_ref& robj);
    ~array_ref();

    // modifiers
    array_ref& operator=(const array_ref& robj);

    // I have to put this function inline because of the buggy
    // Visual C++ 6.0.
    template <typename _Ta>
    array_ref&
    operator=(const array_base<_Ty, _Ta>& robj)
    {
        _Base::m_sz = robj.size();
        ::memcpy(_Base::m_asz, 
                 robj.sizes(), 
                 OPTNET_ARRAY_DIMS * sizeof(size_type));

        _Base::m_p = const_cast<_Ty*>(robj.data());

        // Reinitialize the look-up-table.
        _Base::free_lut();
        _Base::init_lut();

        return *this;
    }


    bool            assign (const _Ty* p,
                            size_type s0, 
                            size_type s1, 
                            size_type s2 = 1, 
                            size_type s3 = 1, 
                            size_type s4 = 1
                            );
    bool            reshape(size_type s0, 
                            size_type s1, 
                            size_type s2 = 1, 
                            size_type s3 = 1, 
                            size_type s4 = 1
                            );

};


} // namespace


#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_base/array_ref.cxx>
#   endif

#endif
