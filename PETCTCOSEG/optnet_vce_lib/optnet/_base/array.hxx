/*
 ==========================================================================
 |   
 |   $Id: array.hxx 146 2005-02-08 05:26:28Z kangli $
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

#ifndef ___ARRAY_HXX___
#   define ___ARRAY_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/_base/array_base.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma warning(disable: 4244)
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class array array.hxx "_base/array.hxx"
///  @brief Array template class.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg = net_f_xy>
class array: public array_base<_Ty, _Tg>
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
    array();
    array(const size_type* asz);
    array(size_type s0, 
          size_type s1, 
          size_type s2, 
          size_type s3 = 1, 
          size_type s4 = 1
          );
    array(const array& robj);
    array(const array_base<_Ty, _Tg>& robj);
    array(const array2_base<_Ty>& robj);
    ~array();

    template <typename _Ta>
    void copy_to(array<_Ta>& copy)
    {
        copy.create(_Base::m_asz[0], _Base::m_asz[1],
                    _Base::m_asz[2], _Base::m_asz[3],
                    _Base::m_asz[4]);

        const_iterator it_from = this->begin();
        typename array<_Ta>::iterator it_to = copy.begin();

        for (; it_from != this->end(); ++it_from, ++it_to) {
            *it_to = static_cast<_Ta>(*it_from);
        }
    }

    // modifiers
    array& operator=(const array& robj);
    array& operator=(const array_base<_Ty, _Tg>& robj);
    array& operator=(const array2_base<_Ty>& robj);

    bool            create (const size_type* asz);
    bool            create (size_type s0, 
                            size_type s1, 
                            size_type s2, 
                            size_type s3 = 1, 
                            size_type s4 = 1
                            );
    bool            create_and_fill(const value_type& value,
                                    size_type         s0, 
                                    size_type         s1, 
                                    size_type         s2, 
                                    size_type         s3 = 1, 
                                    size_type         s4 = 1
                                    );
    void            release();
    bool            reshape(size_type s0, 
                            size_type s1, 
                            size_type s2,
                            size_type s3 = 1, 
                            size_type s4 = 1
                            );
};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_base/array.cxx>
#   endif

#endif 
