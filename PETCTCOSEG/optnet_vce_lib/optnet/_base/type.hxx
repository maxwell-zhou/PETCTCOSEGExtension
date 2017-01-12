/*
 ==========================================================================
 |   
 |   $Id: type.hxx 80 2005-01-22 07:49:12Z Administrator $
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

#ifndef ___TYPE_HXX___
#   define ___TYPE_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <cstddef>
#   include <optnet/_base/iterator.hxx>

/// @namespace optnet
namespace optnet {

class dummy_type
{
public:
    typedef int                 value_type;
    typedef value_type*         pointer;
    typedef const value_type*   const_pointer;
#   ifndef __OPTNET_CRAPPY_MSC__
    typedef optnet::detail::normal_iterator<pointer, dummy_type>
                                iterator;
    typedef optnet::detail::normal_iterator<const_pointer, dummy_type>
                                const_iterator;
#   else
    typedef value_type*         iterator;
    typedef const value_type*   const_iterator;
#   endif
    typedef value_type&         reference;
    typedef const value_type&   const_reference;
    typedef ptrdiff_t           difference_type;
    typedef size_t              size_type;
};

} // namespace


#endif 
