/*
 ==========================================================================
 |   
 |   $Id: array_traits.hxx 80 2005-01-22 07:49:12Z Administrator $
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

#ifndef ___ARRAY_TRAITS_HXX___
#   define ___ARRAY_TRAITS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/type.hxx>
#   include <optnet/_base/iterator.hxx>

///  @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class array_traits array_traits.hxx "_base/array_traits.hxx"
///  @brief Array character traits class template.
///////////////////////////////////////////////////////////////////////////
template <typename _T, typename _Container>
class array_traits
{
public:

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that represents the data type used in an array.
    ///////////////////////////////////////////////////////////////////////
    typedef _T                  value_type;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a pointer to an element in an array.
    ///////////////////////////////////////////////////////////////////////
    typedef value_type*         pointer;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a pointer to a const element in an
    ///        array.
    ///////////////////////////////////////////////////////////////////////
    typedef const value_type*   const_pointer;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a random-access iterator that can read
    ///        or modify any element in a vector.
    ///////////////////////////////////////////////////////////////////////
#   ifndef __OPTNET_CRAPPY_MSC__
    typedef optnet::detail::normal_iterator<pointer, _Container>
                                iterator;
#   else
    typedef value_type*         iterator;
#   endif

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a random-access iterator that can read
    ///        a const element in an array.
    ///////////////////////////////////////////////////////////////////////
#   ifndef __OPTNET_CRAPPY_MSC__
    typedef optnet::detail::normal_iterator<const_pointer, _Container>
                                const_iterator;
#   else
    typedef const value_type*   const_iterator;
#   endif

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a reference to an element stored in an
    ///        array.
    ///////////////////////////////////////////////////////////////////////
    typedef value_type&         reference;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides a reference to a const element stored
    ///        in an array for reading and performing const operations.
    ///////////////////////////////////////////////////////////////////////
    typedef const value_type&   const_reference;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that provides the difference between the addresses of
    ///        two elements in an array.
    ///////////////////////////////////////////////////////////////////////
    typedef ptrdiff_t           difference_type;

    ///////////////////////////////////////////////////////////////////////
    /// @brief A type that counts the number of elements in an array.
    ///////////////////////////////////////////////////////////////////////
    typedef size_t              size_type;

};

} // namespace


#endif 
