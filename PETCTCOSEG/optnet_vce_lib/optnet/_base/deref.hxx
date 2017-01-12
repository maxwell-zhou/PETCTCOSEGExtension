/*
==========================================================================
|   
|   $Id: deref.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___DEREF_HXX___
#   define ___DEREF_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class deref
///  @brief Dereference an iterator/pointer.
///////////////////////////////////////////////////////////////////////////
template <typename _It>
struct deref
{
    inline typename _It::reference
        operator()(_It it) { return *it;  }
    inline typename _It::const_reference
        operator()(_It it) const { return *it;  }
};

///  Create a deref object.
template <typename _It>
inline deref<_It> make_deref (_It it)   { return deref<_It>(); }


///////////////////////////////////////////////////////////////////////////
///  @class deref_first
///  @brief Returns the first  member of a tuple.
///////////////////////////////////////////////////////////////////////////
template <typename _TupleIt>
struct deref_first
{
    inline typename _TupleIt::value_type::first_type &
        operator()(_TupleIt it) { return it->first;  }
    inline const typename _TupleIt::value_type::first_type &
        operator()(_TupleIt it) const { return it->first;  }
};

///  Create a deref_first  object.
template <typename _TupleIt>
inline deref_first <_TupleIt> make_deref_first (_TupleIt it)
{ return deref_first <_TupleIt>(); }


///////////////////////////////////////////////////////////////////////////
///  @class deref_second
///  @brief Returns the second member of a tuple.
///////////////////////////////////////////////////////////////////////////
template <typename _TupleIt>
struct deref_second
{
    inline typename _TupleIt::value_type::second_type&
        operator()(_TupleIt it) { return it->second; }
    inline const typename _TupleIt::value_type::second_type&
        operator()(_TupleIt it) const { return it->second; }
};

///  Create a deref_second object.
template <typename _TupleIt>
inline deref_second<_TupleIt> make_deref_second(_TupleIt it)
{ return deref_second<_TupleIt>(); }


} // namespace

#endif // ___DEREF_HXX___
