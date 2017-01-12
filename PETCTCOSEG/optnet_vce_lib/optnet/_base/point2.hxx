/*
 ==========================================================================
 |   
 |   $Id: point2.hxx 186 2005-02-27 02:04:22Z kangli $
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

#ifndef ___POINT2_HXX___
#   define ___POINT2_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/point.hxx>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class point
///  @brief 2-D point type.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
class point2 : public point<_T, 2>
{
    typedef point<_T, 2> _Base;

public:
    typedef typename _Base::value_type       value_type;
    typedef typename _Base::reference        reference;
    typedef typename _Base::const_reference  const_reference;
    typedef typename _Base::pointer          pointer;
    typedef typename _Base::const_pointer    const_pointer;
    typedef typename _Base::size_type        size_type;

    point2() :
        _Base()
    {
    }

    point2(const _Base& rhs) :
        _Base(rhs)
    {
    }

    point2(const value_type& x, const value_type& y)
    {
        _Base::v[0] = x;
        _Base::v[1] = y;
    }

    inline point2& operator=(const _Base& rhs)
    {
        _Base::operator=(rhs);
        return *this;
    }

    inline point2 offset(const value_type& ox,
                         const value_type& oy
                         ) const
    {
        return point2(_Base::v[0] + ox,
                      _Base::v[1] + oy
                      );
    }

};

} // namespace

#endif // ___POINT2_HXX___
