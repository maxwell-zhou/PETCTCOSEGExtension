/*
 ==========================================================================
 |   
 |   $Id: indexer.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___INDEXER_HXX___
#   define ___INDEXER_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/tags.hxx>
#   include <cassert>

//
// Note: Microsoft VC6 does not support class template partial
//       specialization. For compatibility, we are passing all
//       parameters as function arguments. This makes the im-
//       plementation much more clumsy.
//

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class indexer
///  @brief Array indexer template.
///
///  This template should never be instantiated in real programs.
///  Use its specializations instead.
///////////////////////////////////////////////////////////////////////////
template <typename _Tg>
class indexer
{
public:
    template <typename _Ty> 
    inline _Ty& size_0(const _Ty* s) const 
    { 
        assert(false); 
        return _Ty(); 
    }
    template <typename _Ty> 
    inline _Ty& size_1(const _Ty* s) const
    {
        assert(false);
        return _Ty();
    }
    template <typename _Ty> 
    inline _Ty& size_2(const _Ty* s) const
    {
        assert(false);
        return _Ty();
    }
    template <typename _Ty> 
    inline _Ty& size_3(const _Ty* s) const
    { 
        assert(false);
        return _Ty();
    }
    template <typename _Ty> 
    inline _Ty& size_4(const _Ty* s) const
    {
        assert(false);
        return _Ty();
    }
    
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti* s
                     ) const
    {
        assert(false);
        return _Ty();
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti* s
                     ) const
    {
        assert(false);
        return _Ty();
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti* s
                     ) const
    {
        assert(false);
        return _Ty();
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_xy>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = f(i0, i1).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_xy>
{
public:
    template <typename _Ty> 
        inline _Ty& size_0(_Ty* s) const { return s[0]; }
    template <typename _Ty> 
        inline _Ty& size_1(_Ty* s) const { return s[1]; }
    template <typename _Ty> 
        inline _Ty& size_2(_Ty* s) const { return s[2]; }
    template <typename _Ty> 
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty> 
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti*
                     ) const
    {
       return p[i2][i1][i0];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti*
                     ) const
    {
       return p[i3][i2][i1][i0];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti*
                     ) const
    {
       return p[i4][i3][i2][i1][i0];
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_yz>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = f(i1, i2).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_yz>
{ 
public:
    template <typename _Ty>
        inline _Ty& size_0(_Ty* s) const { return s[1]; }
    template <typename _Ty>
        inline _Ty& size_1(_Ty* s) const { return s[2]; }
    template <typename _Ty>
        inline _Ty& size_2(_Ty* s) const { return s[0]; }
    template <typename _Ty>
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty>
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti*
                     ) const
    {
       return p[i1][i0][i2];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti*
                     ) const
    {
       return p[i3][i1][i0][i2];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti*
                     ) const
    {
       return p[i4][i3][i1][i0][i2];
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_zx>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = f(z, x).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_zx>
{ 
public:
    template <typename _Ty> 
        inline _Ty& size_0(_Ty* s) const { return s[2]; }
    template <typename _Ty> 
        inline _Ty& size_1(_Ty* s) const { return s[0]; }
    template <typename _Ty> 
        inline _Ty& size_2(_Ty* s) const { return s[1]; }
    template <typename _Ty> 
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty> 
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti*
                     ) const
    {
       return p[i0][i2][i1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti*
                     ) const
    {
       return p[i3][i0][i2][i1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti*
                     ) const
    {
       return p[i4][i3][i0][i2][i1];
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_xy_zflipped>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = s2 - 1 - f(i0, i1).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_xy_zflipped>
{
public:
    template <typename _Ty> 
        inline _Ty& size_0(_Ty* s) const { return s[0]; }
    template <typename _Ty> 
        inline _Ty& size_1(_Ty* s) const { return s[1]; }
    template <typename _Ty> 
        inline _Ty& size_2(_Ty* s) const { return s[2]; }
    template <typename _Ty> 
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty> 
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti* s
                     ) const
    {
       return p[s[2] - i2 - 1][i1][i0];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti* s
                     ) const
    {
       return p[i3][s[2] - i2 - 1][i1][i0];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti* s
                     ) const
    {
       return p[i4][i3][s[2] - i2 - 1][i1][i0];
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_yz_xflipped>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = s0 - 1 - f(i1, i2).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_yz_xflipped>
{ 
public:
    template <typename _Ty>
        inline _Ty& size_0(_Ty* s) const { return s[1]; }
    template <typename _Ty>
        inline _Ty& size_1(_Ty* s) const { return s[2]; }
    template <typename _Ty>
        inline _Ty& size_2(_Ty* s) const { return s[0]; }
    template <typename _Ty>
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty>
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti* s
                     ) const
    {
       return p[i1][i0][s[0] - i2 - 1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti* s
                     ) const
    {
       return p[i3][i1][i0][s[0] - i2 - 1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti* s
                     ) const
    {
       return p[i4][i3][i1][i0][s[0] - i2 - 1];
    }
};

///////////////////////////////////////////////////////////////////////////
///  @class indexer<net_f_zx_yflipped>
///  @brief Template specialization of indexer.
///
///  Array indexer for the net in the form Net = s1 - 1 - f(z, x).
///////////////////////////////////////////////////////////////////////////
template <>
class indexer<net_f_zx_yflipped>
{ 
public:
    template <typename _Ty> 
        inline _Ty& size_0(_Ty* s) const { return s[2]; }
    template <typename _Ty> 
        inline _Ty& size_1(_Ty* s) const { return s[0]; }
    template <typename _Ty> 
        inline _Ty& size_2(_Ty* s) const { return s[1]; }
    template <typename _Ty> 
        inline _Ty& size_3(_Ty* s) const { return s[3]; }
    template <typename _Ty> 
        inline _Ty& size_4(_Ty* s) const { return s[4]; }

    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***   p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     const _Ti* s
                     ) const
    {
       return p[i0][s[1] - i2 - 1][i1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty****  p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     const _Ti* s
                     ) const
    {
       return p[i3][i0][s[1] - i2 - 1][i1];
    }
    template <typename _Ty, typename _Ti>
    inline _Ty& data(_Ty***** p,
                     _Ti i0,
                     _Ti i1,
                     _Ti i2,
                     _Ti i3,
                     _Ti i4,
                     const _Ti* s
                     ) const
    {
       return p[i4][i3][i0][s[1] - i2 - 1][i1];
    }
};

} // namespace


#endif 
