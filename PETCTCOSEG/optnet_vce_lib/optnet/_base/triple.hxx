/* 
 ==========================================================================
 |   
 |   $Id: triple.hxx 21 2005-01-14 15:52:31Z kangli $
 |
 |   Written by Kang Li <kangl@cmu.edu>
 |   College of Engineering Imaging Group
 |   The University of Iowa
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
 
#ifndef ___TRIPLE_HXX___
#   define ___TRIPLE_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/config.h>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class triple
///  @brief A class that provides for the ability to treat three objects
///         as a single object.
///////////////////////////////////////////////////////////////////////////
template <typename _T1, typename _T2, typename _T3>
struct triple {
    
    typedef _T1 first_type;     /// Type of the first element of triple.
    typedef _T2 second_type;    /// Type of the second element of triple.
    typedef _T3 third_type;     /// Type of the third element of triple.

    _T1 first;                  /// The first element of triple.
    _T2 second;                 /// The second element of triple.
    _T3 third;                  /// The third element of triple.
  
    ///////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    ///  @remarks This constructor initializes first element of the triple
    ///           to the default of first_type, second element to default
    ///           of second_type and third to default of third_type. 
    ///////////////////////////////////////////////////////////////////////
    triple() : first(_T1()), second(_T2()), third(_T3()) {}

    ///////////////////////////////////////////////////////////////////////
    ///  Constructor.
    ///  @param a Value initializing the first element of triple. 
    ///  @param b Value initializing the second element of triple. 
    ///  @param c Value initializing the third element of triple. 
    ///////////////////////////////////////////////////////////////////////
    triple(const _T1& a, const _T2& b, const _T3& c) : 
        first(a), second(b), third(c) {}

#ifdef __OPTNET_MEMBER_TEMPLATES__
    ///////////////////////////////////////////////////////////////////////
    ///  Copy constructor.
    ///  @param r A triple whose values are to be used to initialize the
    ///           elements of another triple.
    ///////////////////////////////////////////////////////////////////////
    template <typename _U1, typename _U2, typename _U3>
    triple(const triple<_U1, _U2, _U3>& r) : 
        first(r.first), second(r.second), third(r.third) {}
#endif

};

template <typename _T1, typename _T2, typename _T3>
inline bool operator== (const triple<_T1, _T2, _T3>& __x,
                        const triple<_T1, _T2, _T3>& __y
                        )
{ 
    return __x.first  == __y.first  &&
           __x.second == __y.second &&
           __x.third  == __y.third; 
}

template <typename _T1, typename _T2, typename _T3>
inline bool operator<  (const triple<_T1, _T2, _T3>& __x,
                        const triple<_T1, _T2, _T3>& __y
                        )
{ 
    return (__x.first < __y.first)
        || (!(__y.first < __x.first) && (__x.second < __y.second))
        || (
               !(__y.first  < __x.first)
            && !(__y.second < __x.second)
            &&  (__x.third  < __y.third)
            );
}


///////////////////////////////////////////////////////////////////////////
///  A template helper function used to construct triple objects,
///  where the component types are based on the datatypes passed
///  as parameters.
///
///  @param a Value initializing the first element of triple. 
///  @param b Value initializing the second element of triple. 
///  @param c Value initializing the third element of triple. 
///  @return  The constructed triple object: triple<_T1, _T2, _T3>(a, b, c).
///
///////////////////////////////////////////////////////////////////////////
template <typename _T1, typename _T2, typename _T3>
inline triple<_T1, _T2, _T3> make_triple(const _T1& a,
                                         const _T2& b,
                                         const _T3& c
                                         )
{
    return triple<_T1, _T2, _T3>(a, b, c);
}

} // namespace

#endif // ___TRIPLE_HXX___
