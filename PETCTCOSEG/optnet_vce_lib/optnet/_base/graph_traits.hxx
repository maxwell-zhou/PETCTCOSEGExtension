/*
 ==========================================================================
 |   
 |   $Id: graph_traits.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___GRAPH_TRAITS_HXX___
#   define ___GRAPH_TRAITS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/_base/type.hxx>


/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class graph_traits graph_traits.hxx "_base/graph_traits.hxx"
///  @brief Graph character traits class template.
///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _NodeCnt,
          typename _FwdArcCnt = dummy_type,
          typename _RevArcCnt = dummy_type>
class graph_traits
{
public:

    typedef _Cap                    capacity_type;

    typedef typename 
        _NodeCnt::value_type        node_type;
    typedef typename 
        _NodeCnt::reference         node_reference;
    typedef typename 
        _NodeCnt::const_reference   node_const_reference;
    typedef typename 
        _NodeCnt::iterator          node_iterator;
    typedef typename 
        _NodeCnt::const_iterator    node_const_iterator;
    typedef typename 
        _NodeCnt::pointer           node_pointer;
    typedef typename 
        _NodeCnt::const_pointer     node_const_pointer;

    typedef typename 
        _FwdArcCnt::value_type      forward_arc_type;
    typedef typename 
        _FwdArcCnt::reference       forward_arc_reference;
    typedef typename 
        _FwdArcCnt::const_reference forward_arc_const_reference;
    typedef typename 
        _FwdArcCnt::iterator        forward_arc_iterator;
    typedef typename 
        _FwdArcCnt::const_iterator  forward_arc_const_iterator;
    typedef typename 
        _FwdArcCnt::pointer         forward_arc_pointer;
    typedef typename 
        _FwdArcCnt::const_pointer   forward_arc_const_pointer;

    typedef typename
        _RevArcCnt::value_type      reverse_arc_type;
    typedef typename
        _RevArcCnt::reference       reverse_arc_reference;
    typedef typename
        _RevArcCnt::const_reference reverse_arc_const_reference;
    typedef typename 
        _RevArcCnt::iterator        reverse_arc_iterator;
    typedef typename 
        _RevArcCnt::const_iterator  reverse_arc_const_iterator;
    typedef typename 
        _RevArcCnt::pointer         reverse_arc_pointer;
    typedef typename 
        _RevArcCnt::const_pointer   reverse_arc_const_pointer;

    typedef ptrdiff_t               difference_type;
    typedef size_t                  size_type;
};

} // namespace


#endif 
