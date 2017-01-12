/*
 ==========================================================================
 |   
 |   $Id: levelset_node.hxx 36 2005-01-18 00:11:45Z kangli $
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

#ifndef ___LEVELSET_NODE_HXX___
#   define ___LEVELSET_NODE_HXX___

namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class levelset_node
///  @brief Represents a node in a descretized level set function.
///////////////////////////////////////////////////////////////////////////
template <typename _Time, typename _Point>
struct levelset_node
{
    levelset_node() : value(_Time()) {}

    levelset_node(const levelset_node& rhs) :
        value(rhs.value), index(rhs.index) {}

    levelset_node& operator=(const levelset_node& rhs)
    {
        value = rhs.value;
        index = rhs.index;
        return *this;
    }

    inline bool
        index_equals(const levelset_node& rhs) const
            { return index == rhs.index; }

    inline bool operator==(const levelset_node& rhs) const
                { return value == rhs.value; }
    inline bool operator> (const levelset_node& rhs) const
                { return value >  rhs.value; }
    inline bool operator>=(const levelset_node& rhs) const
                { return value >= rhs.value; }
    inline bool operator< (const levelset_node& rhs) const
                { return value <  rhs.value; }
    inline bool operator<=(const levelset_node& rhs) const
                { return value <= rhs.value; }

    _Time   value;
    _Point  index;
};

} // namespace

#endif // ___LEVELSET_NODE_HXX___
