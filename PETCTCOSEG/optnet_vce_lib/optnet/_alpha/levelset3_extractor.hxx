/*
 ==========================================================================
 |   
 |   $Id: levelset3_extractor.hxx 130 2005-02-06 04:27:30Z kangli $
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

#ifndef ___LEVELSET3_EXTRACTOR_HXX___
#   define ___LEVELSET3_EXTRACTOR_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <limits>
#   include <cmath>

#   ifdef max
#       undef max
#   endif

namespace optnet {

template <typename _Ty,
          typename _Node,
          typename _NodeContainer,
          typename _Tg = net_f_xy>
class levelset3_extractor
{
public:
    
    typedef _Ty                         value_type;
    typedef value_type&                 reference;
    typedef const value_type&           const_reference;

    typedef array_base<value_type, _Tg> array_base_type;
    typedef array_ref<value_type, _Tg>  array_ref_type;
    typedef array<value_type, _Tg>      array_type;

    typedef _Node                       node_type;
    typedef node_type&                  node_reference;
    typedef const node_type&            node_const_reference;

    typedef _NodeContainer              node_container;

    typedef size_t                      size_type;

    enum side_type
    {
        REMOTE =  0,
        INNER  =  1,
        OUTER  = -1
    };
    

    ///////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    ///////////////////////////////////////////////////////////////////////
    levelset3_extractor() :
        m_levelset_value(value_type())
    {
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Constructor.
    ///
    ///  @param levelset_value  The value of the level set front.
    ///////////////////////////////////////////////////////////////////////
    levelset3_extractor(const value_type& levelset_value) :
        m_levelset_value(levelset_value)
    {
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Extract the nodes adjacent to the the specified level set front.
    ///
    ///  @param levelset     The embedded level set funtion.
    ///  @param inner_nodes  The vector containing the nodes that are on
    ///                      or immediately inside the level set front.
    ///  @param outer_nodes  The vector containing the nodes that are
    ///                      immediately outside the level set front.
    ///
    ///////////////////////////////////////////////////////////////////////
    void extract(array_base_type& levelset,
                 node_container & inner_nodes,
                 node_container & outer_nodes
                 )
    {
        size_type  i0, i1, i2;
        value_type out;

        for (i2 = 0; i2 < levelset.size_2(); ++i2) {
            for (i1 = 0; i1 < levelset.size_1(); ++i1) {
                for (i0 = 0; i0 < levelset.size_0(); ++i0) {

                    side_type side = redistance(levelset, i0, i1, i2, out);
                    
                    if (side == INNER) {
                        node_type node;
                        node.index.v[0] = i0;
                        node.index.v[1] = i1;
                        node.index.v[2] = i2;
                        node.value      = out;
                        inner_nodes.push_back(node);
                    }
                    else if (side == OUTER) {
                        node_type node;
                        node.index.v[0] = i0;
                        node.index.v[1] = i1;
                        node.index.v[2] = i2;
                        node.value      = out;
                        outer_nodes.push_back(node);
                    }

                } // i0
            } // i1
        } // i2
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Sets the value of the level set front.
    ///
    ///  @param levelset_value The value of the level set front.
    ///////////////////////////////////////////////////////////////////////
    inline void set_levelset_value(const value_type& levelset_value)
    {
        m_levelset_value = levelset_value;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the value of the level set front.
    ///
    ///  @return The value of the level set front.
    ///////////////////////////////////////////////////////////////////////
    inline const_reference get_levelset_value() const
    {
        return m_levelset_value;
    }


private:

    value_type  m_levelset_value;

    ///////////////////////////////////////////////////////////////////////
    inline side_type redistance(array_base_type& levelset,
                                size_type        i0,
                                size_type        i1,
                                size_type        i2,
                                value_type&      out
                                )
    {
        static const double LARGE_DISTANCE
            = std::numeric_limits<double>::max();

        bool        inside;
        double      distance, neighbors[3];
        value_type  neighvalue, value = levelset(i0, i1, i2);

        if (value == m_levelset_value) {
            // nodes with value equal to the level set value are
            // classified as inner nodes.
            return INNER;
        }
        
        value  -= m_levelset_value;
        inside  = (value >= (value_type)(0));
        
        neighbors[0] = neighbors[1] = neighbors[2] = LARGE_DISTANCE;

        if (i0 > 0) {
            neighvalue = levelset(i0 - 1, i1, i2) - m_levelset_value;
            if ((neighvalue <  0 &&  inside) ||
                (neighvalue >= 0 && !inside)) {
                distance = value / (double)(value - neighvalue);
                neighbors[0] = distance;
            }
        }

        if (i0 + 1 < levelset.size_0()) {
            neighvalue = levelset(i0 + 1, i1, i2) - m_levelset_value;
            if ((neighvalue <  0 &&  inside) ||
                (neighvalue >= 0 && !inside)) {
                distance = value / (double)(value - neighvalue);
                if (neighbors[0] > distance) {
                    neighbors[0] = distance;
                }
            }
        }

        if (i1 > 0) {
            neighvalue = levelset(i0, i1 - 1, i2) - m_levelset_value;
            if ((neighvalue <  0 &&  inside) ||
                (neighvalue >= 0 && !inside)) {
                distance = value / (double)(value - neighvalue);
                neighbors[1] = distance;
            }
        }

        if (i1 + 1 < levelset.size_1()) {
            neighvalue = levelset(i0, i1 + 1, i2) - m_levelset_value;
            if ((neighvalue <  0 &&  inside) ||
                (neighvalue >= 0 && !inside)) {
                distance = value / (double)(value - neighvalue);
                if (neighbors[1] > distance) {
                    neighbors[1] = distance;
                }
            }
        }

        if (i2 > 0) {
            neighvalue = levelset(i0, i1, i2 - 1) - m_levelset_value;
            if ((neighvalue <  0 &&  inside) ||
                (neighvalue >= 0 && !inside)) {
                distance = value / (double)(value - neighvalue);
                neighbors[2] = distance;
            }
        }

        if (i2 + 1 < levelset.size_2()) {
            neighvalue = levelset(i0, i1, i2 + 1) - m_levelset_value;
            if ((neighvalue <  0 &&  inside) ||
                (neighvalue >= 0 && !inside)) {
                distance = value / (double)(value - neighvalue);
                if (neighbors[2] > distance) {
                    neighbors[2] = distance;
                }
            }
        }

        std::sort(neighbors, neighbors + 3);

        if (neighbors[0] >= LARGE_DISTANCE) return REMOTE;

        distance = 1.0 / (neighbors[0] * neighbors[0]);
        if (neighbors[1] < LARGE_DISTANCE)
            distance += 1.0 / (neighbors[1] * neighbors[1]);
        if (neighbors[2] < LARGE_DISTANCE)
            distance += 1.0 / (neighbors[2] * neighbors[2]);
        distance = sqrt(1.0 / distance);

        if (inside) {
            out = (value_type)(+distance);
            return INNER;
        }
        else {
            out = (value_type)(-distance);
            return OUTER;
        }
    }

};

} // namespace

#endif // ___LEVELSET3_EXTRACTOR_HXX___
