/*
 ==========================================================================
 |   
 |   $Id: levelset2_extractor.hxx 44 2005-01-19 03:18:25Z kangli $
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

#ifndef ___LEVELSET2_EXTRACTOR_HXX___
#   define ___LEVELSET2_EXTRACTOR_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/array2.hxx>
#   include <optnet/_base/array2_ref.hxx>
#   include <cmath>

namespace optnet {

template <typename _Ty,
          typename _Node,
          typename _NodeContainer>
class levelset2_extractor
{
public:
    
    typedef _Ty                         value_type;
    typedef value_type&                 reference;
    typedef const value_type&           const_reference;

    typedef array2_base<value_type>     array_base_type;
    typedef array2_ref<value_type>      array_ref_type;
    typedef array2<value_type>          array_type;

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
    levelset2_extractor() :
        m_levelset_value(value_type())
    {
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Constructor.
    ///
    ///  @param levelset_value  The value of the level set front.
    ///////////////////////////////////////////////////////////////////////
    levelset2_extractor(const value_type& levelset_value) :
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
        size_type  i0, i1;
        value_type out;

        for (i1 = 0; i1 < levelset.size_1(); ++i1) {
            for (i0 = 0; i0 < levelset.size_0(); ++i0) {

                side_type side = redistance(levelset, i0, i1, out);
                
                if (side == INNER) {
                    node_type node;
                    node.index.v[0] = i0;
                    node.index.v[1] = i1;
                    node.value      = out;
                    inner_nodes.push_back(node);
                }
                else if (side == OUTER) {
                    node_type node;
                    node.index.v[0] = i0;
                    node.index.v[1] = i1;
                    node.value      = out;
                    outer_nodes.push_back(node);
                }

            } // i0
        } // i1
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
                                value_type&      out
                                )
    {
        bool        inside, assigned[2];
        double      distance, neighbors[2];
        value_type  neighvalue, value = levelset(i0, i1);

        if (value == m_levelset_value) {
            // nodes with value equal to the level set value are
            // classified as inner nodes.
            return INNER;
        }
        
        value  -= m_levelset_value;
        inside  = (value >= (value_type)(0));
        
        // assigned will be true if the neighborhood is crossing the
        // levelset contour
        assigned[0] = assigned[1] = false;

        if (i0 > 0) {
            neighvalue = levelset(i0 - 1, i1) - m_levelset_value;
            if ((neighvalue <  0 &&  inside) ||
                (neighvalue >= 0 && !inside)) {
                distance = value / (double)(value - neighvalue);
                neighbors[0] = distance;
                assigned[0] = true;
            }
        }

        if (i0 + 1 < levelset.size_0()) {
            neighvalue = levelset(i0 + 1, i1) - m_levelset_value;
            if ((neighvalue <  0 &&  inside) ||
                (neighvalue >= 0 && !inside)) {
                distance = value / (double)(value - neighvalue);
                if (!assigned[0] || (neighbors[0] > distance)) {
                    neighbors[0] = distance;
                    assigned[0] = true;
                }
            }
        }

        if (i1 > 0) {
            neighvalue = levelset(i0, i1 - 1) - m_levelset_value;
            if ((neighvalue <  0 &&  inside) ||
                (neighvalue >= 0 && !inside)) {
                distance = value / (double)(value - neighvalue);
                neighbors[1] = distance;
                assigned[1] = true;
            }
        }

        if (i1 + 1 < levelset.size_1()) {
            neighvalue = levelset(i0, i1 + 1) - m_levelset_value;
            if ((neighvalue <  0 &&  inside) ||
                (neighvalue >= 0 && !inside)) {
                distance = value / (double)(value - neighvalue);
                if (!assigned[1] || (neighbors[1] > distance)) {
                    neighbors[1] = distance;
                    assigned[1] = true;
                }
            }
        }

        distance = 0;

        if (assigned[0] && assigned[1]) {
            double a = neighbors[0] * neighbors[0];
            double b = neighbors[1] * neighbors[1];
            distance = sqrt(a * b / (a + b));
        }
        else if (assigned[0]) {
            distance = neighbors[0];
        }
        else if (assigned[1]) {
            distance = neighbors[1];
        }
        else {
            return REMOTE;
        }

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

#endif // ___LEVELSET2_EXTRACTOR_HXX___
