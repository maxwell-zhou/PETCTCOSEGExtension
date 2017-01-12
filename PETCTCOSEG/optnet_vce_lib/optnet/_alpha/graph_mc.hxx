/*
 ==========================================================================
 |   
 |   $Id: graph_mc.hxx 21 2005-01-14 15:52:31Z kangli $
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

/*
 ==========================================================================
  - Purpose:
     
     This file implements a multicolumn graph data structure using the
     forward-star representation.

 ==========================================================================
 */

#ifndef ___GRAPH_MC_HXX___
#   define ___GRAPH_MC_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/array2.hxx>
#   include <optnet/_base/array2_ref.hxx>
#   include <optnet/_base/chunk_list.hxx>
#   include <optnet/_base/graph_traits.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class graph_mc
///  @brief A multicolumn graph class with forward-star representation.
///////////////////////////////////////////////////////////////////////////
template <typename _Cap>
class graph_mc
{
    struct fwd_arc;
    struct rev_arc;

    struct node
    {
        fwd_arc*        p_first_out_arc; // Pointer to the first out arc.
        rev_arc*        p_first_in_arc;  // Pointer to the first in arc.
        fwd_arc*        p_parent_arc;    // Pointer to the parent arc.
        size_t          dist_id;
        size_t          dist;
        _Cap            cap;
        unsigned char   tag;
    };

    // Forward arc struct.
    struct fwd_arc
    {
        ptrdiff_t       shift;          // Tail node + shift = head node.
        _Cap            rev_cap;        // Reverse residual capacity.
    };

    // Reverse arc struct.
    struct rev_arc
    {
        fwd_arc*        p_fwd;          // Pointer to the forward arc.
    };


public:

    typedef array2<node>                        node_container;
    typedef chunk_list<fwd_arc>                 forward_arc_container;
    typedef chunk_list<rev_arc>                 reverse_arc_container;

    typedef graph_traits<
        _Cap, 
        node_container, 
        forward_arc_container, 
        reverse_arc_container>                  _Traits;

    typedef typename
        _Traits::node_type                      node_type;
    typedef typename 
        _Traits::node_reference                 node_reference;
    typedef typename 
        _Traits::node_const_reference           node_const_reference;
    typedef typename 
        _Traits::node_iterator                  node_iterator;
    typedef typename 
        _Traits::node_const_iterator            node_const_iterator;
    typedef typename 
        _Traits::node_pointer                   node_pointer;
    typedef typename 
        _Traits::node_const_pointer             node_const_pointer;

    typedef typename 
        _Traits::forward_arc_type               forward_arc_type;
    typedef typename 
        _Traits::forward_arc_reference          forward_arc_reference;
    typedef typename 
        _Traits::forward_arc_const_reference    forward_arc_const_reference;
    typedef typename 
        _Traits::forward_arc_iterator           forward_arc_iterator;
    typedef typename 
        _Traits::forward_arc_const_iterator     forward_arc_const_iterator;
    typedef typename 
        _Traits::forward_arc_pointer            forward_arc_pointer;
    typedef typename 
        _Traits::forward_arc_const_pointer      forward_arc_const_pointer;

    typedef typename 
        _Traits::reverse_arc_type               reverse_arc_type;
    typedef typename 
        _Traits::reverse_arc_reference          reverse_arc_reference;
    typedef typename 
        _Traits::reverse_arc_const_reference    reverse_arc_const_reference;
    typedef typename 
        _Traits::reverse_arc_iterator           reverse_arc_iterator;
    typedef typename 
        _Traits::reverse_arc_const_iterator     reverse_arc_const_iterator;
    typedef typename 
        _Traits::reverse_arc_pointer            reverse_arc_pointer;
    typedef typename 
        _Traits::reverse_arc_const_pointer      reverse_arc_const_pointer;

    typedef typename _Traits::capacity_type     capacity_type;
    typedef typename _Traits::difference_type   difference_type;
    typedef typename _Traits::size_type         size_type;

    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    graph_mc();

    ///////////////////////////////////////////////////////////////////////
    ///  Construct a graph_mc object of the given size.
    ///
    ///  @param  colsize  The number of nodes in a column.
    ///  @param  numcols  The number of columns.
    ///
    ///////////////////////////////////////////////////////////////////////
    graph_mc(size_type colsize, size_type numcols);

    ///////////////////////////////////////////////////////////////////////
    /// Default destructor.
    ///////////////////////////////////////////////////////////////////////
    virtual ~graph_mc() {}

    ///////////////////////////////////////////////////////////////////////
    ///  Create a graph of the given size.
    ///
    ///  @param  colsize  The number of nodes in a column.
    ///  @param  numcols  The number of columns.
    ///
    ///  @returns Returns true if the graph is successfully created,
    ///           false otherwise.
    ///
    ///////////////////////////////////////////////////////////////////////
    virtual bool create(size_type colsize, size_type numcols);

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the number of columns in the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type num_columns() const { return m_nodes.size_1(); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the number of nodes in a column.
    ///////////////////////////////////////////////////////////////////////
    inline size_type nodes_per_column() const { return m_nodes.size_0(); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the total number of nodes in the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size() const { return m_nodes.size(); }

    ///////////////////////////////////////////////////////////////////////
    ///  Add an arc connecting from the specified tail node to the
    ///  specified head node.
    ///
    ///  @param  tail_hi  The height-index of the tail node.
    ///  @param  tail_ci  The column-index of the tail node.
    ///  @param  head_hi  The height-index of the head node.
    ///  @param  head_ci  The column-index of the head node.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_arc(size_type tail_hi,
                        size_type tail_ci,
                        size_type head_hi,
                        size_type head_ci
                        )
    {
        node_pointer p_tail_node = &(m_nodes(tail_hi, tail_ci));
        node_pointer p_head_node = &(m_nodes(head_hi, head_ci));

        add_arc(p_tail_node, p_head_node);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Add an arc connecting from the specified tail node to the
    ///  specified head node.
    ///
    ///  @param  p_tail_node  A pointer to the tail node.
    ///  @param  p_head_node  A pointer to the head node.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_arc(node_pointer p_tail_node,
                        node_pointer p_head_node
                        )
    {
        assert(m_prepared == false);

        // Allocate new forward and reverse arcs and append them to the
        // corresponding arc lists.
        rev_arc* p_rev_arc   = m_rev_arcs.append();
        fwd_arc* p_fwd_arc   = m_fwd_arcs.append();

        // Temporarily store the tail and head nodes in the p_fwd and
        // shift fields of the reverse and forward arcs.
        p_rev_arc->p_fwd     = reinterpret_cast<fwd_arc *>(p_tail_node);
        p_fwd_arc->shift     = reinterpret_cast<ptrdiff_t>(p_head_node);

        // The reverse residual capacity is initialized to be zero.
        p_fwd_arc->rev_cap   = static_cast<capacity_type>(0);
        
        // Temporarily store the number of outgoing and incoming arcs
        // in the p_first_out_arc and p_first_in_arc fields.
        p_tail_node->p_first_out_arc = reinterpret_cast<fwd_arc*>(
            reinterpret_cast<ptrdiff_t>
                (p_tail_node->p_first_out_arc) + 1
            );
        p_head_node->p_first_in_arc = reinterpret_cast<rev_arc*>(
            reinterpret_cast<ptrdiff_t>
                (p_head_node->p_first_in_arc) + 1
            );
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Remove all arcs that were added.
    ///////////////////////////////////////////////////////////////////////
    inline void clear_arcs()
    {
        if (m_fwd_arcs.count() > 0 || m_rev_arcs.count() > 0) {

            m_fwd_arcs.clear();
            m_rev_arcs.clear();

            for (node_iterator it = m_nodes.begin();
                 it != m_nodes.end();
                 ++it) {

                it->p_first_out_arc
                    = reinterpret_cast<fwd_arc*>(0);
                it->p_first_in_arc
                    = reinterpret_cast<rev_arc*>(0);
            }
        }
        
        // Set "prepared" flag to false.
        // Call prepare() to make the graph usable.
        m_prepared = false;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Converts the arcs added by 'add_arc()' calls to a forward-star
    ///  representaion, i.e., the arcs are sorted in the order of their
    ///  tail nodes.
    ///
    ///  @remarks After calling this function, the graph is ready for
    ///           use, but at the mean time, one cannot call the
    ///           function 'add_arc()' any more.
    ///////////////////////////////////////////////////////////////////////
    void prepare();


protected:

    // Constants (Initialized in graph_mc.cxx)
    static const unsigned char  SPECIAL_IN;
    static const unsigned char  SPECIAL_OUT;

    node_container              m_nodes;
    forward_arc_container       m_fwd_arcs;
    reverse_arc_container       m_rev_arcs;
    bool                        m_prepared;
};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_alpha/graph_mc.cxx>
#   endif

#endif
