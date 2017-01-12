/*
 ==========================================================================
 |   
 |   $Id: graph_bk.hxx 21 2005-01-14 15:52:31Z kangli $
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
     
     This file implements a graph data structure using the forward-star
     representation.

 ==========================================================================
 */

#ifndef ___GRAPH_BK_HXX___
#   define ___GRAPH_BK_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/chunk_list.hxx>
#   include <optnet/_base/graph_traits.hxx>

/// @namespace optnet
namespace optnet {

    /// @namespace optnet::xtra
    /// @brief The namespace that contains third-party algorithms.
    namespace xtra {

///////////////////////////////////////////////////////////////////////////
///  @class graph_bk
///  @brief A graph class using the forward-star representation scheme.
///         This class is prepared for implementing the original
///         Boykov-Kolmogorov algorithm.
///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg = net_f_xy>
class graph_bk
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
        _Cap            cap;            // Residual capacity.
        _Cap            rev_cap;        // Reverse residual capacity.
    };

    // Reverse arc struct.
    struct rev_arc
    {
        fwd_arc*        p_fwd;          // Pointer to the forward arc.
    };


public:

    typedef optnet::array<node, _Tg>            node_container;
    typedef optnet::chunk_list<fwd_arc>         forward_arc_container;
    typedef optnet::chunk_list<rev_arc>         reverse_arc_container;

    typedef optnet::graph_traits<
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
    graph_bk();

    ///////////////////////////////////////////////////////////////////////
    ///  Construct a graph_bk object of the given size.
    ///
    ///  @param  s0   The size of the first  dimension of the graph. 
    ///  @param  s1   The size of the second dimension of the graph. 
    ///  @param  s2   The size of the third  dimension of the graph. 
    ///  @param  s3   The size of the fourth dimension of the graph
    ///               (default: 1).
    ///  @param  s4   The size of the fifth  dimension of the graph
    ///               (default: 1).
    ///
    ///////////////////////////////////////////////////////////////////////
    graph_bk(size_type s0, 
             size_type s1, 
             size_type s2, 
             size_type s3 = 1, 
             size_type s4 = 1
             );

    ///////////////////////////////////////////////////////////////////////
    /// Default destructor.
    ///////////////////////////////////////////////////////////////////////
    virtual ~graph_bk() {}

    ///////////////////////////////////////////////////////////////////////
    ///  Create a graph of the given size.
    ///
    ///  @param  s0   The size of the first  dimension of the graph. 
    ///  @param  s1   The size of the second dimension of the graph. 
    ///  @param  s2   The size of the third  dimension of the graph. 
    ///  @param  s3   The size of the fourth dimension of the graph
    ///               (default: 1).
    ///  @param  s4   The size of the fifth  dimension of the graph
    ///               (default: 1).
    ///
    ///  @returns Returns true if the graph is successfully created,
    ///           false otherwise.
    ///
    ///  @remarks The created graph can have up to five dimensions.
    ///
    ///////////////////////////////////////////////////////////////////////
    virtual bool create(size_type s0, 
                        size_type s1, 
                        size_type s2, 
                        size_type s3 = 1,
                        size_type s4 = 1
                        );

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the total number of nodes in the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size()   const { return m_nodes.size();   }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the first  dimension of the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size_0() const { return m_nodes.size_0(); }
    
    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the second dimension of the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size_1() const { return m_nodes.size_1(); }
    
    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the third  dimension of the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size_2() const { return m_nodes.size_2(); }
    
    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the fourth dimension of the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size_3() const { return m_nodes.size_3(); }
    
    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the fifth  dimension of the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size_4() const { return m_nodes.size_4(); }

    ///////////////////////////////////////////////////////////////////////
    ///  Add an arc connecting from the specified tail node to the
    ///  specified head node.
    ///
    ///  @param  tail_i0 The first  index of the tail node.
    ///  @param  tail_i1 The second index of the tail node.
    ///  @param  tail_i2 The third  index of the tail node.
    ///  @param  head_i0 The first  index of the head node.
    ///  @param  head_i1 The second index of the head node.
    ///  @param  head_i2 The third  index of the head node.
    ///  @param  cap     Capacity of the arc.
    ///  @param  rev_cap Reverse capacity of the arc.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_arc(size_type     tail_i0,
                        size_type     tail_i1,
                        size_type     tail_i2,
                        size_type     head_i0,
                        size_type     head_i1,
                        size_type     head_i2,
                        capacity_type cap,
                        capacity_type rev_cap
                        )
    {
        node_pointer p_tail_node = &(m_nodes(tail_i0,
                                             tail_i1,
                                             tail_i2
                                             ));
        node_pointer p_head_node = &(m_nodes(head_i0,
                                             head_i1,
                                             head_i2
                                             ));

        add_arc_helper(p_tail_node, p_head_node, cap, rev_cap);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Add an arc connecting from the specified tail node to the
    ///  specified head node.
    ///
    ///  @param  tail_i0 The first  index of the tail node.
    ///  @param  tail_i1 The second index of the tail node.
    ///  @param  tail_i2 The third  index of the tail node.
    ///  @param  tail_i3 The fourth index of the tail node.
    ///  @param  head_i0 The first  index of the head node.
    ///  @param  head_i1 The second index of the head node.
    ///  @param  head_i2 The third  index of the head node.
    ///  @param  head_i3 The fourth index of the head node.
    ///  @param  cap     Capacity of the arc.
    ///  @param  rev_cap Reverse capacity of the arc.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_arc(size_type tail_i0,
                        size_type tail_i1,
                        size_type tail_i2,
                        size_type tail_i3,
                        size_type head_i0,
                        size_type head_i1,
                        size_type head_i2,
                        size_type head_i3,
                        capacity_type cap,
                        capacity_type rev_cap
                        )
    {
        node_pointer p_tail_node = &(m_nodes(tail_i0,
                                             tail_i1,
                                             tail_i2,
                                             tail_i3
                                             ));
        node_pointer p_head_node = &(m_nodes(head_i0,
                                             head_i1,
                                             head_i2,
                                             head_i3
                                             ));

        add_arc_helper(p_tail_node, p_head_node, cap, rev_cap);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Add an arc connecting from the specified tail node to the
    ///  specified head node.
    ///
    ///  @param  tail_i0 The first  index of the tail node.
    ///  @param  tail_i1 The second index of the tail node.
    ///  @param  tail_i2 The third  index of the tail node.
    ///  @param  tail_i3 The fourth index of the tail node.
    ///  @param  tail_i4 The fifth  index of the tail node.
    ///  @param  head_i0 The first  index of the head node.
    ///  @param  head_i1 The second index of the head node.
    ///  @param  head_i2 The third  index of the head node.
    ///  @param  head_i3 The fourth index of the head node.
    ///  @param  head_i4 The fifth  index of the head node.
    ///  @param  cap     Capacity of the arc.
    ///  @param  rev_cap Reverse capacity of the arc.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_arc(size_type     tail_i0,
                        size_type     tail_i1,
                        size_type     tail_i2,
                        size_type     tail_i3,
                        size_type     tail_i4,
                        size_type     head_i0,
                        size_type     head_i1,
                        size_type     head_i2,
                        size_type     head_i3,
                        size_type     head_i4,
                        capacity_type cap,
                        capacity_type rev_cap
                        )
    {
        node_pointer p_tail_node = &(m_nodes(tail_i0,
                                             tail_i1,
                                             tail_i2,
                                             tail_i3,
                                             tail_i4
                                             ));
        node_pointer p_head_node = &(m_nodes(head_i0,
                                             head_i1,
                                             head_i2,
                                             head_i3,
                                             head_i4
                                             ));

        add_arc_helper(p_tail_node, p_head_node, cap, rev_cap);
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

    ///////////////////////////////////////////////////////////////////////
    ///  This is an helper function to the 'add_arc()' calls.
    ///////////////////////////////////////////////////////////////////////
    inline void add_arc_helper(node_pointer  p_tail_node,
                               node_pointer  p_head_node,
                               capacity_type cap,
                               capacity_type rev_cap
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

        // Set arc capacity and reverse capacity.
        p_fwd_arc->cap       = cap;
        p_fwd_arc->rev_cap   = rev_cap;
        
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


    static const unsigned char  SPECIAL_IN, SPECIAL_OUT;

    node_container              m_nodes;
    forward_arc_container       m_fwd_arcs;
    reverse_arc_container       m_rev_arcs;
    bool                        m_prepared;
};

    } // namespace
} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_xtra/graph_bk.cxx>
#   endif

#endif
