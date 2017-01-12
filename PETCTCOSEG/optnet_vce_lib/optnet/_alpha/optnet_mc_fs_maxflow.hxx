/*
 ==========================================================================
 |   
 |   $Id: optnet_mc_fs_maxflow.hxx 21 2005-01-14 15:52:31Z kangli $
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

      This file implements Boykov--Kolmogorov's max-flow/min-cut algorithm
      on a forward-star multicolumn graph.
      
      The implementation makes one or more simplifications, which are only
      applicable to the Optimal-Net algorithm, as listed below:

         * The residual capacity of non-st arcs is assumed to be +infinity
           and will not be decremented.


  - Reference(s):

    [1] Yuri Boykov and Vladimir Kolmogorov
        An Experimental Comparison of Min-Cut/Max-Flow Algorithms for 
            Energy Minimization in Vision
        IEEE Trans. on Pattern Analysis and Machine Intelligence, 2004
        URL: http://www.csd.uwo.ca/faculty/yuri/Abstracts/pami04-abs.html
 ==========================================================================
 */

#ifndef ___OPTNET_MC_FS_MAXFLOW_HXX___
#   define ___OPTNET_MC_FS_MAXFLOW_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4127)
#   endif

#   include <optnet/_alpha/graph_mc.hxx>
#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <queue>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class optnet_mc_fs_maxflow
///  @brief Implementation of the Boykov-Kolmogorov max-flow algorithm on
///         a forward-star represented graph.
///////////////////////////////////////////////////////////////////////////
template <typename _Cap>
class optnet_mc_fs_maxflow
    : public graph_mc<_Cap>
{
    typedef graph_mc<_Cap> _Base;

public:

    typedef typename
        _Base::node_container                   node_container;
    typedef typename 
        _Base::forward_arc_container            forward_arc_container;
    typedef typename 
        _Base::reverse_arc_container            reverse_arc_container;

    typedef typename
        _Base::node_type                        node_type;
    typedef typename 
        _Base::node_reference                   node_reference;
    typedef typename 
        _Base::node_const_reference             node_const_reference;
    typedef typename 
        _Base::node_iterator                    node_iterator;
    typedef typename 
        _Base::node_const_iterator              node_const_iterator;
    typedef typename 
        _Base::node_pointer                     node_pointer;
    typedef typename 
        _Base::node_const_pointer               node_const_pointer;

    typedef typename 
        _Base::forward_arc_type                 forward_arc_type;
    typedef typename 
        _Base::forward_arc_reference            forward_arc_reference;
    typedef typename 
        _Base::forward_arc_const_reference      forward_arc_const_reference;
    typedef typename 
        _Base::forward_arc_iterator             forward_arc_iterator;
    typedef typename 
        _Base::forward_arc_const_iterator       forward_arc_const_iterator;
    typedef typename 
        _Base::forward_arc_pointer              forward_arc_pointer;
    typedef typename 
        _Base::forward_arc_const_pointer        forward_arc_const_pointer;

    typedef typename 
        _Base::reverse_arc_type                 reverse_arc_type;
    typedef typename 
        _Base::reverse_arc_reference            reverse_arc_reference;
    typedef typename 
        _Base::reverse_arc_const_reference      reverse_arc_const_reference;
    typedef typename 
        _Base::reverse_arc_iterator             reverse_arc_iterator;
    typedef typename 
        _Base::reverse_arc_const_iterator       reverse_arc_const_iterator;
    typedef typename 
        _Base::reverse_arc_pointer              reverse_arc_pointer;
    typedef typename 
        _Base::reverse_arc_const_pointer        reverse_arc_const_pointer;

    typedef typename _Base::capacity_type       capacity_type;
    typedef typename _Base::difference_type     difference_type;
    typedef typename _Base::size_type           size_type;


    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    optnet_mc_fs_maxflow();

    ///////////////////////////////////////////////////////////////////////
    ///  Construct a optnet_mc_fs_maxflow object with the underlying graph
    ///  being created according to the given size information.
    ///
    ///  @param  colsize  The number of nodes in a column.
    ///  @param  numcols  The number of columns.
    ///
    ///////////////////////////////////////////////////////////////////////
    optnet_mc_fs_maxflow(size_type colsize, size_type numcols);

    ///////////////////////////////////////////////////////////////////////
    ///  Construct a optnet_mc_fs_maxflow object with the underlying graph
    ///  being created according to the given size information.
    ///
    ///  @param  colsize  The number of nodes in a column.
    ///  @param  numcols  The number of columns.
    ///
    ///////////////////////////////////////////////////////////////////////
    virtual bool create(size_type colsize, size_type numcols);

    ///////////////////////////////////////////////////////////////////////
    ///  Solve the maximum-flow/minimum s-t cut problem.
    ///
    ///  @returns The maximum flow value.
    ///////////////////////////////////////////////////////////////////////
    capacity_type solve();

    ///////////////////////////////////////////////////////////////////////
    ///  Add arc(s) connecting a node to the source and/or the sink node.
    ///
    ///  @param  s    The capacity of the arc from the source node.
    ///  @param  t    The capacity of the arc to the sink node.
    ///  @param  hi   The height-index of the node.
    ///  @param  ci   The column-index of the node.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_st_arc(capacity_type s,
                           capacity_type t,
                           size_type     hi,
                           size_type     ci
                           )
    {
        m_nodes(hi, ci).cap = s - t;
        m_preflow += (s < t) ? s : t; // min(s, t)
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Add arc(s) connecting a node to the source and/or the sink node.
    ///
    ///  @param  s       The capacity of the arc from the source node.
    ///  @param  t       The capacity of the arc to the sink node.
    ///  @param  p_node  A pointer to the node.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_st_arc(capacity_type s,
                           capacity_type t,
                           node_pointer  p_node
                           )
    {
        assert(0 != p_node);
        p_node->cap = s - t;
        m_preflow += (s < t) ? s : t; // min(s, t)
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Determines if the given node is in the source set of the cut.
    ///
    ///  @param  hi   The height-index of the node.
    ///  @param  ci   The column-index of the node.
    ///
    ///  @return Returns true if the given node is the source set,
    ///          false otherwise.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline bool in_source_set(size_type hi, size_type ci) const
    {
        return (!(m_nodes(hi, ci).tag & IS_SINK))
            && (m_nodes(hi, ci).p_parent_arc != 0);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Determines if the given node is in the source set of the cut.
    ///
    ///  @param  p_node  A pointer to the node.
    ///
    ///  @return Returns true if the given node is the source set,
    ///          false otherwise.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline bool in_source_set(node_const_pointer p_node) const
    {
        return (!p_node->tag & IS_SINK)
            && (p_node->p_parent_arc != 0);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Set the initial flow value.
    ///
    ///  @param flow The initial flow value.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void set_initial_flow(const capacity_type& flow)
    {
        m_preflow = flow;
    }

private:

    typedef std::deque<node_pointer> node_queue;


    void maxflow_init();
    void maxflow_augment(node_pointer   p_s_start_node,
                         node_pointer   p_t_start_node,
                         capacity_type* p_mid_fwd_cap,
                         capacity_type* p_mid_rev_cap);
    void maxflow_adopt_source_orphan(node_pointer p_orphan);
    void maxflow_adopt_sink_orphan(node_pointer p_orphan);

    inline void activate(node_pointer p_node)
    {
        if (0 == (p_node->tag & IS_ACTIVE)) {  // Not active yet.
            m_active_nodes.push_back(p_node);
            p_node->tag |= IS_ACTIVE;
        }
    }

    inline node_pointer neighbor_node_fwd(node_pointer    p_node,
                                          difference_type shift)
    {
        return reinterpret_cast<node_pointer>
                   (reinterpret_cast<char*>(p_node) + shift);
    }
    
    inline node_pointer neighbor_node_rev(node_pointer    p_node,
                                          difference_type shift)
    {
        return reinterpret_cast<node_pointer>
                   (reinterpret_cast<char*>(p_node) - shift);
    }

    // Constants (Initialized in optnet_mc_fs_maxflow.cxx)
    static const unsigned char       PARENT_REV;// Parent arc is reverse.
    static const unsigned char       IS_ACTIVE; // The node is active.
    static const unsigned char       IS_SINK;   // The node belongs to
                                                //   the sink tree.
    static forward_arc_const_pointer TERMINAL;  // Parent is terminal node.
    static forward_arc_const_pointer ORPHAN;    // No parent.

    size_type     m_dist_id;
    node_queue    m_active_nodes, m_orphan_nodes;
    capacity_type m_preflow, m_flow;

};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_alpha/optnet_mc_fs_maxflow.cxx>
#   endif

#endif
