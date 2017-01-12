/*
 ==========================================================================
 |   
 |   $Id: optnet_pr_maxflow.hxx 21 2005-01-14 15:52:31Z kangli $
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

      This file implements Goldberg's push-relabel max-flow/min-cut algo-
      rithm on a forward-star graph.
      
      The implementation makes one or more simplifications that are only
      applicable to the Optimal-Net algorithm, as listed below:

       * The residual capacity of non-st arcs is assumed to be +infinity
         and will not be decremented.


  - Reference(s):

    [1] Boris V. Cherkassky and Andrew V. Goldberg
        On Implementing Push-Relabel Method for the Maximum-Flow Problem
        Algorithmica, vol. 19, no. 4, pp 390-410, September, 1994.
        URL: http://citeseer.ist.psu.edu/cherkassky94implementing.html
    [2] Andrew V. Goldberg and Robert E. Tarjan
        A New Approach to the Maximum-Flow Problem
        Journal of the ACM (JACM), vol. 35, issue 4, pp 921-940, 1988
        URL: http://portal.acm.org/citation.cfm?id=61051&dl=ACM&coll=portal
 ==========================================================================
 */

#ifndef ___OPTNET_PR_MAXFLOW_HXX___
#   define ___OPTNET_PR_MAXFLOW_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/_pr/graph_pr.hxx>
#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <vector>
#   include <deque>


namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class optnet_pr_maxflow
///  @brief Implementation of the Push-Relabel max-flow algorithm on
///         a forward-star represented graph.
///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg = net_f_xy>
class optnet_pr_maxflow
    : public graph_pr<_Cap, _Tg>
{
    typedef graph_pr<_Cap, _Tg> _Base;

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
    optnet_pr_maxflow();

    ///////////////////////////////////////////////////////////////////////
    ///  Construct a optnet_pr_maxflow object with the underlying graph
    ///  being created according to the given size information.
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
    optnet_pr_maxflow(size_type s0,
                      size_type s1,
                      size_type s2,
                      size_type s3 = 1,
                      size_type s4 = 1
                      );
    
    ///////////////////////////////////////////////////////////////////////
    ///  Create a graph for the maximum-flow problem according to the
    ///  given size information.
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
    ///  @param  i0   The first  index of the node. 
    ///  @param  i1   The second index of the node. 
    ///  @param  i2   The third  index of the node. 
    ///  @param  i3   The fourth index of the node (default: 0)
    ///  @param  i4   The fifth  index of the node (default: 0)
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_st_arc(capacity_type s,
                           capacity_type t,
                           size_type     i0,
                           size_type     i1,
                           size_type     i2,
                           size_type     i3 = 0,
                           size_type     i4 = 0
                           )
    {
        m_nodes(i0, i1, i2, i3, i4).cap = s - t;
        m_preflow += (s < t) ? s : t;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Determines if the given node is in the source set of the graph.
    ///
    ///  @param  i0   The first  index of the node. 
    ///  @param  i1   The second index of the node. 
    ///  @param  i2   The third  index of the node. 
    ///  @param  i3   The fourth index of the node (default: 0). 
    ///  @param  i4   The fifth  index of the node (default: 0). 
    ///
    ///  @return Returns true if the given node is the source set,
    ///          false otherwise.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline bool in_source_set(size_type i0,
                              size_type i1,
                              size_type i2,
                              size_type i3 = 0,
                              size_type i4 = 0
                              ) const
    {
        return (0 == (m_nodes(i0, i1, i2, i3, i4).tag & IS_SINK));
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

    // Internal convenient types.
    typedef std::deque<node_pointer>    node_queue;

    struct  _Layer
    {
        node_pointer                    p_first_active;
        node_pointer                    p_first_inactive;
    };

    typedef std::vector<_Layer>         layer_vector;


    // Helper.
    void      maxflow_compute_preflow();            // First phase.
    void      maxflow_convert_preflow_to_flow();    // Second phase.
    void      maxflow_grow_source();

    // Core functions.
    void      maxflow_init();

    void      maxflow_push_fwd(node_pointer        u,
                               node_pointer        v,
                               forward_arc_pointer p
                               );
    void      maxflow_push_rev(node_pointer        u,
                               node_pointer        v,
                               reverse_arc_pointer q
                               );
    void      maxflow_relabel(node_pointer u);
    void      maxflow_discharge(node_pointer u);

    // Heuristics.
    void      maxflow_gap_relabel(size_type empty_dist);
    void      maxflow_global_update();

    ///////////////////////////////////////////////////////////////////////
    //
    inline void push_active  (node_pointer u, _Layer& layer)
    {
        u->p_next = layer.p_first_active;
        layer.p_first_active = u;

        if (u->dist > m_max_active) m_max_active = u->dist;
        if (u->dist < m_min_active) m_min_active = u->dist;
    }

    ///////////////////////////////////////////////////////////////////////
    //
    inline void pop_active   (node_pointer u, _Layer& layer)
    {
        assert(layer.p_first_active == u);
        layer.p_first_active = u->p_next;
    }

    ///////////////////////////////////////////////////////////////////////
    //
    inline void push_inactive(node_pointer u, _Layer& layer)
    {
        node_pointer p_next    = layer.p_first_inactive;
        p_next->p_prev         = u;
        u->p_next              = p_next;
        u->p_prev              = &m_dummy_node;
        layer.p_first_inactive = u;
    }

    ///////////////////////////////////////////////////////////////////////
    //
    inline void pop_inactive (node_pointer u, _Layer& layer)
    {
        node_pointer p_prev, p_next = u->p_next;
        if (layer.p_first_inactive == u) {
            layer.p_first_inactive = p_next;
            p_next->p_prev = &m_dummy_node;
        }
        else {
            p_prev = u->p_prev;
            p_prev->p_next = p_next;
            p_next->p_prev = p_prev;
        }
    }

    ///////////////////////////////////////////////////////////////////////
    //
    inline node_pointer neighbor_node_fwd(node_pointer    p_node,
                                          difference_type shift)
    {
        return reinterpret_cast<node_pointer>
                   (reinterpret_cast<char*>(p_node) + shift);
    }
    
    ///////////////////////////////////////////////////////////////////////
    //
    inline node_pointer neighbor_node_rev(node_pointer    p_node,
                                          difference_type shift)
    {
        return reinterpret_cast<node_pointer>
                   (reinterpret_cast<char*>(p_node) - shift);
    }


    ///////////////////////////////////////////////////////////////////////
    static const unsigned char CURRENT_REV;
    static const unsigned char IS_TERMINAL;
    static const unsigned char IS_SINK;

    double              m_global_update_freq;
    node_type           m_dummy_node;
    size_type           m_alpha, m_beta, m_n, m_nm, m_global_counter,
                        m_max_dist, m_max_active, m_min_active;
    capacity_type       m_preflow, m_flow;
    layer_vector        m_layers;
    node_queue          m_Q;
};


} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_pr/optnet_pr_maxflow.cxx>
#   endif

#endif
