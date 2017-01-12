/*
 ==========================================================================
 |   
 |   $Id: optnet_pr_maxflow.cxx 127 2005-02-06 00:46:11Z kangli $
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

#ifndef ___OPTNET_PR_MAXFLOW_CXX___
#   define ___OPTNET_PR_MAXFLOW_CXX___

#   include <optnet/_pr/optnet_pr_maxflow.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
optnet_pr_maxflow<_Cap, _Tg>::optnet_pr_maxflow() :
    m_global_update_freq(0.1), m_alpha(6), m_beta(12), // Default parameters.
    m_preflow(0)
{
    //
    // Experiments show that the global relabelling heuristic does not
    // really help speeding-up the algorithm, so we do not do it 
    // as frequently as suggested by the 'classic' literature.
    //
    // One can fine-tune the global relabel frequency by adjusting the
    // m_global_update_freq variable.
    //
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
optnet_pr_maxflow<_Cap, _Tg>::optnet_pr_maxflow(size_type s0,
                                                size_type s1,
                                                size_type s2,
                                                size_type s3,
                                                size_type s4
                                                ) :
    _Base(s0, s1, s2, s3, s4),
    m_global_update_freq(0.1), m_alpha(6), m_beta(12), // Default parameters.
    m_preflow(0)
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
bool
optnet_pr_maxflow<_Cap, _Tg>::create(size_type s0,
                                     size_type s1,
                                     size_type s2,
                                     size_type s3,
                                     size_type s4
                                     )
{
    // Clear pre-calculated flow.
    m_preflow = 0;

    // Create and initialize the graph.
    return _Base::create(s0, s1, s2, s3, s4);
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
typename optnet_pr_maxflow<_Cap, _Tg>::capacity_type
optnet_pr_maxflow<_Cap, _Tg>::solve()
{
    maxflow_init();

    maxflow_compute_preflow();          // Phase 1.
    maxflow_convert_preflow_to_flow();  // Phase 2.

    maxflow_grow_source();

    return m_flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_pr_maxflow<_Cap, _Tg>::maxflow_compute_preflow()
{
    // This is the first phase of the push-relabel algorithm,
    // which is also the core phase. At the end of the first
    // phase, the excess at the sink is equal to the maximum
    // flow value.
    while (m_max_active >= m_min_active) {

        _Layer& layer = m_layers[m_max_active];
        node_pointer u = layer.p_first_active;

        if (u == &m_dummy_node) --m_max_active;
        else {
            
            pop_active(u, layer);
            maxflow_discharge(u);

            if (m_global_counter * m_global_update_freq > m_nm) {
                
                // Global update heuristic.
                maxflow_global_update();
                m_global_counter = 0;
            }
        }

    } // while
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_pr_maxflow<_Cap, _Tg>::maxflow_convert_preflow_to_flow()
{
    static const size_type WHITE = 2;
    static const size_type GRAY  = 1;
    static const size_type BLACK = 0;

    forward_arc_pointer p_fwd_arc, p_cur_arc, p_first_out_arc,
                        p_last_out_arc;
    reverse_arc_pointer p_rev_arc, p_first_in_arc,
                        p_last_in_arc;
    node_pointer        u, v, r, restart;
    node_pointer        tos = 0, bos = 0; // Top of Stack, Bottom of Stack.
    

    ///////////////////////////////////////////////////////////////////////
    // We should handle self loops here, but for our graph we don't
    // have any self-loop (do we???), so we happily skip this step.
    
    ///////////////////////////////////////////////////////////////////////
    // STEP 1: Initialize.
    for (u = &*_Base::m_nodes.begin(); u != &*_Base::m_nodes.end(); ++u) {

        u->p_current_arc = (u->tag & _Base::SPECIAL_OUT) ? 
                                u->p_first_out_arc + 1 : 
                                u->p_first_out_arc;

        u->tag   &= ~CURRENT_REV;
        u->dist   = WHITE; // white
        u->p_prev = u;
    }

    ///////////////////////////////////////////////////////////////////////
    // STEP 2: Eliminate flow cycles and topologically order the vertices.
    for (u = &*_Base::m_nodes.begin(); u != &*_Base::m_nodes.end(); ++u) {

        if (u->dist == WHITE && u->excess > 0) {
            
            r = u;
            r->dist = GRAY; // gray

            while (true) {

                // This sucks but we have to do it every time.
                if (u->tag & _Base::SPECIAL_OUT) {
                    p_first_out_arc = u->p_first_out_arc + 1;
                    p_last_out_arc
                        = reinterpret_cast<forward_arc_pointer>
                            (u->p_first_out_arc->shift);
                }
                else {
                    p_first_out_arc = u->p_first_out_arc;
                    p_last_out_arc  = (u + 1)->p_first_out_arc;
                }
                
                if (u->tag & _Base::SPECIAL_IN) {
                    p_first_in_arc  = u->p_first_in_arc + 1;
                    p_last_in_arc
                        = reinterpret_cast<reverse_arc_pointer>
                            (u->p_first_in_arc->p_fwd);
                }
                else {
                    p_first_in_arc  = u->p_first_in_arc;
                    p_last_in_arc   = (u + 1)->p_first_in_arc;
                }

                if (0 == (u->tag & CURRENT_REV)) {
                    // Skip the forward arcs.
                    u->p_current_arc = reinterpret_cast
                                        <forward_arc_pointer>
                                            (p_first_in_arc);
                    u->tag |= CURRENT_REV;
                }

                while (u->p_current_arc != 
                    reinterpret_cast<forward_arc_pointer>(p_last_in_arc)) {
                    
                    p_rev_arc = reinterpret_cast
                                    <reverse_arc_pointer>(u->p_current_arc);
                    p_fwd_arc = p_rev_arc->p_fwd;

                    if (p_fwd_arc->rev_cap > 0) {

                        v = neighbor_node_rev(u, p_fwd_arc->shift);

                        if (WHITE == v->dist) {
                            v->dist   = GRAY;
                            v->p_prev = u;
                            u         = v;
                            break; // while
                        }
                        else if (GRAY == v->dist) {

                            // Find minimum flow on the cycle.
                            capacity_type delta = p_fwd_arc->rev_cap;
                            while (true) {
                                if (v->tag & CURRENT_REV) { 
                                    p_cur_arc
                                        = reinterpret_cast
                                            <reverse_arc_pointer>
                                                (v->p_current_arc)->p_fwd;

                                    if (p_cur_arc->rev_cap < delta)
                                        delta = p_cur_arc->rev_cap;

                                    if (v != u) {
                                        v = neighbor_node_rev
                                                (v, p_cur_arc->shift);
                                    }
                                    else break;
                                }
                                else { // !CURRENT_REV
                                    p_cur_arc = v->p_current_arc;

                                    if (v != u) {
                                        v = neighbor_node_rev
                                                (v, p_cur_arc->shift);
                                    }
                                    else break;
                                }
                            } // while
                            
                            // Remove delta flow units.
                            v = u;
                            while (true) {
                                if (v->tag & CURRENT_REV) { //  CURRENT_REV
                                    p_cur_arc
                                        = reinterpret_cast
                                            <reverse_arc_pointer>
                                                (v->p_current_arc)->p_fwd;
                                    p_cur_arc->rev_cap -= delta;
                                    v = neighbor_node_rev
                                            (v, p_cur_arc->shift);
                                    if (v == u) break;
                                }
                                else {                      // !CURRENT_REV
                                    p_cur_arc = v->p_current_arc;
                                    p_cur_arc->rev_cap += delta;
                                    v = neighbor_node_fwd
                                            (u, p_cur_arc->shift);
                                    if (v == u) break;
                                }
                            } // while

                            // Back-out of DFS to the first saturated arc.
                            restart = u;
                            if (u->tag & CURRENT_REV) {     //  CURRENT_REV
                                v = neighbor_node_rev(
                                        u,
                                        reinterpret_cast
                                            <reverse_arc_pointer>
                                                (u->p_current_arc)->p_fwd->shift
                                    );
                            }
                            else {                          // !CURRENT_REV
                                v = neighbor_node_fwd
                                        (u, u->p_current_arc->shift);
                            }

                            while (v != u) {

                                if (v->tag & CURRENT_REV) {
                                    p_cur_arc 
                                        = reinterpret_cast
                                            <reverse_arc_pointer>
                                                (v->p_current_arc)->p_fwd;

                                    if (v->dist == WHITE ||
                                        p_cur_arc->rev_cap == 0) {
                                        neighbor_node_rev
                                            (v, p_cur_arc->shift)->dist = WHITE;
                                        if (v->dist != WHITE)
                                            restart = v;
                                    }
                                }
                                else {
                                    p_cur_arc = v->p_current_arc;
                                    if (v->dist == WHITE) {
                                        neighbor_node_fwd
                                            (v, p_cur_arc->shift)->dist = WHITE;
                                    }
                                }
                                
                                // Proceed to the next v.
                                if (v->tag & CURRENT_REV) {
                                    v = neighbor_node_rev
                                            (v, p_cur_arc->shift);
                                }
                                else { // !cur_rev
                                    v = neighbor_node_fwd
                                            (v, p_cur_arc->shift);
                                } // cur_rev

                            } // while (v != u)

                            if (restart != u) {
                                u = restart;
                                if (u->tag & CURRENT_REV) {
                                    u->p_current_arc
                                        = reinterpret_cast<forward_arc_pointer>
                                            (reinterpret_cast<reverse_arc_pointer>
                                                (u->p_current_arc) + 1);
                                }
                                else
                                    ++(u->p_current_arc);
                                break;
                            }
                        } // gray
                    } // if (p_fwd_arc->rev_cap > 0)

                    u->p_current_arc
                        = reinterpret_cast<forward_arc_pointer>
                              (reinterpret_cast<reverse_arc_pointer>
                                  (u->p_current_arc) + 1);

                } // while
                
                if (u->p_current_arc == 
                    reinterpret_cast<forward_arc_pointer>(p_last_in_arc)) {

                    u->dist = BLACK;

                    if (0 == bos) { bos = u; tos = u; }
                    else {
                        u->p_next = tos;
                        tos = u;
                    }

                    if (u != r) {
                        u = u->p_prev;
                        if (u->tag & CURRENT_REV) {
                            u->p_current_arc
                                = reinterpret_cast<forward_arc_pointer>
                                    (reinterpret_cast<reverse_arc_pointer>
                                        (u->p_current_arc) + 1);
                        }
                        else
                            ++(u->p_current_arc);
                    }
                    else break;
                }

            } // while (true)
        } // if (u->dist ...
    } // for (u ...
    
    ///////////////////////////////////////////////////////////////////////
    // STEP 3: Return excess flows.
    if (0 != bos) {

        for (u = tos; ; u = u->p_next) {

            // This sucks but we have to do it every time.
            if (u->tag & _Base::SPECIAL_OUT) {
                p_first_out_arc = u->p_first_out_arc + 1;
                p_last_out_arc
                    = reinterpret_cast<forward_arc_pointer>
                        (u->p_first_out_arc->shift);
            }
            else {
                p_first_out_arc = u->p_first_out_arc;
                p_last_out_arc  = (u + 1)->p_first_out_arc;
            }
            
            if (u->tag & _Base::SPECIAL_IN) {
                p_first_in_arc  = u->p_first_in_arc + 1;
                p_last_in_arc
                    = reinterpret_cast<reverse_arc_pointer>
                        (u->p_first_in_arc->p_fwd);
            }
            else {
                p_first_in_arc  = u->p_first_in_arc;
                p_last_in_arc   = (u + 1)->p_first_in_arc;
            }
            
            p_rev_arc = p_first_in_arc;

            while (u->excess > 0 &&
                   p_rev_arc != p_last_in_arc) {

                p_fwd_arc = p_rev_arc->p_fwd;

                if (p_fwd_arc->rev_cap > 0) {
                    v = neighbor_node_rev(u, p_fwd_arc->shift);
                    maxflow_push_rev(u, v, p_rev_arc);
                }

                ++p_rev_arc;
            } // while

            if (u == bos)
                break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_pr_maxflow<_Cap, _Tg>::maxflow_grow_source()
{
    node_pointer u, v;

    for (u = &*_Base::m_nodes.begin(); u != &*_Base::m_nodes.end(); ++u) {
        if (u->cap > 0 && u->excess > 0) { // Is source??
            u->tag &= ~IS_SINK;
            m_Q.push_back(u);
        }
        else {
            u->tag |=  IS_SINK;
        }
    } // for

    while (!m_Q.empty()) {

        forward_arc_pointer p_fwd_arc, p_first_out_arc, p_last_out_arc;
        reverse_arc_pointer p_rev_arc, p_first_in_arc, p_last_in_arc;

        u = m_Q.front(); m_Q.pop_front();

        if (u->tag & _Base::SPECIAL_OUT) {
            p_first_out_arc = u->p_first_out_arc + 1;
            p_last_out_arc
                = (forward_arc_pointer)(u->p_first_out_arc->shift);
        }
        else {
            p_first_out_arc = u->p_first_out_arc;
            p_last_out_arc  = (u + 1)->p_first_out_arc;
        }
        
        if (u->tag & _Base::SPECIAL_IN) {
            p_first_in_arc  = u->p_first_in_arc + 1;
            p_last_in_arc
                = (reverse_arc_pointer)(u->p_first_in_arc->p_fwd);
        }
        else {
            p_first_in_arc  = u->p_first_in_arc;
            p_last_in_arc   = (u + 1)->p_first_in_arc;
        }
            
        // Process outgoing arcs.
        for (p_fwd_arc  = p_first_out_arc;
             p_fwd_arc != p_last_out_arc;
             ++p_fwd_arc) {
            
            v = neighbor_node_fwd(u, p_fwd_arc->shift);

            if (v->tag & IS_SINK) {
                v->tag &= ~IS_SINK;
                m_Q.push_back(v);
            } // if

        } // for

        // Process incoming arcs.
        for (p_rev_arc  = p_first_in_arc;
             p_rev_arc != p_last_in_arc;
             ++p_rev_arc) {

            p_fwd_arc = p_rev_arc->p_fwd;
            v = neighbor_node_rev(u, p_fwd_arc->shift);

            if ((v->tag & IS_SINK) && p_fwd_arc->rev_cap > 0) {
                v->tag &= ~IS_SINK;
                m_Q.push_back(v);
            } // if

        } // for
    }
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_pr_maxflow<_Cap, _Tg>::maxflow_init()
{
    _Base::prepare();

    m_global_counter = 0;
    m_n              = _Base::m_nodes.size() + 2;
    m_nm             = m_alpha * m_n + _Base::m_fwd_arcs.count()
                                     + _Base::m_nodes.size();
    m_flow           = m_preflow;
    m_max_dist       = m_n - 1;
    m_min_active     = m_n;
    m_max_active     = 0;
    
    m_layers.resize(m_n);

    for (typename layer_vector::iterator l = m_layers.begin();
         l != m_layers.end();
         ++l) {
         // Empty the active and inactive node lists.
         l->p_first_active   = &m_dummy_node;
         l->p_first_inactive = &m_dummy_node;
    }

    for (node_pointer u = &*_Base::m_nodes.begin(); u != &*_Base::m_nodes.end(); ++u) {
        
        u->p_current_arc = (u->tag & _Base::SPECIAL_OUT) ? 
                                u->p_first_out_arc + 1 : 
                                u->p_first_out_arc;

        u->tag &= ~CURRENT_REV;
        u->dist = 1;        
        
        if (u->cap > 0) {
            u->excess   = u->cap;
            u->tag     &= ~IS_SINK;
            u->tag     |= IS_TERMINAL;
            push_active(u, m_layers[1]);
        }
        else { // u->cap <= 0
            
            u->excess   = 0;
            
            if (u->cap < 0) {
                u->tag |= IS_SINK;
                u->tag |= IS_TERMINAL;
            }
            else {
                u->tag &= ~IS_SINK;
                u->tag &= ~IS_TERMINAL;
            }
            push_inactive(u, m_layers[1]);
        }
    } // for
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_pr_maxflow<_Cap, _Tg>::maxflow_push_fwd(node_pointer        u,
                                               node_pointer        v,
                                               forward_arc_pointer p
                                               )
{
    capacity_type d = u->excess;
    p->rev_cap += d;
    u->excess  -= d;
    v->excess  += d;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_pr_maxflow<_Cap, _Tg>::maxflow_push_rev(node_pointer        u,
                                               node_pointer        v,
                                               reverse_arc_pointer q
                                               )
{
    forward_arc_pointer p = q->p_fwd;
    capacity_type d = u->excess < p->rev_cap ? 
                        u->excess : p->rev_cap;
    p->rev_cap -= d;
    u->excess  -= d;
    v->excess  += d;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_pr_maxflow<_Cap, _Tg>::maxflow_relabel(node_pointer u)
{
    forward_arc_pointer p_min_arc = 0;
    forward_arc_pointer p_fwd_arc, p_first_out_arc, p_last_out_arc;
    reverse_arc_pointer p_rev_arc, p_first_in_arc, p_last_in_arc;
    size_type           min_dist = m_n;
    bool                min_is_fwd = true;


    m_global_counter += m_beta;
    
    if ((u->tag & IS_TERMINAL) && (u->tag & IS_SINK)) u->dist = 1;
    else {

        u->dist = min_dist;

        if (u->tag & _Base::SPECIAL_OUT) {
            p_first_out_arc = u->p_first_out_arc + 1;
            p_last_out_arc
                = (forward_arc_pointer)(u->p_first_out_arc->shift);
        }
        else {
            p_first_out_arc = u->p_first_out_arc;
            p_last_out_arc  = (u + 1)->p_first_out_arc;
        }
    
        if (u->tag & _Base::SPECIAL_IN) {
            p_first_in_arc  = u->p_first_in_arc + 1;
            p_last_in_arc
                = (reverse_arc_pointer)(u->p_first_in_arc->p_fwd);
        }
        else {
            p_first_in_arc  = u->p_first_in_arc;
            p_last_in_arc   = (u + 1)->p_first_in_arc;
        }

        // Process outgoing arcs.
        for (p_fwd_arc = p_first_out_arc;
             p_fwd_arc < p_last_out_arc;
             ++p_fwd_arc) {
        
            node_pointer v = neighbor_node_fwd(u, p_fwd_arc->shift);

            if (v->dist < min_dist) {
                min_dist   = v->dist;
                p_min_arc  = p_fwd_arc;
                min_is_fwd = true;
            } // if

            ++m_global_counter;

        } // for

        // Process incoming arcs.
        for (p_rev_arc = p_first_in_arc;
             p_rev_arc < p_last_in_arc;
             ++p_rev_arc) {

            p_fwd_arc = p_rev_arc->p_fwd;
            node_pointer v = neighbor_node_rev(u, p_fwd_arc->shift);

            if (p_fwd_arc->rev_cap > 0 && v->dist < min_dist) {
                min_dist   = v->dist;
                p_min_arc  = reinterpret_cast<forward_arc_pointer>(p_rev_arc);
                min_is_fwd = false;
            } // if

            ++m_global_counter;

        } // for

        if (++min_dist < m_n) {

            u->dist = min_dist;
            u->p_current_arc = p_min_arc;

            if (min_is_fwd) u->tag &= ~CURRENT_REV;
            else            u->tag |=  CURRENT_REV;

            if (min_dist > m_max_dist) m_max_dist = min_dist;
        }
    }
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_pr_maxflow<_Cap, _Tg>::maxflow_discharge(node_pointer u)
{
    forward_arc_pointer p_first_out_arc, p_last_out_arc, p_fwd_arc = 0;
    reverse_arc_pointer p_first_in_arc, p_last_in_arc, p_rev_arc = 0;
    bool                found, found_rev;

    assert(u->excess > 0);

    if (u->tag & _Base::SPECIAL_OUT) {
        p_first_out_arc = u->p_first_out_arc + 1;
        p_last_out_arc
            = (forward_arc_pointer)(u->p_first_out_arc->shift);
    }
    else {
        p_first_out_arc = u->p_first_out_arc;
        p_last_out_arc  = (u + 1)->p_first_out_arc;
    }
    
    if (u->tag & _Base::SPECIAL_IN) {
        p_first_in_arc  = u->p_first_in_arc + 1;
        p_last_in_arc
            = (reverse_arc_pointer)(u->p_first_in_arc->p_fwd);
    }
    else {
        p_first_in_arc  = u->p_first_in_arc;
        p_last_in_arc   = (u + 1)->p_first_in_arc;
    }

    // Main loop.
    while (true) {

        node_pointer v;

        found     = false;
        found_rev = false;

        if ((u->tag & IS_TERMINAL) && 
            (u->tag & IS_SINK) && 
            (u->dist == 1)) {
        
            assert(u->cap < 0);
            assert(u->excess > 0);

            // Have admissible arc to the sink?
            capacity_type flow = -u->cap;
            if (u->excess < flow) flow = u->excess;

            u->excess  -= flow;
            u->cap     += flow;
            m_flow     += flow;

            if (u->cap == 0) {
                u->tag &= ~IS_TERMINAL;
                u->tag &= ~IS_SINK;
            }

            if (u->excess == 0) {
                push_inactive(u, m_layers[u->dist]);
                break;
            }

        }
                
        if (u->tag & CURRENT_REV) { // Current arc is reverse.
            for (p_rev_arc  = reinterpret_cast<reverse_arc_pointer>(u->p_current_arc);
                 p_rev_arc != p_last_in_arc;
                 ++p_rev_arc) {

                p_fwd_arc = p_rev_arc->p_fwd;
                v = neighbor_node_rev(u, p_fwd_arc->shift);

                if (p_fwd_arc->rev_cap) { // if 1
                    if (u->dist == v->dist + 1) { // Is admissible?

                        if (v->excess == 0) {
                            pop_inactive(v, m_layers[v->dist]);
                            push_active(v, m_layers[v->dist]);
                        }

                        // Push.
                        maxflow_push_rev(u, v, p_rev_arc);

                        assert(v->excess > 0);

                        if (u->excess == 0) {
                            found_rev = true;
                            found = true;
                            break;
                        }
                    }
                } // if 1
            } // for
        }
        else {

            for (p_fwd_arc  = u->p_current_arc;
                 p_fwd_arc != p_last_out_arc;
                 ++p_fwd_arc) {

                v = neighbor_node_fwd(u, p_fwd_arc->shift);

                if (u->dist == v->dist + 1) { // Is admissible?

                    if (v->excess == 0) {
                        pop_inactive(v, m_layers[v->dist]);
                        push_active(v, m_layers[v->dist]);
                    }
                    // Push.
                    maxflow_push_fwd(u, v, p_fwd_arc);

                    assert(v->excess > 0);

                    if (u->excess == 0) {
                        found = true;
                        break;
                    }
                }
            } // for

            if (!found) {

                for (p_rev_arc  = p_first_in_arc;
                     p_rev_arc != p_last_in_arc;
                     ++p_rev_arc) {

                    p_fwd_arc = p_rev_arc->p_fwd;
                    v = neighbor_node_rev(u, p_fwd_arc->shift);

                    if (p_fwd_arc->rev_cap) { // if 3
                        if (u->dist == v->dist + 1) { // Is admissible?

                            if (v->excess == 0) {
                                pop_inactive(v, m_layers[v->dist]);
                                push_active(v, m_layers[v->dist]);
                            }
                            // Push.
                            maxflow_push_rev(u, v, p_rev_arc);

                            assert(v->excess > 0);

                            if (u->excess == 0) {
                                found_rev = true;
                                found = true;
                                break;
                            }
                        }
                    }
                }

            } // if (!found)

        } // if (u->tag & CURRENT_REV)


        // Post-processing.
        _Layer&     layer = m_layers[u->dist];
        size_type   du = u->dist;

        if (!found) {

            // Node u must be relabeled.
            maxflow_relabel(u);

            if (layer.p_first_active == &m_dummy_node && 
                layer.p_first_inactive == &m_dummy_node) {
                maxflow_gap_relabel(du);
            }
            
            if (u->dist == m_n)
                break;
        }
        else {
            // Node u is no longer active.
            if (found_rev) {
                u->p_current_arc  = reinterpret_cast<forward_arc_pointer>(p_rev_arc);
                u->tag           |=  CURRENT_REV;
            }
            else {
                u->p_current_arc  = p_fwd_arc;
                u->tag           &= ~CURRENT_REV;
            }
            push_inactive(u, layer);
            break;
        }

    } // while
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_pr_maxflow<_Cap, _Tg>::maxflow_gap_relabel(size_type empty_dist)
{
    typename layer_vector::iterator l;

    // Distance of layer before the current layer.
    size_type r = empty_dist - 1;

    // Set the distance of the nodes beyond the gap to "infinity".
    for (l  = m_layers.begin() + empty_dist + 1; 
         l != m_layers.begin() + m_max_dist;
         ++l) {

        for (node_pointer u = l->p_first_inactive;
             u != &m_dummy_node;
             u = u->p_next) {

            // Update distance.
            u->dist = m_n;
        }
        
        l->p_first_inactive = &m_dummy_node;
    }

    m_max_active = r;
    m_max_dist = r;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_pr_maxflow<_Cap, _Tg>::maxflow_global_update()
{
    node_pointer u, v;

    for (size_type i = 0; i <= m_max_dist; ++i) {
        // Empty the active and inactive node lists.
        m_layers[i].p_first_active   = &m_dummy_node;
        m_layers[i].p_first_inactive = &m_dummy_node;
    }

    m_min_active = m_n;
    m_max_active = 0;
    m_max_dist   = 1;

    for (u = &*_Base::m_nodes.begin(); u != &*_Base::m_nodes.end(); ++u) {
        
        if ((u->tag & IS_TERMINAL) && 
            (u->tag & IS_SINK)) {

            assert(u->cap < 0);
            
            u->p_current_arc = (u->tag & _Base::SPECIAL_OUT) ? 
                                    u->p_first_out_arc + 1 : 
                                    u->p_first_out_arc;
            u->tag &= ~CURRENT_REV;
            u->dist = 1;

            if (u->excess > 0)
                push_active(u, m_layers[1]);
            else 
                push_inactive(u, m_layers[1]);

            m_Q.push_back(u);
        }
        else {
            u->tag &= ~IS_TERMINAL;
            u->tag &= ~IS_SINK;
            u->dist = m_n;
        }
    }

    while (!m_Q.empty()) {

        size_type           new_dist;
        forward_arc_pointer p_fwd_arc, p_first_out_arc, p_last_out_arc;
        reverse_arc_pointer p_rev_arc, p_first_in_arc, p_last_in_arc;

        u = m_Q.front();
        m_Q.pop_front();

        new_dist = u->dist + 1;
        assert(new_dist < m_n);

        if (u->tag & _Base::SPECIAL_OUT) {
            p_first_out_arc = u->p_first_out_arc + 1;
            p_last_out_arc
                = (forward_arc_pointer)(u->p_first_out_arc->shift);
        }
        else {
            p_first_out_arc = u->p_first_out_arc;
            p_last_out_arc  = (u + 1)->p_first_out_arc;
        }
        
        if (u->tag & _Base::SPECIAL_IN) {
            p_first_in_arc  = u->p_first_in_arc + 1;
            p_last_in_arc
                = (reverse_arc_pointer)(u->p_first_in_arc->p_fwd);
        }
        else {
            p_first_in_arc  = u->p_first_in_arc;
            p_last_in_arc   = (u + 1)->p_first_in_arc;
        }
            
        // Process outgoing arcs.
        for (p_fwd_arc  = p_first_out_arc;
             p_fwd_arc != p_last_out_arc;
             ++p_fwd_arc) {
            
            v = neighbor_node_fwd(u, p_fwd_arc->shift);

            if (v->dist == m_n && p_fwd_arc->rev_cap > 0) {

                v->dist = new_dist;
                v->tag |=  IS_SINK;
                v->tag &= ~CURRENT_REV;

                v->p_current_arc = (v->tag & _Base::SPECIAL_OUT) ? 
                                        v->p_first_out_arc + 1 : 
                                        v->p_first_out_arc;

                if (v->excess > 0)
                    push_active(v, m_layers[new_dist]);
                else 
                    push_inactive(v, m_layers[new_dist]);

                if (v->dist > m_max_dist)
                    m_max_dist = v->dist;

                m_Q.push_back(v);

            } // if
        } // for

        // Process incoming arcs.
        for (p_rev_arc  = p_first_in_arc;
             p_rev_arc != p_last_in_arc;
             ++p_rev_arc) {

            p_fwd_arc = p_rev_arc->p_fwd;
            v = neighbor_node_rev(u, p_fwd_arc->shift);

            if (v->dist == m_n) {

                v->dist = new_dist;
                v->tag |=  IS_SINK;
                v->tag &= ~CURRENT_REV;

                v->p_current_arc = (v->tag & _Base::SPECIAL_OUT) ? 
                                        v->p_first_out_arc + 1 : 
                                        v->p_first_out_arc;

                if (v->excess > 0)
                    push_active(v, m_layers[new_dist]);
                else 
                    push_inactive(v, m_layers[new_dist]);

                if (v->dist > m_max_dist)
                    m_max_dist = v->dist;

                m_Q.push_back(v);

            } // if
        } // for
    } // while
}

//
// Constants
//
template<typename _Cap, typename _Tg>
    const unsigned char optnet_pr_maxflow<_Cap, _Tg>::CURRENT_REV = 0x04;

template<typename _Cap, typename _Tg>
    const unsigned char optnet_pr_maxflow<_Cap, _Tg>::IS_TERMINAL = 0x02;

template<typename _Cap, typename _Tg>
    const unsigned char optnet_pr_maxflow<_Cap, _Tg>::IS_SINK     = 0x01;


} // namespace

#endif
