/*
 ==========================================================================
 |   
 |   $Id: optnet_pr_3d_multi.cxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___OPTNET_PR_3D_MULTI_CXX___
#   define ___OPTNET_PR_3D_MULTI_CXX___

#   include <optnet/_base/except.hxx>
#   include <optnet/_pr/optnet_pr_3d_multi.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <deque>


namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
optnet_pr_3d_multi<_Cost, _Cap, _Tg>::optnet_pr_3d_multi() :
    m_pcost(0)
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d_multi<_Cost, _Cap, _Tg>::create(
                                        const cost_array_base_type& cost
                                        )
{
    if (!m_graph.create(cost.size_0(), 
                        cost.size_1(), 
                        cost.size_2(),
                        cost.size_3()
                        )) {
        throw_exception(std::runtime_error(
            "optnet_pr_3d_multi::create: Could not create graph."
            ));
    }

    m_intra.resize(m_graph.size_3());

    // Set up default parameters.
    for (size_type i = 0; i < m_graph.size_3(); ++i) {
    m_intra[i].tonode[0].clear();
    m_intra[i].tonode[1].clear();
        m_intra[i].dropped   = false;
        m_intra[i].climbed   = false;
        m_intra[i].circle[0] = false;
        m_intra[i].circle[1] = false;
        m_intra[i].margin[0] = 0;
        m_intra[i].margin[1] = 0;
        m_intra[i].smooth[0] = 1;
        m_intra[i].smooth[1] = 1;
    }

    m_inter.clear();

    // Temporarily save a pointer to the cost vector.
    m_pcost = &cost;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d_multi<_Cost, _Cap, _Tg>::set_params(size_type k,
                                                 int       smooth0,
                                                 int       smooth1,
                                                 bool      circle0,
                                                 bool      circle1
                                                 )
{
    // Check arguments.
    if (k >= m_intra.size() ||  smooth0 < 0 || smooth1 < 0) {
        throw_exception(
            std::invalid_argument(
            "optnet_pr_3d_multi::set_params: Invalid argument."
        ));
    }

    m_intra[k].smooth[0] = smooth0;
    m_intra[k].smooth[1] = smooth1;
    m_intra[k].circle[0] = circle0;
    m_intra[k].circle[1] = circle1;

}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d_multi<_Cost, _Cap, _Tg>::set_relation(size_type k0,
                                                   size_type k1,
                                                   int       r0,
                                                   int       r1
                                                   )
{
    // Check arguments.
    if ((k0 >= m_graph.size_3()) || (k1 >= m_graph.size_3()) ||
        (r0 >= 0 && r1 >= 0)) {
        throw_exception(
            std::invalid_argument(
            "optnet_pr_3d_multi::set_relation: Invalid argument."
        ));
    }
    
    // The relations poses constraints that some nodes definitely cannot
    //  be on the final surfaces. 
    //
    // Note that:
    //  if r0 < 0 || r1 < 0 but not both : non-crossing case
    //  if r0 < 0 && r1 < 0              : crossing case
    //  if r0 > 0 && r1 > 0              : invalid
    //
    if (r0 * r1 <= 0) {
        
        // Negative number specifies maximum distance, possitive
        // number specifies minimum distance.
        if (r0 + r1 > 0) {
            throw_exception(
                std::invalid_argument(
                "optnet_pr_3d_multi::set_relation: Invalid surface interrelation."
            ));
        }

        // For lower margin.
        m_intra[k0].tonode[0].push_back(std::pair<size_t, int>(k1, r0));
        m_intra[k1].tonode[0].push_back(std::pair<size_t, int>(k0, r1));

        // For upper margin.
        m_intra[k1].tonode[1].push_back(std::pair<size_t, int>(k0, r0));
        m_intra[k0].tonode[1].push_back(std::pair<size_t, int>(k1, r1));

        if (r0 > 0) {
            m_intra[k1].climbed = true;
            m_intra[k0].dropped = true;
        }
        else if (r1 > 0) {
            m_intra[k0].climbed = true;
            m_intra[k1].dropped = true;
        }
    }

    // Create a new relation and 
    //  append it to the relation vector.
    _Inter relation;

    relation.k[0] = k0;
    relation.k[1] = k1;
    relation.r[0] = r0;
    relation.r[1] = r1;
    
    m_inter.push_back(relation);

}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d_multi<_Cost, _Cap, _Tg>::solve(net_base_type& net, 
                                            capacity_type* pflow
                                            )
{
    size_type       i0, i1, i2, i3;
    capacity_type   flow;

    if (net.size_0() != m_graph.size_0() || 
        net.size_1() != m_graph.size_1() ||
        net.size_2() != m_graph.size_3()
        ) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "optnet_pr_3d_multi::solve: The net size must match the graph size."
        ));
    }

    // Compute the lower and upper margin of the graph nodes in
    // each 3-D subgraph.
    get_bounds_of_subgraphs();

    // Assign the cost of graph nodes based on the input cost
    // vector. We also perform the "translation operation"
    // here to guaranttee a non-empty solution.
    transform_costs();

    // Build the arcs of the graphs.
    build_arcs();

    // Calculate max-flow/min-cut.
    flow = m_graph.solve();

    // Recover optimal surface(s) from source set.
    for (i3 = 0; i3 < m_graph.size_3(); ++i3) {
        for (i1 = 0; i1 < m_graph.size_1(); ++i1) {
            for (i0 = 0; i0 < m_graph.size_0(); ++i0) {
                for (i2 = m_intra[i3].margin[0];
                     i2 < m_graph.size_2() - m_intra[i3].margin[1];
                     ++i2) {
                    
                    // Find upper envelope.
                    if (!m_graph.in_source_set(i0, i1, i2, i3))
                        break;
                }
                net(i0, i1, i3) = (int)(i2 - 1);
            } // for i0
        } // for i1
    } // for i3

    if (0 != pflow)
        *pflow = flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d_multi<_Cost, _Cap, _Tg>::get_bounds_of_subgraphs()
{
    bool                stop[2];
    int                 r, new_bound;
    size_type           i, k, n = m_intra.size();
    std::deque<_Intra*> Q[2];

    //
    // Hopefully this is also a correct way of detecting illegal surface
    // relation specifications. 
    //

    //
    // STEP 1: Initialize.
    //
    for (i = 0; i < n; ++i) {

        _Intra* p = &m_intra[i];
        
        if (p->climbed) p->visits[0] = 0;
        else { // else 1
            Q[0].push_back(p);
            p->visits[0] = 1;
        } // else 1

        if (p->dropped) p->visits[1] = 0;
        else { // else 2
            Q[1].push_back(p);
            p->visits[1] = 1;
        } // else 2

    } // for

    //
    // STEP 2: Compute the margin of the graphs.
    //
    if (Q[0].empty() || Q[1].empty()) {
        // The specified relations are invalid.
        throw_exception(
            std::invalid_argument(
            "optnet_pr_3d_multi::get_bounds_of_subgraphs: Invalid surface interrelation."
        ));
    }
    else {
        stop[0] = false;
        stop[1] = false;
        
        //
        // Lower margin.
        do {
            _Intra* p = Q[0].front();
            Q[0].pop_front();

            for (i = 0; i < p->tonode[0].size(); ++i) {
                k = p->tonode[0][i].first;
                r = p->tonode[0][i].second;
                new_bound = p->margin[0] + r;
                if (m_intra[k].margin[0] < new_bound) {
                    if (++m_intra[k].visits[0] < (int)n) {
                        m_intra[k].margin[0] = new_bound;
                        Q[0].push_back(&m_intra[k]);
                    }
                    else {
                        stop[0] = true;
                        break;
                    } // else
                } // if
            } // for

        } // do
        while (!Q[0].empty() && !stop[0]);

        //
        // Upper margin.
        do {
            _Intra* p = Q[1].front();
            Q[1].pop_front();

            for (i = 0; i < p->tonode[1].size(); ++i) {
                k = p->tonode[1][i].first;
                r = p->tonode[1][i].second;
                new_bound = p->margin[1] + r;
                if (m_intra[k].margin[1] < new_bound) {
                    if (++m_intra[k].visits[1] < (int)n) {
                        m_intra[k].margin[1] = new_bound;
                        Q[1].push_back(&m_intra[k]);
                    }
                    else {
                        stop[1] = true;
                        break;
                    } // else
                } // if
            } // for

        } // do
        while (!Q[1].empty() && !stop[1]);

        if (stop[0] || stop[1]) {
            // There is a loop in the relations specified.
            throw_exception(
                std::invalid_argument(
                "optnet_pr_3d_multi::get_bounds_of_subgraphs: Invalid surface interrelation."
            ));
        } // if
    }
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d_multi<_Cost, _Cap, _Tg>::transform_costs()
{
    assert(m_pcost != 0);

    size_type            i0, i1, i2, i3;
    capacity_type        cap_sum0 = 0;

    const size_type& s0 = m_pcost->size_0();
    const size_type& s1 = m_pcost->size_1();
    const size_type& s2 = m_pcost->size_2();
    const size_type& s3 = m_pcost->size_3();

    m_graph.set_initial_flow(0);

    //
    // Calculate the cost sum of the zero-plane.
    //
    for (i3 = 0; i3 < s3; ++i3) {
        {
            i2 = m_intra[i3].margin[0]; // The lower bound. 
                                        // Nodes below the lower bound
                                        //   will be ignored.
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {
                    cap_sum0 += (capacity_type)(*m_pcost)(i0, i1, i2, i3);
                } // for i0
            } // for i1
        }
    } // for i3

    //
    // Construct the s-t graph "G_st".
    //
    for (i3 = 0; i3 < s3; ++i3) {
        // Nodes below margin[0] and above s2-margin[1] will
        // be ignored.
        for (i2 = m_intra[i3].margin[0] + 1;
             i2 < s2 - m_intra[i3].margin[1];
             ++i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {

                    capacity_type cap
                        = (capacity_type)(*m_pcost)(i0, i1, i2, i3)
                        - (capacity_type)(*m_pcost)(i0, i1, i2 - 1, i3);

                    if (cap >= 0) // non-negative -> connect to t
                        m_graph.add_st_arc(0, +cap, i0, i1, i2, i3);
                    else          // negative     -> connect to s
                        m_graph.add_st_arc(-cap, 0, i0, i1, i2, i3);

                } // for i0
            } // for i1
        } // for i2
    } // for i3

    // Zero-plane.
    if (cap_sum0 >= 0) {
        for (i3 = 0; i3 < s3; ++i3) {
            i2 = m_intra[i3].margin[0];
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {

                    capacity_type cap
                        = (capacity_type)(*m_pcost)(i0, i1, i2, i3);

                    // Translation operation.
                    if (i3 == 0 && i1 == 0 && i0 == 0)
                        cap -= (cap_sum0 + 1);

                    if (cap >= 0) // non-negative -> connect to t
                        m_graph.add_st_arc(0, +cap, i0, i1, i2, i3);
                    else          // negative     -> connect to s
                        m_graph.add_st_arc(-cap, 0, i0, i1, i2, i3);

                } // for i0
            } // for i1
        } // for i3
    }
    else {
        for (i3 = 0; i3 < s3; ++i3) {
            i2 = m_intra[i3].margin[0];
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {

                    capacity_type cap
                        = (capacity_type)(*m_pcost)(i0, i1, i2, i3);

                    if (cap >= 0) // non-negative -> connect to t
                        m_graph.add_st_arc(0, +cap, i0, i1, i2, i3);
                    else          // negative     -> connect to s
                        m_graph.add_st_arc(-cap, 0, i0, i1, i2, i3);

                } // for i0
            } // for i1
        } // for i3
    }

}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d_multi<_Cost, _Cap, _Tg>::build_arcs()
{
    int ii, i0, i1, i2, i3, s0, s1, s2, s3;

    s0 = (int)m_graph.size_0();
    s1 = (int)m_graph.size_1();
    s2 = (int)m_graph.size_2();
    s3 = (int)m_graph.size_3();

    // Clear all constructed arcs.
    m_graph.clear_arcs();

    //
    // Construct intra-surface arcs here.
    //

    for (i3 = 0; i3 < s3; ++i3) {

        const bool& circle0 = m_intra[i3].circle[0];
        const bool& circle1 = m_intra[i3].circle[1];
        const int & bounds0 = m_intra[i3].margin[0];
        const int & bounds1 = m_intra[i3].margin[1];
        const int & smooth0 = m_intra[i3].smooth[0];
        const int & smooth1 = m_intra[i3].smooth[1];

        // -- Intra-column (vertical) arcs.
        for (i2 = s2 - bounds1 - 1;
             i2 > bounds0;
             --i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1,          i2 - 1,          i3);
                } // for i0 
            } // for i1
        } // for i2

        // -- Inter-column arcs (dir-0).
        for (i2 = s2 - bounds1 - 1;
             i2 > smooth0 + bounds0;
             --i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 1; i0 < s0 - 1; ++i0) {
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0 - 1,      i1,          i2 - smooth0,    i3);
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0 + 1,      i1,          i2 - smooth0,    i3);
                } {
                    m_graph.add_arc(0,           i1,          i2,          i3,          1,           i1,          i2 - smooth0,    i3);
                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          s0 - 2,      i1,          i2 - smooth0,    i3);
                } // for i0
            } // for i1
        } // for i2

        // -- Inter-column arcs (dir-1).
        for (i2 = s2 - bounds1 - 1;
             i2 > smooth1 + bounds0;
             --i2) {
            for (i0 = 0; i0 < s0; ++i0) {
                for (i1 = 1; i1 < s1 - 1; ++i1) {
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1 - 1,      i2 - smooth1,    i3);
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1 + 1,      i2 - smooth1,    i3);
                } {
                    m_graph.add_arc(i0,          0,           i2,          i3,          i0,          1,           i2 - smooth1,    i3);
                    m_graph.add_arc(i0,          s1 - 1,      i2,          i3,          i0,          s1 - 2,      i2 - smooth1,    i3);
                } // for i0
            } // for i1
        } // for i2

        // -- The zero-plane
        { // i2 = bounds0
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 1; i0 < s0 - 1; ++i0) {
                    m_graph.add_arc(i0,          i1,          bounds0,     i3,          i0 - 1,      i1,          bounds0,         i3);
                    m_graph.add_arc(i0,          i1,          bounds0,     i3,          i0 + 1,      i1,          bounds0,         i3);
                } {
                    m_graph.add_arc(0,           i1,          bounds0,     i3,          1,           i1,          bounds0,         i3);
                    m_graph.add_arc(s0 - 1,      i1,          bounds0,     i3,          s0 - 2,      i1,          bounds0,         i3);
                } // for i0
            } // for i1
        }

        { // i2 = bounds0
            for (i0 = 0; i0 < s0; ++i0) {
                for (i1 = 1; i1 < s1 - 1; ++i1) {
                    m_graph.add_arc(i0,          i1,          bounds0,     i3,          i0,          i1 - 1,      bounds0,         i3);
                    m_graph.add_arc(i0,          i1,          bounds0,     i3,          i0,          i1 + 1,      bounds0,         i3);
                } {
                    m_graph.add_arc(i0,          0,           bounds0,     i3,          i0,          1,           bounds0,         i3);
                    m_graph.add_arc(i0,          s1 - 1,      bounds0,     i3,          i0,          s1 - 2,      bounds0,         i3);
                } // for i0
            } // for i1
        }

        // -- Circular graph connections.
        if (circle0) {
            for (i2 = s2 - bounds1 - 1;
                 i2 > smooth0 + bounds0;
                 --i2) {
                for (i1 = 0; i1 < s1; ++i1) {
                    m_graph.add_arc(0,           i1,          i2,          i3,          s0 - 1,      i1,          i2 - smooth0,    i3);
                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          0,           i1,          i2 - smooth0,    i3);
                }
            } { // i2 = bounds0
                for (i1 = 0; i1 < s1; ++i1) {
                    m_graph.add_arc(0,           i1,          bounds0,     i3,          s0 - 1,      i1,          bounds0,         i3);
                    m_graph.add_arc(s0 - 1,      i1,          bounds0,     i3,          0,           i1,          bounds0,         i3);
                }
            }
        } // if

        if (circle1) {
            for (i2 = s2 - bounds1 - 1;
                 i2 > smooth1 + bounds0;
                 --i2) {
                for (i0 = 0; i0 < s0; ++i0) {
                    m_graph.add_arc(i0,          0,           i2,          i3,          i0,          s1 - 1,      i2 - smooth1,    i3);
                    m_graph.add_arc(i0,          s1 - 1,      i2,          i3,          i0,          0,           i2 - smooth1,    i3);
                }                                                          
            } { // i2 = 0                                                  
                for (i0 = 0; i0 < s0; ++i0) {                              
                    m_graph.add_arc(i0,          0,           bounds0,     i3,          i0,          s1 - 1,      bounds0,         i3);
                    m_graph.add_arc(i0,          s1 - 1,      bounds0,     i3,          i0,          0,           bounds0,         i3);
                }
            }
        } // if

    } // for i3
    
    // Construct inter-surface arcs here.
    for (i3 = 0; i3 < (int)m_inter.size(); ++i3) {

        const size_type& k0 = m_inter[i3].k[0];
        const size_type& k1 = m_inter[i3].k[1];
        const int&       r0 = m_inter[i3].r[0];
        const int&       r1 = m_inter[i3].r[1];

        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {

                // k0 --> k1
                for (i2 = m_intra[k0].margin[0];
                     i2 < s2 - m_intra[k0].margin[1];
                     ++i2) {
                    
                    ii = i2 + r0;
                    if (ii <= m_intra[k1].margin[0] || 
                        ii >= s2 - m_intra[k1].margin[1]) continue;
                    m_graph.add_arc(i0, i1, i2, k0, i0, i1, ii, k1);
                }

                // k1 --> k0
                for (i2 = m_intra[k1].margin[0];
                     i2 < s2 - m_intra[k1].margin[1];
                     ++i2) {
                    
                    ii = i2 + r1;
                    if (ii <= m_intra[k0].margin[0] || 
                        ii >= s2 - m_intra[k0].margin[1]) continue;
                    m_graph.add_arc(i0, i1, i2, k1, i0, i1, ii, k0);
                }

            } // for i0
        } // for i1
    } // for i3

    // Interconnect "zero-plane"s.
    for (i3 = 1; i3 < s3; ++i3) {

        m_graph.add_arc(0, 0, m_intra[i3    ].margin[0], i3,     
                        0, 0, m_intra[i3 - 1].margin[0], i3 - 1
                        );
        m_graph.add_arc(0, 0, m_intra[i3 - 1].margin[0], i3 - 1,
                        0, 0, m_intra[i3    ].margin[0], i3
                        );
    }

}


} // namespace

#endif
