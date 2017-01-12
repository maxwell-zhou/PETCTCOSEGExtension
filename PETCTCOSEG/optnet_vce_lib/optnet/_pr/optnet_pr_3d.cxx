/*
 ==========================================================================
 |   
 |   $Id: optnet_pr_3d.cxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___OPTNET_PR_3D_CXX___
#   define ___OPTNET_PR_3D_CXX___

#   include <optnet/_base/except.hxx>
#   include <optnet/_pr/optnet_pr_3d.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
optnet_pr_3d<_Cost, _Cap, _Tg>::optnet_pr_3d() :
    m_pcost(0)
{
    // Initialize smoothness and circleity parameters.
    m_smooth[0] = 1;
    m_smooth[1] = 1;
    m_circle[0] = false;
    m_circle[1] = false;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d<_Cost, _Cap, _Tg>::create(
                                    const cost_array_base_type& cost
                                    )
{
    if (!m_graph.create(cost.size_0(), 
                        cost.size_1(), 
                        cost.size_2()
                        )) {
        throw_exception(std::runtime_error(
            "optnet_pr_3d::create: Could not create graph."
            ));
    }

    // Temporarily save a pointer to the cost vector.
    m_pcost = &cost;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d<_Cost, _Cap, _Tg>::set_params(int  smooth0,
                                           int  smooth1,
                                           bool circle0,
                                           bool circle1
                                           )
{
    if (smooth0 < 0 || smooth1 < 0) {
        throw_exception(std::invalid_argument(
            "optnet_pr_3d::set_params: Invalid argument."
            ));
    }

    m_smooth[0] = smooth0;
    m_smooth[1] = smooth1;
    m_circle[0] = circle0;
    m_circle[1] = circle1;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d<_Cost, _Cap, _Tg>::solve(net_base_type& net,
                                      capacity_type* pflow
                                      )
{
    size_type       i0, i1, i2;
    capacity_type   flow;

    if (net.size_0() != m_graph.size_0() || 
        net.size_1() != m_graph.size_1()
        ) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "optnet_pr_3d::solve: The output surface size must match the graph size."
            ));
    }

    // Assign the cost of graph nodes based on the input cost
    // vector. We also perform the "translation operation"
    // here to guaranttee a non-empty solution.
    transform_costs();

    // Build the arcs of the graphs.
    build_arcs();

    // Calculate max-flow/min-cut.
    flow = m_graph.solve();

    // Recover optimal surface(s) from source set.
    for (i1 = 0; i1 < m_graph.size_1(); ++i1) {
        for (i0 = 0; i0 < m_graph.size_0(); ++i0) {
            for (i2 = 0; i2 < m_graph.size_2(); ++i2) {
                
                // Find upper envelope.
                if (!m_graph.in_source_set(i0, i1, i2))
                    break;
            }
            net(i0, i1, 0) = (int)(i2 - 1);
        }
    }

    if (0 != pflow)
        *pflow = flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d<_Cost, _Cap, _Tg>::transform_costs()
{
    assert(m_pcost != 0);

    size_type            i0, i1, i2;
    capacity_type        cap_sum0 = 0;

    const size_type& s0 = m_pcost->size_0();
    const size_type& s1 = m_pcost->size_1();
    const size_type& s2 = m_pcost->size_2();

    m_graph.set_initial_flow(0);

    //
    // Calculate the cost sum of the zero-plane.
    //
    { // s2 = 0
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
                cap_sum0 += (capacity_type)(*m_pcost)(i0, i1, 0);
            }
        }
    }

    //
    // Construct the s-t graph "G_st".
    //
    for (i2 = 1; i2 < s2; ++i2) {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {

                capacity_type cap = (capacity_type)(*m_pcost)(i0, i1, i2)
                                  - (capacity_type)(*m_pcost)(i0, i1, i2 - 1);

                if (cap >= 0) // non-negative -> connect to t
                    m_graph.add_st_arc(0, +cap, i0, i1, i2);
                else          // negative     -> connect to s
                    m_graph.add_st_arc(-cap, 0, i0, i1, i2);

            } // for i0
        } // for i1
    } // for i2

    // s2 = 0
    if (cap_sum0 >= 0) {

        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {

                capacity_type cap = (capacity_type)(*m_pcost)(i0, i1, 0);

                // Perform translation operation.
                if ((0 == i1) && (0 == i0)) cap -= (cap_sum0 + 1);

                if (cap >= 0) // non-negative -> connect to t
                    m_graph.add_st_arc(0, +cap, i0, i1, 0);
                else          // negative     -> connect to s
                    m_graph.add_st_arc(-cap, 0, i0, i1, 0);

            } // for i0
        } // for i1

    }
    else {
        
        // Need not perform the translation operation.

        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {

                capacity_type cap = (capacity_type)(*m_pcost)(i0, i1, 0);

                if (cap >= 0) // non-negative -> connect to t
                    m_graph.add_st_arc(0, +cap, i0, i1, 0);
                else          // negative     -> connect to s
                    m_graph.add_st_arc(-cap, 0, i0, i1, 0);

            } // for i0
        } // for i1

    }
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_pr_3d<_Cost, _Cap, _Tg>::build_arcs()
{
    int i0, i1, i2, s0, s1, s2;

    s0 = (int)m_graph.size_0();
    s1 = (int)m_graph.size_1();
    s2 = (int)m_graph.size_2();

    // Clear all constructed arcs.
    m_graph.clear_arcs();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Construct graph arcs based on the given parameters.
    // -- Sorry, the lines here are a bit too long. ;-)

    // Intra-column (vertical) arcs.
    for (i2 = s2 - 1; i2 > 0; --i2) {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
                m_graph.add_arc(i0,          i1,          i2,          i0,          i1,          i2 - 1          );
            }
        }
    }

    // Inter-column arcs (dir-0).
    for (i2 = s2 - 1; i2 > m_smooth[0]; --i2) {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 1; i0 < s0 - 1; ++i0) {
                m_graph.add_arc(i0,          i1,          i2,          i0 - 1,      i1,          i2 - m_smooth[0]);
                m_graph.add_arc(i0,          i1,          i2,          i0 + 1,      i1,          i2 - m_smooth[0]);
            } {
                m_graph.add_arc(0,           i1,          i2,          1,           i1,          i2 - m_smooth[0]);
                m_graph.add_arc(s0 - 1,      i1,          i2,          s0 - 2,      i1,          i2 - m_smooth[0]);
            }
        }
    }

    // Inter-column arcs (dir-1).
    for (i2 = s2 - 1; i2 > m_smooth[1]; --i2) {
        for (i0 = 0; i0 < s0; ++i0) {
            for (i1 = 1; i1 < s1 - 1; ++i1) {
                m_graph.add_arc(i0,          i1,          i2,          i0,          i1 - 1,      i2 - m_smooth[1]);
                m_graph.add_arc(i0,          i1,          i2,          i0,          i1 + 1,      i2 - m_smooth[1]);
            } {
                m_graph.add_arc(i0,          0,           i2,          i0,          1,           i2 - m_smooth[1]);
                m_graph.add_arc(i0,          s1 - 1,      i2,          i0,          s1 - 2,      i2 - m_smooth[1]);
            }
        }
    }

    // The zero-plane.
    {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 1; i0 < s0 - 1; ++i0) {
                m_graph.add_arc(i0,          i1,          0,           i0 - 1,      i1,          0               );
                m_graph.add_arc(i0,          i1,          0,           i0 + 1,      i1,          0               );
            } {
                m_graph.add_arc(0,           i1,          0,           1,           i1,          0               );
                m_graph.add_arc(s0 - 1,      i1,          0,           s0 - 2,      i1,          0               );
            }
        }

        for (i0 = 0; i0 < s0; ++i0) {
            for (i1 = 1; i1 < s1 - 1; ++i1) {
                m_graph.add_arc(i0,          i1,          0,           i0,          i1 - 1,      0               );
                m_graph.add_arc(i0,          i1,          0,           i0,          i1 + 1,      0               );
            } {
                m_graph.add_arc(i0,          0,           0,           i0,          1,           0               );
                m_graph.add_arc(i0,          s1 - 1,      0,           i0,          s1 - 2,      0               );
            }
        }
    }

    // Circular graph connections.
    if (m_circle[0]) {
        for (i2 = s2 - 1; i2 > m_smooth[0]; --i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                m_graph.add_arc(0,           i1,          i2,          s0 - 1,      i1,          i2 - m_smooth[0]);
                m_graph.add_arc(s0 - 1,      i1,          i2,          0,           i1,          i2 - m_smooth[0]);
            }
        } { // i2 = 0
            for (i1 = 0; i1 < s1; ++i1) {
                m_graph.add_arc(0,           i1,          0,           s0 - 1,      i1,          0               );
                m_graph.add_arc(s0 - 1,      i1,          0,           0,           i1,          0               );
            }
        }
    } // if

    if (m_circle[1]) {
        for (i2 = s2 - 1; i2 > m_smooth[1]; --i2) {
            for (i0 = 0; i0 < s0; ++i0) {
                m_graph.add_arc(i0,          0,           i2,          i0,          s1 - 1,      i2 - m_smooth[1]);
                m_graph.add_arc(i0,          s1 - 1,      i2,          i0,          0,           i2 - m_smooth[1]);
            }
        } { // i2 = 0
            for (i0 = 0; i0 < s0; ++i0) {
                m_graph.add_arc(i0,          0,           0,           i0,          s1 - 1,      0               );
                m_graph.add_arc(i0,          s1 - 1,      0,           i0,          0,           0               );
            }
        }
    } // if
}

} // namespace

#endif
