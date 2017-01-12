/*
 ==========================================================================
 |   
 |   $Id: optnet_mc_fs_3d.cxx 141 2005-02-07 05:59:03Z kangli $
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

#ifndef ___OPTNET_MC_FS_3D_CXX___
#   define ___OPTNET_MC_FS_3D_CXX___

#   include <optnet/_alpha/optnet_mc_fs_3d.hxx>
#   include <optnet/_utils/interp_bspline3.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
optnet_mc_fs_3d<_Cost, _Cap, _Tg>::optnet_mc_fs_3d() :
    m_sc(1), m_pgs(0)
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void optnet_mc_fs_3d<_Cost, _Cap, _Tg>::create(
                                        const cost_array_base_type& cost,
                                        const point_type&           orig,
                                        graphspec_type&             spec,
                                        real_value_type             rmin,
                                        real_value_type             rmax,
                                        size_type                   npc,
                                        int                         sc
                                        )
{
    using namespace optnet::utils;

    // Interpolator type.
    typedef interp_bspline3<cost_type,
                            real_value_type,
                            _Tg>            interp_type;

    size_type s0, s1, s2;
    int       i, j, numcols, numadjs;

    assert(npc >= 2); // at least two nodes per column
    assert(sc  >= 0);

    s0 = cost.size_0();
    s1 = cost.size_1();
    s2 = cost.size_2();

    if (orig[0] >= s0 || orig[1] >= s1 || orig[2] >= s2) {
        throw_exception(std::invalid_argument(
            "optnet_mc_fs_3d::create: Origin out of bound."
            ));
    }

    if (rmin < 0) rmin = 0;
    if (rmax < 0) {
        // Automatically determine the upper physical column
        // boundary based on the origin.
        real_value_type r0, r1;

        r0 = min(orig[0], s0 - orig[0]);
        r1 = min(orig[1], s1 - orig[1]);
        r0 = min(r0, r1);
        r1 = min(orig[2], s2 - orig[2]);
        r0 = min(r0, r1);

        if (r0 < 1 || r0 < rmin) {
            throw_exception(std::invalid_argument(
                "optnet_mc_fs_3d::create: Invalid upper physical column boundary."
                ));
        } // if

        rmax = r0;
    } // if

    numcols = (int)spec.num_columns();
    numadjs = (int)spec.num_adjacencies();

    // ====================================================================
    // Create multicolumn graph.
    m_graph.create(npc, numcols);

    // ====================================================================
    // Create intra-column arcs.
    for (i = 0; i < numcols; ++i) {
        for (j = 1; j < (int)npc; ++j) {
            // Add arcs with infinite capacities.
            m_graph.add_arc(j, i, j - 1, i);
        } // for j
    } // for i

    // ====================================================================
    // Create inter-column arcs.
    for (i = 0; i < numadjs; ++i) {

        typename graphspec_type::adjacency_const_reference
                    adj = spec.get_adjacency(i);
        int c0 = adj.id[0];
        int c1 = adj.id[1];

        // Add arcs with infinite capacities.
        for (j = sc + 1; j < (int)npc; ++j) {
            m_graph.add_arc(j, c0, j - sc, c1);
            m_graph.add_arc(j, c1, j - sc, c0);
        } // for j

        m_graph.add_arc(0, c0, 0, c1);
        m_graph.add_arc(0, c1, 0, c0);

    } // for i

    // ====================================================================
    // Create t-links.

    // m_vfact[j] * col[i].pos is the true physical position of node j of
    // column i in the volume.
    m_vfact.resize(npc);

    real_value_type frac =  (rmax - rmin) / (real_value_type)(npc - 1);
    for (j = 0; j < (int)npc; ++j)
        m_vfact[j] = j * frac + rmin;

    // Initialize interpolator.
    interp_type interp3(cost, 1); // trilinear interpolation

    // Loop through all columns.
    for (i = 0; i < numcols; ++i) {
        typename graphspec_type::column_const_reference
                    col = spec.get_column(i);
        real_value_type cv0, cv1, pos[3];
        capacity_type cap;

        pos[0] = real_value_type(col.pos[0] * m_vfact[0]) + orig[0];
        pos[1] = real_value_type(col.pos[1] * m_vfact[0]) + orig[1];
        pos[2] = real_value_type(col.pos[2] * m_vfact[0]) + orig[2];
                
        cv0 = interp3.interp(pos[0], pos[1], pos[2]);

        // Add t-link with capacity 1 to the source.
        m_graph.add_st_arc(1, 0, 0, i);

        // Loop through all nodes in the current column.
        for (j = 1; j < (int)npc; ++j) {
            pos[0] = real_value_type(col.pos[0] * m_vfact[j]) + orig[0];
            pos[1] = real_value_type(col.pos[1] * m_vfact[j]) + orig[1];
            pos[2] = real_value_type(col.pos[2] * m_vfact[j]) + orig[2];

            cv1 = interp3.interp(pos[0], pos[1], pos[2]);
            cap = (capacity_type)(cv1 - cv0);
            
            if (cap >= 0) // non-negative -> connect to t
                m_graph.add_st_arc(0, +cap, j, i);
            else          // negative     -> connect to s
                m_graph.add_st_arc(-cap, 0, j, i);

            cv0 = cv1;
        } // for j
    } // for i

    // Save the boundary and origin.
    spec.set_boundary(0.0, (double)(s0 - 1), 
                      0.0, (double)(s1 - 1),
                      0.0, (double)(s2 - 1));
    spec.set_origin(orig);

    // Save the parameters.
    m_sc = sc;
    m_pgs = &spec;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void optnet_mc_fs_3d<_Cost, _Cap, _Tg>::solve(capacity_type* pflow)
{
    capacity_type   flow;
    size_type       i, j, npc, numcols;

    if (NULL == m_pgs) {
        throw_exception(std::runtime_error(
            "optnet_mc_fs_3d::solve: The graph was not created properly."
            ));
    }

    numcols = m_graph.num_columns();
    npc     = m_graph.nodes_per_column();

    assert(numcols != 0);
    assert(npc != 0);

    // Calculate max-flow/min-cut.
    flow = m_graph.solve();
    
    for (i = 0; i < numcols; ++i) {
        // Find upper envelope for each column.
        for (j = 0; j < npc; ++j)
            if (!m_graph.in_source_set(j, i)) break;

        if (j - 1 < npc)
            m_pgs->set_fact(i, 0, m_vfact[j - 1]);
        else
            m_pgs->set_fact(i, 0, -1);
    } // for

    if (0 != pflow)
        *pflow = flow;
}

} // namespace

#endif // ___OPTNET_MC_FS_3D_CXX___
