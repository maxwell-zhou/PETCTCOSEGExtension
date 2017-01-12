/*
 ==========================================================================
 |   
 |   $Id: optnet_mc_fs_3d_double_ps.cxx 170 2005-02-14 05:55:28Z kangli $
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

#ifndef ___OPTNET_MC_FS_3D_DOUBLE_PS_CXX___
#   define ___OPTNET_MC_FS_3D_DOUBLE_PS_CXX___

#   include <optnet/_alpha/optnet_mc_fs_3d_double_ps.hxx>
#   include <optnet/_utils/interp_bspline3.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
optnet_mc_fs_3d_double_ps<_Cost, _Cap, _Tg>::
optnet_mc_fs_3d_double_ps() :
    m_pgs(0), m_dmin(0), m_dmax(0)
{
    m_asc[0] = 0;
    m_asc[1] = 0;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void optnet_mc_fs_3d_double_ps<_Cost, _Cap, _Tg>::
create(graphspec_type& spec,
       real_value_type rmin,
       real_value_type rmax,
       size_type       npc
       )
{
    assert(npc >= 2);   // At least two nodes per column.

    m_pgs = &spec;      // Save a pointer to the graphspec.
    
    int numcols = (int)spec.num_columns();

    // ====================================================================
    // Create a multicolumn graph.
    m_graph.create(npc, numcols * 2);

    // ====================================================================
    // col[i].pos + col[i].normal m_vfact[j] is the true physical position
    // of node j of column i in the volume.
    m_vfact.resize(npc);

    real_value_type frac =  (rmax - rmin) / (real_value_type)(npc - 1);
    for (int j = 0; j < (int)npc; ++j)
        m_vfact[j] = j * frac + rmin;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void optnet_mc_fs_3d_double_ps<_Cost, _Cap, _Tg>::
set_cost_function(size_type index,
                  const     cost_array_base_type& cost,
                  int       smoothness
                  )
{
    using namespace optnet::utils;

    assert(index == 0 || index == 1);
    assert(smoothness >= 0);
    assert(NULL != m_pgs);

    int i, j;

    m_asc[index] = smoothness;

    // Interpolator type.
    typedef interp_bspline3<cost_type,
                            real_value_type,
                            _Tg>            interp_type;

    int npc        = (int)m_graph.nodes_per_column();
    int numcols    = (int)m_pgs->num_columns();
    int col_offset = (int)index * numcols; // column offset

    // ====================================================================
    // Create t-links.

    // Initialize interpolator.
    interp_type interp3(cost, 1); // trilinear interpolation

    // Loop through all columns.
    for (i = 0; i < numcols; ++i) {

        typename graphspec_type::column_const_reference
                    col = m_pgs->get_column(i);

        real_value_type cv0, cv1, pos[3];
        capacity_type cap;

        pos[0] = real_value_type(col.normal[0] * m_vfact[0] + col.pos[0]);
        pos[1] = real_value_type(col.normal[1] * m_vfact[0] + col.pos[1]);
        pos[2] = real_value_type(col.normal[2] * m_vfact[0] + col.pos[2]);
                
        cv0 = interp3.interp(pos[0], pos[1], pos[2]);

        // Add t-link with capacity 1 to the source.
        m_graph.add_st_arc(1, 0, 0, i + col_offset);

        // Loop through all nodes in the current column.
        for (j = 1; j < (int)npc; ++j) {
            pos[0] = real_value_type
                        (col.normal[0] * m_vfact[j] + col.pos[0]);
            pos[1] = real_value_type
                        (col.normal[1] * m_vfact[j] + col.pos[1]);
            pos[2] = real_value_type
                        (col.normal[2] * m_vfact[j] + col.pos[2]);

            cv1 = interp3.interp(pos[0], pos[1], pos[2]);
            cap = (capacity_type)(cv1 - cv0);
            
            if (cap >= 0) // non-negative -> connect to t
                m_graph.add_st_arc(0, +cap, j, i + col_offset);
            else          // negative     -> connect to s
                m_graph.add_st_arc(-cap, 0, j, i + col_offset);

            cv0 = cv1;
        } // for j
    } // for i

    // Save image boundary.
    size_type s0 = cost.size_0();
    size_type s1 = cost.size_1();
    size_type s2 = cost.size_2();

    m_pgs->set_boundary(0.0, (double)(s0 - 1), 
                        0.0, (double)(s1 - 1),
                        0.0, (double)(s2 - 1));
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void optnet_mc_fs_3d_double_ps<_Cost, _Cap, _Tg>::
set_relation(int dmin, int dmax)
{
    assert(NULL != m_pgs);
    assert(dmin >= 0 && dmax >= 0);

    int i, j, ii, jj;

    // Surface 0 is assumed to be the inner surface, when the surface
    // normals are pointing outwards.

    int dist01 = dmin, dist10 = -dmax;

    // For surface 1, the nodes on plane dmin are strongly connected.
    // The nodes below that are ignored.
    // For surface 0, the nodes on 0-th plane are strongly connected.
    // The nodes on or above plane npc-dmin are ignored.

    int numcols = (int)m_pgs->num_columns();
    int numadjs = (int)m_pgs->num_adjacencies();
    int npc     = (int)m_graph.nodes_per_column();

    // ====================================================================
    // Construct intra-column arcs.
    for (i = 0; i < numcols; ++i) {

        // -- surface 0:
        for (j = 1; j < npc - dmin; ++j) {
            // Add arcs with infinite capacities.
            m_graph.add_arc(j, i, j - 1, i);
        } // for j

        // -- surface 1:
        for (j = dmin + 1; j < npc; ++j) {
            // Add arcs with infinite capacities.
            m_graph.add_arc(j, i + numcols, j - 1, i + numcols);
        } // for j

    } // for i

    // ====================================================================
    // Construct inter-column arcs.
    for (i = 0; i < numadjs; ++i) {

        typename graphspec_type::adjacency_const_reference
                    adj = m_pgs->get_adjacency(i);

        int c00 = adj.id[0];
        int c01 = adj.id[0] + numcols;
        int c10 = adj.id[1];
        int c11 = adj.id[1] + numcols;

        // -- surface 0:
        //    Add arcs with infinite capacities.
        for (j = m_asc[0] + 1; j < npc - dmin; ++j) {
            m_graph.add_arc(j, c00, j - m_asc[0], c10);
            m_graph.add_arc(j, c10, j - m_asc[0], c00);
        } // for j

        // -- surface 1:
        //    Add arcs with infinite capacities.
        for (j = m_asc[1] + dmin + 1; j < npc; ++j) {
            m_graph.add_arc(j, c01, j - m_asc[1], c11);
            m_graph.add_arc(j, c11, j - m_asc[1], c01);
        } // for j

        // base sets:
        m_graph.add_arc(   0, c00,    0, c10);
        m_graph.add_arc(   0, c10,    0, c00);
        m_graph.add_arc(dmin, c01, dmin, c11);
        m_graph.add_arc(dmin, c11, dmin, c01);

    } // for i

    // ====================================================================
    // Construct inter-surface arcs.
    for (i = 0; i < numcols; ++i) {
        
        ii = i + numcols;

        // -- surface 0 -> surface 1:
        for (j = 1; j < npc - dmin; ++j) {
            // Add arcs with infinite capacities.
            jj = j + dist01;
            if (jj <= dmin || jj >= npc) continue;
            m_graph.add_arc(j, i, j + dist01, ii);
        } // for j

        // -- surface 1 -> surface 0:
        for (j = dmin + 1; j < npc; ++j) {
            // Add arcs with infinite capacities.
            jj = j + dist10;
            if (jj <= 0 || jj >= npc - dmin) continue;
            m_graph.add_arc(j, ii, j + dist10, i);
        } // for j

    } // for i

    // Save parameters.
    m_dmin = dmin;
    m_dmax = dmax;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void optnet_mc_fs_3d_double_ps<_Cost, _Cap, _Tg>::
solve(capacity_type* pflow)
{
    capacity_type   flow;
    int             i, j, npc, numcols;

    if (NULL == m_pgs) {
        throw_exception(std::runtime_error(
            "optnet_mc_fs_3d_double_ps::solve: The graph was not created properly."
            ));
    }

    numcols = (int)m_pgs->num_columns();
    npc     = (int)m_graph.nodes_per_column();

    assert(numcols != 0);
    assert(npc != 0);

    // Calculate max-flow/min-cut.
    flow = m_graph.solve();
    
    for (i = 0; i < numcols; ++i) {
        int lo, hi, ii;

        // -- surface 0:
        //    Find upper envelope for each column.
        lo = 0;
        hi = npc - m_dmin;
        for (j = lo; j < hi; ++j)
            if (!m_graph.in_source_set(j, i)) break;

        if (j != lo)
            m_pgs->set_fact(i, 0, m_vfact[j - 1]);
        else
            m_pgs->set_fact(i, 0, -1);

        // -- surface 1:
        ii =  i + numcols;
        //    Find upper envelope for each column.
        lo = m_dmin;
        hi = npc;
        for (j = lo; j < hi; ++j)
            if (!m_graph.in_source_set(j, ii)) break;

        if (j != lo)
            m_pgs->set_fact(i, 1, m_vfact[j - 1]);
        else
            m_pgs->set_fact(i, 1, -1);
    } // for

    m_pgs->positions_from_facts();
    m_pgs->set_column_type(6);

    if (0 != pflow)
        *pflow = flow;
}

} // namespace

#endif // ___OPTNET_MC_FS_3D_DOUBLE_PS_CXX___
