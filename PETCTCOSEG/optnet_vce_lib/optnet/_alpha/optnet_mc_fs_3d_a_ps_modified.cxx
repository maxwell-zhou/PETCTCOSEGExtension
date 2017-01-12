/*
 ==========================================================================
 |   
 |   $Id: optnet_mc_fs_3d_a_ps.cxx 141 2005-02-07 05:59:03Z kangli $
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

#ifndef ___OPTNET_MC_FS_3D_A_PS_CXX___
#   define ___OPTNET_MC_FS_3D_A_PS_CXX___

#   include <optnet/_alpha/optnet_mc_fs_3d_a_ps_modified.hxx>
#   include <optnet/_utils/interp_bspline3.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _image, typename _Cap, typename _Tg>
optnet_mc_fs_3d_a_ps<_image, _Cap, _Tg>::optnet_mc_fs_3d_a_ps() :
    m_sc(1), m_pgs(0)
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _image, typename _Cap, typename _Tg>
void optnet_mc_fs_3d_a_ps<_image, _Cap, _Tg>::create(
                                        const image_base_type& image,
                                        graphspec_type&        spec,
                                        real_value_type        rmin,
                                        real_value_type        rmax,
                                        size_type              npc,
                                        int                    sc,
                                        bool                   dark
                                        )
{
    using namespace optnet::utils;

    // Interpolator type.
    typedef interp_bspline3<voxel_type,
                            real_value_type,
                            _Tg>            interp_type;

    size_type s0, s1, s2;
    int       i, j, numcols, numadjs;

    assert(npc >= 2); // at least two nodes per column
    assert(sc  >= 0);

    s0 = image.size_0();
    s1 = image.size_1();
    s2 = image.size_2();

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

    // m_vfact[j] * col[i].normal + col[i].to is the true physical
    // position of node j of column i in the volume.
    m_vfact.resize(npc);

    real_value_type frac =  (rmax - rmin) / (real_value_type)(npc - 1);
    for (j = 0; j < (int)npc; ++j)
        m_vfact[j] = j * frac + rmin;

    // Initialize interpolator.
    interp_type interp3(image, 3); // cubic b-spline interpolation

    // Loop through all columns.
    for (i = 0; i < numcols; ++i) {
        typename graphspec_type::column_const_reference
                    col = spec.get_column(i);
        //real_value_type cv0, cv1, d10, d11, pos[3];
		real_value_type cv0, cv1, cv2, origin0, origin1, origin2, pos[3];
        capacity_type cap;
		bool flag = false;
		float threshold = 10;

        pos[0] = real_value_type(
            col.pos[0] + col.normal[0] * m_vfact[0]);
        pos[1] = real_value_type(
            col.pos[1] + col.normal[1] * m_vfact[0]);
        pos[2] = real_value_type(
            col.pos[2] + col.normal[2] * m_vfact[0]);
                
        cv0 = interp3.interp(pos[0], pos[1], pos[2]);
		origin0 = cv0;
		cv0 = cv0 + 10;
        //d10 = (real_value_type)(std::numeric_limits<real_value_type>::max() / 2.0);

        // Add t-link with capacity 1 to the source.
        m_graph.add_st_arc(1, 0, 0, i);
		//
		if ( npc > 1 )
		{
			pos[0] = real_value_type(col.pos[0] + col.normal[0] * m_vfact[1]);
            pos[1] = real_value_type(col.pos[1] + col.normal[1] * m_vfact[1]);
            pos[2] = real_value_type(col.pos[2] + col.normal[2] * m_vfact[1]);
            cv1 = interp3.interp(pos[0], pos[1], pos[2]);
			origin1 = cv1;
			cv1 = cv1 + 10;
		}


        // Loop through all nodes in the current column.
        for (j = 1; j < (int)npc; ++j) {
            
            //d11 = dark ? cv0 - cv1 : cv1 - cv0; // first derivative
           
            if ( j < (int)npc - 1 )
			{
				pos[0] = real_value_type( col.pos[0] + col.normal[0] * m_vfact[j + 1] );
                pos[1] = real_value_type( col.pos[1] + col.normal[1] * m_vfact[j + 1] );
                pos[2] = real_value_type( col.pos[2] + col.normal[2] * m_vfact[j + 1]);
			}

			
			cv2 = interp3.interp(pos[0], pos[1], pos[2]);
			origin2 = cv2;
			cv2 = cv2 + 1.0 / j * 10;

			
			if ( origin2 - origin1 > threshold && origin1 - origin0 > threshold )
				cv1 = cv1 + 100;
			else if ( origin2 - origin1 < -threshold && origin1 - origin0 < -threshold && flag == false)
			{
				 cv1 = (( cv1 - 100 ) > 0 ? ( cv1 - 100 ) : 0 );
				 flag = true;
			} 
			


            cap = (capacity_type)(cv1 - cv0);
            
            if (cap >= 0) // non-negative -> connect to t
                m_graph.add_st_arc(0, +cap, j, i);
            else          // negative     -> connect to s
                m_graph.add_st_arc(-cap, 0, j, i);

			cv0 = cv1;
			origin0 = origin1;
			cv1 = cv2;
			origin1 = origin2;

            //cv0 = cv1;
            //d10 = d11;
        } // for j
	} // for i
	
    // Save the boundary and origin.
    spec.set_boundary(0.0, (double)(s0 - 1), 
                      0.0, (double)(s1 - 1),
                      0.0, (double)(s2 - 1));

    // Save the parameters.
    m_sc = sc;
    m_pgs = &spec;
}

///////////////////////////////////////////////////////////////////////////
template <typename _image, typename _Cap, typename _Tg>
void optnet_mc_fs_3d_a_ps<_image, _Cap, _Tg>::solve(capacity_type* pflow)
{
    capacity_type   flow;
    size_type       i, j, npc, numcols;

    if (NULL == m_pgs) {
        throw_exception(std::runtime_error(
            "optnet_mc_fs_3d_a_ps::solve: The graph was not created properly."
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

#endif // ___OPTNET_MC_FS_3D_A_PS_CXX___
