/*
 ==========================================================================
 |   
 |   $Id: igc_bk.cxx 80 2005-01-22 07:49:12Z Administrator $
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

#ifndef ___IGC_BK_CXX___
#   define ___IGC_BK_CXX___

#   include <cstdio>
#   include <cstdlib>
#   include <optnet/_base/except.hxx>
#   include <optnet/_base/secure_s.hxx>
#   include <optnet/_xtra/igc_bk.hxx>

namespace optnet { namespace xtra {

///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::igc_bk() :
    m_pimage(0), m_lambda(50), m_beta(0.5),
    m_fg_num_clusters(0), m_bg_num_clusters(0),
    m_inf(std::numeric_limits<capacity_type>::max())
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
void igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::create(
    const image_base_type& image
    )
{
    if (!m_graph.create(image.size_0(), 
                        image.size_1(), 
                        image.size_2()
                        )) {
        throw_exception(std::runtime_error(
            "igc_bk::create: Could not create graph."
            ));
    }
    // Save a pointer to the input image.
    m_pimage = &image;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
void igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::set_trimap(
    const trimap_base_type& trimap
    )
{
    assert(0 != m_pimage);

    if (trimap.size_0() != m_pimage->size_0() || 
        trimap.size_1() != m_pimage->size_1() ||
        trimap.size_2() != m_pimage->size_2()
        ) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "igc_bk::set_trimap: The trimap size must match the image size."
        ));
    }

    // Clear all seed points.
    m_fg_seeds.clear();
    m_bg_seeds.clear();

    for (size_type i = 0; i < trimap.size(); ++i) {
        if (trimap[i] == FOREGROUND) { // foreground voxel
            m_fg_seeds.insert(index_voxel_pair(i, (*m_pimage)[i]));
        } else if (trimap[i] == BACKGROUND) { // background voxel
            m_bg_seeds.insert(index_voxel_pair(i, (*m_pimage)[i]));
        }
    } // for
}

///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
void igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::solve(
    int            nh,
    capacity_type* pflow
    )
{
    if (nh != 6 && nh != 18 && nh != 26) {
        throw_exception(std::invalid_argument(
            "igc_bk::solve: Unsupported neighborhood system"
        ));
    }
    
    assert(0 != m_pimage);

    capacity_type flow;

    // Do k-means clustering of the seeds.
    cluster_seeds();

    // Construct graph arcs.
    build_arcs(nh);

    // Calculate max-flow/min-cut.
    // Once this returns, user can call is_foreground to determine if a
    //   voxel is foreground or background.
    flow = m_graph.solve();

    // Return the maximum-flow value (optional).
    if (0 != pflow)
        *pflow = flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
void igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::solve(
    mask_base_type& mask,
    int             nh,
    capacity_type*  pflow
    )
{
    size_type i0, i1, i2;

    if (mask.size_0() != m_graph.size_0() || 
        mask.size_1() != m_graph.size_1() ||
        mask.size_2() != m_graph.size_2()
        ) {
        // Throw an invalid_argument exception.
        char errmsg[1024];
        secure_sprintf(
            errmsg, sizeof(errmsg),
            "igc_bk::solve: The mask size must match the image size (%d/%d/%d).",
            m_graph.size_0(),
            m_graph.size_1(),
            m_graph.size_2()
        );
        throw_exception(std::invalid_argument(errmsg));
    }

    // Perform segmentation.
    solve(nh, pflow);

    // Generate image mask that indicates foreground and background.
    for (i2 = 0; i2 < m_graph.size_2(); ++i2) {
        for (i1 = 0; i1 < m_graph.size_1(); ++i1) {
            for (i0 = 0; i0 < m_graph.size_0(); ++i0) {

                if (m_graph.in_source_set(i0, i1, i2)) {
                    // foreground voxel
                    mask(i0, i1, i2) = FOREGROUND;
                }
                else {
                    // background voxel
                    mask(i0, i1, i2) = BACKGROUND;
                }

            } // for i0
        } // for i1
    } // for i2

}

///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
void igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::cluster_seeds()
{
    // Cluster foreground seeds.
    m_fg_num_clusters = kmeans(m_fg_seeds, m_fg_clusters);
    
    // Cluster background seeds.
    m_bg_num_clusters = kmeans(m_bg_seeds, m_bg_clusters);

}

///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
int igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::kmeans(
    const index_voxel_map&  seeds,
    real_voxel_type         clusters[BK_CLUSTERS]
    )
{
    //
    // NOTE:
    //   This is a most naive implementation of k-means,
    //   ...may be replaced in the future.
    //
    static const double EPS = 1e-6;

    int             num_clusters = 0;
    typename        index_voxel_map::const_iterator it;
    real_voxel_type clusters2[BK_CLUSTERS];

    // Do k-means clustering.
    int i, j, k;

    // Initialize: randomly select initial centroids.
    srand(2004); // fixed seed
    for (it = seeds.begin(); it != seeds.end(); ++it) {
        bool used = false;
        int  use  = rand() % 2; // binomial
        if (0 != use) {
            for (i = 0; i < num_clusters; ++i) {
                if (clusters[i] == it->second) {
                    used = true;
                    break;
                }
            }
            if (!used) {
                clusters[num_clusters] = it->second;
                if (++num_clusters == BK_CLUSTERS)
                    break;
            }
        }
    } // for

    // Begin clustering
    if (num_clusters >= BK_CLUSTERS) {

        // If the number of clusters is less than expected, do nothing.
        // Else...
        
        real_voxel_type sums[BK_CLUSTERS];
        size_type       nsum[BK_CLUSTERS];

        int    min_i = -1;
        double min_dist = std::numeric_limits<double>::max();
        bool   diff;

        do {
            diff = false;

            for (i = 0; i < num_clusters; ++i) {
                clusters2[i] = clusters[i];
                sums[i]      = real_voxel_type();
                nsum[i]      = 0;
            }

            // Loop over all seeds.
            for (it = seeds.begin(); it != seeds.end(); ++it) {

                for (i = 0; i < num_clusters; ++i) {
                    double dist = distance_euclidean2(clusters[i], (real_voxel_type)(it->second));
                    if (dist < min_dist) {
                        min_dist = dist;
                        min_i    = i;
                    }
                } // for (i...
                
                assert(min_i < num_clusters);
                sums[min_i] += it->second;
                nsum[min_i] += 1;
            } // for (it...

            j = 0;
            k = num_clusters;
            for (i = 0; i < num_clusters; ++i) {
                if (0 != nsum[i]) {
                    clusters[j++] = sums[i] / nsum[i];
                }
            } // for
            
            if (j == num_clusters) {
                for (i = 0; i < num_clusters; ++i) {
                    if (distance_manhattan(clusters[i], clusters2[i]) > EPS) {
                        diff = true;
                        break;
                    } // if
                } // for
            } else {
                num_clusters = j;
                diff = true;
            }

        } while (diff);

    } // if

    assert(num_clusters >= 1);

    return num_clusters;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
void igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::build_arcs(int nh)
{

#   define BUILD_N_LINKS_6  \
    build_n_link(i0, i1, i2, i0 - 1, i1,     i2    ); \
    build_n_link(i0, i1, i2, i0 + 1, i1,     i2    ); \
    build_n_link(i0, i1, i2, i0,     i1 - 1, i2    ); \
    build_n_link(i0, i1, i2, i0,     i1 + 1, i2    ); \
    build_n_link(i0, i1, i2, i0,     i1,     i2 - 1); \
    build_n_link(i0, i1, i2, i0,     i1,     i2 + 1)

#   define BUILD_N_LINKS_18 BUILD_N_LINKS_6;  \
    build_n_link(i0, i1, i2, i0 - 1, i1 - 1, i2    ); \
    build_n_link(i0, i1, i2, i0 + 1, i1 - 1, i2    ); \
    build_n_link(i0, i1, i2, i0 - 1, i1 + 1, i2    ); \
    build_n_link(i0, i1, i2, i0 + 1, i1 + 1, i2    ); \
    build_n_link(i0, i1, i2, i0 - 1, i1,     i2 - 1); \
    build_n_link(i0, i1, i2, i0 + 1, i1,     i2 - 1); \
    build_n_link(i0, i1, i2, i0 - 1, i1,     i2 + 1); \
    build_n_link(i0, i1, i2, i0 + 1, i1,     i2 + 1); \
    build_n_link(i0, i1, i2, i0,     i1 - 1, i2 - 1); \
    build_n_link(i0, i1, i2, i0,     i1 + 1, i2 - 1); \
    build_n_link(i0, i1, i2, i0,     i1 - 1, i2 + 1); \
    build_n_link(i0, i1, i2, i0,     i1 + 1, i2 + 1)

#   define BUILD_N_LINKS_26 BUILD_N_LINKS_18; \
    build_n_link(i0, i1, i2, i0 - 1, i1 - 1, i2 - 1); \
    build_n_link(i0, i1, i2, i0 + 1, i1 - 1, i2 - 1); \
    build_n_link(i0, i1, i2, i0 - 1, i1 + 1, i2 - 1); \
    build_n_link(i0, i1, i2, i0 + 1, i1 + 1, i2 - 1); \
    build_n_link(i0, i1, i2, i0 - 1, i1 - 1, i2 + 1); \
    build_n_link(i0, i1, i2, i0 + 1, i1 - 1, i2 + 1); \
    build_n_link(i0, i1, i2, i0 - 1, i1 + 1, i2 + 1); \
    build_n_link(i0, i1, i2, i0 + 1, i1 + 1, i2 + 1)

    size_type i0, i1, i2;

    // Clear all constructed arcs.
    m_graph.clear_arcs();

    //
    // Begin building arcs and assigning arc weights.
    //
    switch (nh) {
    case 6:
        for (i2 = 0; i2 < m_graph.size_2(); ++i2) {
            for (i1 = 0; i1 < m_graph.size_1(); ++i1) {
                for (i0 = 0; i0 < m_graph.size_0(); ++i0) {
                    build_t_link(i0, i1, i2);   // Build t-link.
                    BUILD_N_LINKS_6;            // Build n-link ( 6-neighbor)
                } // for i0
            } // for i1
        } // for i2
        break;
    case 18:
        for (i2 = 0; i2 < m_graph.size_2(); ++i2) {
            for (i1 = 0; i1 < m_graph.size_1(); ++i1) {
                for (i0 = 0; i0 < m_graph.size_0(); ++i0) {
                    build_t_link(i0, i1, i2);   // Build t-link.
                    BUILD_N_LINKS_18;           // Build n-link (18-neighbor)
                } // for i0
            } // for i1
        } // for i2
        break;
    case 26:
        for (i2 = 0; i2 < m_graph.size_2(); ++i2) {
            for (i1 = 0; i1 < m_graph.size_1(); ++i1) {
                for (i0 = 0; i0 < m_graph.size_0(); ++i0) {
                    build_t_link(i0, i1, i2);   // Build t-link.
                    BUILD_N_LINKS_26;           // Build n-link (26-neighbor)
                } // for i0
            } // for i1
        } // for i2
        break;
    default:
        ; // do nothing
    } // switch
}


///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
void igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::build_t_link(
    size_type i0, size_type i1, size_type i2
    )
{
    size_type index = m_pimage->offset(i0, i1, i2);
    
    // Build t-links.
    //
    typename index_voxel_map::iterator it_fg = m_fg_seeds.find(index);
    typename index_voxel_map::iterator it_bg = m_bg_seeds.find(index);

    if (it_fg != m_fg_seeds.end()) {

#   ifndef __BK_INTERACTIVE_NO_CHECK_SEEDS__
        if (it_bg != m_bg_seeds.end()) {
            // The voxel is marked as both foreground
            // and background. Throw an exception.
            throw_exception(std::runtime_error(
                "igc_bk::build_arcs: Ambiguous seed."
            ));
        }
#   endif
        // foreground seeds
        m_graph.add_st_arc(m_inf, (capacity_type)0, i0, i1, i2);
    }
    else if (it_bg != m_bg_seeds.end()) {
        // background seeds
        m_graph.add_st_arc((capacity_type)0, m_inf, i0, i1, i2);
    }
    else {
        // other voxels
        
        capacity_type fg_energy, bg_energy;

        // Compute likelihood energy.
        energy_likelihood(
            (*m_pimage)(i0, i1, i2),
            fg_energy,
            bg_energy
        );
        
        m_graph.add_st_arc(
            bg_energy,  // source
            fg_energy,  // sink
            i0, i1, i2
        );

    } // if
}

///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
void igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::build_n_link(
    size_type tail_i0,
    size_type tail_i1,
    size_type tail_i2,
    size_type head_i0,
    size_type head_i1,
    size_type head_i2
    )
{
    assert(0 != m_pimage);

    if ((tail_i0 >= m_graph.size_0()) ||
        (tail_i1 >= m_graph.size_1()) ||
        (tail_i2 >= m_graph.size_2()) ||
        (head_i0 >= m_graph.size_0()) ||
        (head_i1 >= m_graph.size_1()) ||
        (head_i2 >= m_graph.size_2())
        )
    {
        // Node out of boundary.
        return;
    }
    
    // Compute prior energy.
    capacity_type energy = energy_prior(
        tail_i0, tail_i1, tail_i2,
        head_i0, head_i1, head_i2
    );

    m_graph.add_arc(
        tail_i0, tail_i1, tail_i2,  // tail node
        head_i0, head_i1, head_i2,  // head node
        energy,                     // forward capacity
        0                           // reverse capacity
    );
}


///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
    const unsigned char
        igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::FOREGROUND = 0xFF;

template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
    const unsigned char
        igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::BACKGROUND = 0x00;

template <typename _Voxel, typename _RealVoxel, typename _Cap, typename _Tg>
    const unsigned char
        igc_bk<_Voxel, _RealVoxel, _Cap, _Tg>::UNKNOWN    = 0x7F;

}} // namespace

#endif // ___IGC_BK_CXX___
