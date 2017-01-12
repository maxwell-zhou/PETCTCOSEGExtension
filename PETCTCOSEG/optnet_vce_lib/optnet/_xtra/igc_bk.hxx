/*
 ==========================================================================
 |   
 |   $Id: igc_bk.hxx 22 2005-01-16 05:10:32Z kangli $
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

      This file implements Boykov-Jolly's interactive graph cut.
      
  - Reference(s):

    [1] BOYKOV, Y., AND JOLLY, M. P.
        Interactive graph cuts for optimal boundary & region segmentation
          of objects in n-d images.
        In Proceedings of ICCV 2001.
    [2] Yin Li, Jian Sun, Chi-Keung Tang and Heung-Yeung Shum
        Lazy snapping
        ACM Transactions on Graphics (TOG)
        Volume 23, Issue 3 (August 2004)
        Special Issue: Proceedings of the 2004 SIGGRAPH Conference
        Pages: 303-308
        
 ==========================================================================
 */

#ifndef ___IGC_BK_HXX___
#   define ___IGC_BK_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/point.hxx>
#   include <optnet/_base/distance.hxx>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_xtra/bk_fs_maxflow.hxx>
#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <cmath>
#   include <limits>
#   include <map>

#   ifdef max       // The max macro may interfere with
#       undef max   //   std::numeric_limits::max().
#   endif           //

/// @namespace optnet
namespace optnet {

    /// @namespace optnet::xtra
    namespace xtra {

#   define BK_POINT_DIM 3
#   define BK_CLUSTERS  128

///////////////////////////////////////////////////////////////////////////
///  @class igc_bk
///  @brief Interactive graph cut.
///////////////////////////////////////////////////////////////////////////
template <typename _Voxel,
          typename _RealVoxel,
          typename _Cap,
          typename _Tg = net_f_xy>
class igc_bk
{
    typedef bk_fs_maxflow<_Cap, net_f_xy>   graph_type;

public:

    typedef _Voxel                          voxel_type;
    typedef _RealVoxel                      real_voxel_type;
    typedef _Cap                            capacity_type;

    typedef point<size_t, BK_POINT_DIM>     point_type;

    typedef array_base<voxel_type, _Tg>     image_base_type;
    typedef array_ref<voxel_type, _Tg>      image_ref_type;
    typedef array<voxel_type, _Tg>          image_type;
    
    typedef array_base<unsigned char>       trimap_base_type;
    typedef array_ref<unsigned char>        trimap_ref_type;
    typedef array<unsigned char>            trimap_type;

    typedef array_base<unsigned char>       mask_base_type;
    typedef array_ref<unsigned char>        mask_ref_type;
    typedef array<unsigned char>            mask_type;

    typedef std::pair<size_t, voxel_type>   index_voxel_pair;
    typedef std::map<size_t, voxel_type>    index_voxel_map;

    typedef size_t                          size_type;

    static const unsigned char              FOREGROUND;
    static const unsigned char              BACKGROUND;
    static const unsigned char              UNKNOWN;

    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    igc_bk();

    ///////////////////////////////////////////////////////////////////////
    ///  Create graph.
    ///
    ///  @param image The input image.
    ///
    ///////////////////////////////////////////////////////////////////////
    void create(const image_base_type& image);

    ///////////////////////////////////////////////////////////////////////
    ///  Add a foreground seed point.
    ///
    ///  @param i0  The first index of the seed point.
    ///  @param i1  The second index of the seed point.
    ///  @param i2  The third index of the seed point.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_foreground(size_type i0,
                               size_type i1,
                               size_type i2
                               )
    {
        assert(0 != m_pimage);
        assert(i0 < m_pimage->size_0());
        assert(i1 < m_pimage->size_1());
        assert(i2 < m_pimage->size_2());
        size_type index = m_pimage->offset(i0, i1, i2);
        // insert foreground seed
        m_fg_seeds.insert(index_voxel_pair(index, (*m_pimage)[index]));
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Clear foreground seed points.
    ///////////////////////////////////////////////////////////////////////
    inline void clear_foreground_seeds()
    {
        m_fg_seeds.clear();
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Add a background seed point.
    ///
    ///  @param i0  The first index of the seed point.
    ///  @param i1  The second index of the seed point.
    ///  @param i2  The third index of the seed point.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_background(size_type i0,
                               size_type i1,
                               size_type i2
                               )
    {
        assert(0 != m_pimage);
        assert(i0 < m_pimage->size_0());
        assert(i1 < m_pimage->size_1());
        assert(i2 < m_pimage->size_2());
        size_type index = m_pimage->offset(i0, i1, i2);
        // insert background seed
        m_bg_seeds.insert(index_voxel_pair(index, (*m_pimage)[index]));
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Clear background seed points.
    ///////////////////////////////////////////////////////////////////////
    inline void clear_background()
    {
        m_bg_seeds.clear();
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Set foreground and background seeds by setting a trimap.
    ///
    ///  @param trimap  The input trimap.
    ///
    ///  @remarks The trimap has the same size as the image. The value of
    ///           each element of the trimap is either FOREGROUND,
    ///           BACKGROUND, or UNKNOWN(0).
    ///           This function will clear the existing seeds.
    ///
    ///////////////////////////////////////////////////////////////////////
    void set_trimap(const trimap_base_type& trimap);

    ///////////////////////////////////////////////////////////////////////
    ///  Set parameters.
    ///
    ///  @param lambda  Region cost weight.
    ///  @param beta    Camera noise level.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void set_params(double lambda, double beta)
    {
        m_lambda = lambda;
        m_beta = beta;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Perform segmentation.
    ///
    ///  @param nh     The neighborhood system (6 = default, 18 or 26).
    ///  @param pflow  The output maximum flow value.
    ///
    ///  @remarks The pflow parameter, if not NULL, will return the
    ///           computed maximum-flow value. It is used primarily
    ///           for debugging.
    ///
    ///  @see     is_foreground
    ///
    ///////////////////////////////////////////////////////////////////////
    void solve(int nh = 6, capacity_type* pflow = 0);

    ///////////////////////////////////////////////////////////////////////
    ///  Perform segmentation.
    ///
    ///  @param mask   The output mask that indicates foreground and
    ///                background.
    ///  @param nh     The neighborhood system (6 = default, 18 or 26).
    ///  @param pflow  The output maximum flow value.
    ///
    ///  @remarks The pflow parameter, if not NULL, will return the
    ///           computed maximum-flow value. It is used primarily
    ///           for debugging.
    ///
    ///////////////////////////////////////////////////////////////////////
    void solve(mask_base_type& mask,     // [OUT]
               int             nh = 6,   // [IN]
               capacity_type*  pflow = 0 // [OUT]
               );

    ///////////////////////////////////////////////////////////////////////
    ///  Determines if the given voxel is a foreground voxel.
    ///
    ///  @param  i0   The first  index of the node. 
    ///  @param  i1   The second index of the node. 
    ///  @param  i2   The third  index of the node. 
    ///
    ///  @return Returns true if the given voxel belongs to the foreground.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline bool is_foreground(size_type i0,
                              size_type i1,
                              size_type i2
                              ) const
    {
        return m_graph.in_source_set(i0, i1, i2);
    }


private:

    // Cluster foreground and background seeds respectively to
    //   BK_CLUSTERS bins using K-means.
    void    cluster_seeds();
    int     kmeans(const index_voxel_map& seeds,
                   real_voxel_type clusters[]
                   );

    ///////////////////////////////////////////////////////////////////////
    // Compute the minimum distance between a voxel and the foreground
    //   clusters.
    double distance_to_foreground(const voxel_type& v)
    {
        double min_dist = -1;
        if (m_fg_num_clusters > 0) {
            // Compute minimum distance.
            min_dist = distance_euclidean2(m_fg_clusters[0],
                (real_voxel_type)(v));
            for (int i = 1; i < m_fg_num_clusters; ++i) {
                double d = distance_euclidean2(m_fg_clusters[i], 
                    (real_voxel_type)(v));
                if (d < min_dist) min_dist = d;
            }
        }
        return min_dist;
    }

    ///////////////////////////////////////////////////////////////////////
    // Compute the minimum distance between a voxel and the background
    //   clusters.
    double distance_to_background(const voxel_type& v)
    {
        double min_dist = -1;
        if (m_bg_num_clusters > 0) {
            // Compute minimum distance.
            min_dist = distance_euclidean2(m_bg_clusters[0], 
                (real_voxel_type)(v));
            for (int i = 1; i < m_bg_num_clusters; ++i) {
                double d = distance_euclidean2(m_bg_clusters[i], 
                    (real_voxel_type)(v));
                if (d < min_dist) min_dist = d;
            }
        }
        return min_dist;
    }

    ///////////////////////////////////////////////////////////////////////
    inline void energy_likelihood(const voxel_type& v,
                                  capacity_type&    fg_energy,
                                  capacity_type&    bg_energy
                                  )
    {
        double dist_fg = distance_to_foreground(v);
        double dist_bg = distance_to_background(v);

        assert(dist_fg >= 0);
        assert(dist_bg >= 0);

        double sumdist = dist_fg + dist_bg;
        
        fg_energy = (capacity_type)(dist_fg / sumdist);
        bg_energy = (capacity_type)(dist_bg / sumdist);
    }

    ///////////////////////////////////////////////////////////////////////
    inline capacity_type energy_prior(size_type tail_i0,
                                      size_type tail_i1,
                                      size_type tail_i2,
                                      size_type head_i0,
                                      size_type head_i1,
                                      size_type head_i2
                                      )
    {
        const voxel_type& v1 = (*m_pimage)(tail_i0, tail_i1, tail_i2);
        const voxel_type& v2 = (*m_pimage)(head_i0, head_i1, head_i2);
        double d = (tail_i0 - head_i0) * (tail_i0 - head_i0) +
                   (tail_i1 - head_i1) * (tail_i1 - head_i1) +
                   (tail_i2 - head_i2) * (tail_i2 - head_i2);
        double c = distance_euclidean2(v1, v2);
        c = m_lambda / sqrt(d) * exp(-c * m_beta);
        return (capacity_type)(c);
    }

    // Build graph arcs.
    void    build_arcs(int nh);
    void    build_t_link(size_type i0, size_type i1, size_type i2);
    void    build_n_link(size_type tail_i0,
                         size_type tail_i1,
                         size_type tail_i2,
                         size_type head_i0,
                         size_type head_i1,
                         size_type head_i2
                         );

    ///////////////////////////////////////////////////////////////////////
    // member variables
    const   image_base_type* m_pimage;

    double          m_lambda, m_beta;
    index_voxel_map m_fg_seeds;     // foreground seeds
    index_voxel_map m_bg_seeds;     // background seeds
    real_voxel_type m_fg_clusters[BK_CLUSTERS];
    real_voxel_type m_bg_clusters[BK_CLUSTERS];
    int             m_fg_num_clusters;
    int             m_bg_num_clusters;
    capacity_type   m_inf;
    graph_type      m_graph;
};

    } // namespace
} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_xtra/igc_bk.cxx>
#   endif

#endif // ___IGC_BK_HXX___
