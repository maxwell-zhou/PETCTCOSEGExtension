/*
 ==========================================================================
 |   
 |   $Id: optnet_mc_fs_3d_double_a_ps.hxx 136 2005-02-06 23:00:47Z kangli $
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

#ifndef ___OPTNET_MC_FS_3D_DOUBLE_A_PS_HXX___
#   define ___OPTNET_MC_FS_3D_DOUBLE_A_PS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4127)
#   endif

#   include <optnet/_alpha/graphspec.hxx>
#   include <optnet/_alpha/optnet_mc_fs_maxflow.hxx>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_base/array2.hxx>
#   include <optnet/_base/point3.hxx>
#   include <vector>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class optnet_mc_fs_3d_double_a_ps
///  @brief Implementation of the double-surface Optimal Net algorithm
///         using the Boykov-Kolmogorov max-flow algorithm on a
///         forward-star represented multicolumn graph.
///////////////////////////////////////////////////////////////////////////
template <typename _Voxel, typename _Cap, typename _Tg = net_f_xy>
class optnet_mc_fs_3d_double_a_ps
{
    typedef optnet_mc_fs_maxflow<_Cap>  graph_type;

public:

    typedef _Voxel                          voxel_type;         /// voxel type.
    typedef _Cap                            capacity_type;      /// Graph arc capacity type.
    typedef double                          real_value_type;    /// Real value type (default: double).
    typedef size_t                          size_type;          /// Size type.

    typedef graphspec                       graphspec_type;     /// Graph specification type.

    typedef array_base<voxel_type, _Tg>     image_base_type;    /// Image base class.
    typedef array_ref<voxel_type, _Tg>      image_ref_type;     /// Image reference class.
    typedef array<voxel_type, _Tg>          image_type;         /// Image class.
    
    typedef point3<real_value_type>         point_type;         /// Point type with real-valued coordinates.
    typedef std::vector<real_value_type>    real_value_vector;  /// Real-valued vector.


    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    optnet_mc_fs_3d_double_a_ps();


    //
    //TODO:
    //

    ///////////////////////////////////////////////////////////////////////
    ///  Solve the optimal surfaces problem using the given cost function
    ///  and smoothness constraints.
    ///////////////////////////////////////////////////////////////////////
    void solve(capacity_type* pflow = NULL);


private:
    int                     m_sc;
    graph_type              m_graph;
    real_value_vector       m_vfact;
    graphspec_type*         m_pgs;
};
    
} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_alpha/optnet_mc_fs_3d_double_a_ps.cxx>
#   endif

#endif // ___OPTNET_MC_FS_3D_DOUBLE_A_PS_HXX___
