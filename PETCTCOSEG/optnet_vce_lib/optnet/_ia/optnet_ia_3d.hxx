/*
 ==========================================================================
 |   
 |   $Id: optnet_ia_3d.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___OPTNET_IA_3D_HXX___
#   define ___OPTNET_IA_3D_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_ia/optnet_ia_maxflow_3d.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class optnet_ia_3d
///  @brief Implementation of the single-surface Optimal Net algorithm
///         using the Boykov-Kolmogorov max-flow algorithm on a
///         implicit-arc graph.
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg = net_f_xy>
class optnet_ia_3d
{
    typedef optnet_ia_maxflow_3d<_Cap, net_f_xy>    graph_type;

public:

    typedef typename graph_type::roi_base_type      roi_base_type;
    typedef typename graph_type::roi_ref_type       roi_ref_type;
    typedef typename graph_type::roi_type           roi_type;

    typedef size_t                                  size_type;
    typedef _Cost                                   cost_type;
    typedef _Cap                                    capacity_type;

    typedef array_base<cost_type, _Tg>              cost_array_base_type;
    typedef array_ref<cost_type, _Tg>               cost_array_ref_type;
    typedef array<cost_type, _Tg>                   cost_array_type;
    
    typedef array_base<int>                         net_base_type;
    typedef array_ref<int>                          net_ref_type;
    typedef array<int>                              net_type;

    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    optnet_ia_3d();

    ///////////////////////////////////////////////////////////////////////
    ///  Create graph and set the cost values.
    ///
    ///  @param cost The array that contain the cost values for each voxel.
    ///
    ///  @remarks This function will re-assign the node costs of the
    ///           underlying graph based upon the given cost array.
    ///
    ///////////////////////////////////////////////////////////////////////
    void create(const cost_array_base_type& cost);

    ///////////////////////////////////////////////////////////////////////
    ///  Set the smoothness constraints.
    ///
    ///  @param smooth0 The first smoothness parameter.
    ///  @param smooth1 The second smoothness parameter.
    ///  @param circle0 Enabling/disabling circle graph construction.
    ///  @param circle1 Enabling/disabling circle graph construction.
    ///
    ///  @remarks The actual meanings of the smoothness parameters depend
    ///           on the direction setting of the net. The function will
    ///           rebuild the arcs of the underlying graph.
    ///
    ///////////////////////////////////////////////////////////////////////
    void set_params(int  smooth0 = 1,
                    int  smooth1 = 1,
                    bool circle0 = false,
                    bool circle1 = false
                    );

    ///////////////////////////////////////////////////////////////////////
    ///  Set the region of interest.
    ///
    ///  @param roi The region-of-interest mask array.
    ///
    ///////////////////////////////////////////////////////////////////////
    void set_roi(const roi_base_type& roi);

    ///////////////////////////////////////////////////////////////////////
    ///  Solve the optimal surface problem using the given cost function
    ///  and smoothness constraints.
    ///
    ///  @param net   The resulting optimal "net" surface.
    ///  @param pflow The output maximum flow value.
    ///
    ///  @remarks The pflow parameter, if not NULL, will return the
    ///           computed maximum-flow value. It is used primarily
    ///           for debugging.
    ///
    ///////////////////////////////////////////////////////////////////////
    void solve(net_base_type& net,      // [OUT]
               capacity_type* pflow = 0 // [OUT]
               );


private:

    ///////////////////////////////////////////////////////////////////////
    // Transform the costs of the graph nodes based on the given
    // cost vector.
    void transform_costs();

    ///////////////////////////////////////////////////////////////////////
    const cost_array_base_type* m_pcost;

    graph_type      m_graph;
};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_ia/optnet_ia_3d.cxx>
#   endif

#endif
