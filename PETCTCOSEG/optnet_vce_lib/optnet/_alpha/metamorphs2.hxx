/*
 ==========================================================================
 |   
 |   $Id: metamorphs2.hxx 186 2005-02-27 02:04:22Z kangli $
 |
 |   Written by Kang Li <kangl@cmu.edu>
 |   Department of Electrical and Computer Engineering
 |   Carnegie Mellon University
 |   
 ==========================================================================
 |   This file is m_para_a part of the OptimalNet library.
 ==========================================================================
 | Copyright (c) 2004-2005 Kang Li <kangl@cmu.edu>. All Rights Reserved.
 | 
 | This software is supplied under the terms of m_para_a license agreement or
 | nondisclosure agreement  with the author  and may not be copied  or
 | disclosed except in accordance with the terms of that agreement.
 ==========================================================================
 */

/*
 ==========================================================================
  - Reference(s):

    [1] Xiaolei Huang, Dimitris Metaxas, Ting Chen
        MetaMorphs: Deformable Shape and Texture Models
        Proc. CVPR 2004
 ==========================================================================
 */

#ifndef ___METAMORPHS2_HXX___
#   define ___METAMORPHS2_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>
#   include <optnet/_alpha/contour.hxx>
#   include <optnet/_alpha/cubic_spline_ffd2.hxx>
#   include <optnet/_utils/interp_bspline2.hxx>
#   include <optnet/_alpha/fast_marching2.hxx>
#   include <optnet/_base/secure_s.hxx>
#   include <optnet/_base/io/png.hxx>
#   include <optnet/_base/debug.hxx>
#   include <limits>
#   include <string>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class metamorphs2
///  @brief 2-D MetaMorphs deformable shape and texture model. 
///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Real = double>
class metamorphs2
{
public:
    
    typedef _Ty                             pixel_type;
    typedef pixel_type&                     reference;
    typedef const pixel_type&               const_reference;

    typedef array2_base<pixel_type>         image_base_type;
    typedef array2_ref<pixel_type>          image_ref_type;
    typedef array2<pixel_type>              image_type;

    typedef unsigned char                   mask_value_type;

    typedef array2_base<mask_value_type>    mask_image_base_type;
    typedef array2_ref<mask_value_type>     mask_image_ref_type;
    typedef array2<mask_value_type>         mask_image_type;

    typedef _Real                           real_value_type;

    typedef array2_base<real_value_type>    real_image_base_type;
    typedef array2_ref<real_value_type>     real_image_ref_type;
    typedef array2<real_value_type>         real_image_type;

    typedef size_t                          size_type;


    ///////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    ///////////////////////////////////////////////////////////////////////
    metamorphs2() :
        m_verbose(0), m_debug_image_scale(2),
        m_para_a(1.0), m_para_b(1.0), m_para_c(1.0), m_para_d(1.0),
        m_pimage(0), m_plevelset(0)
    {
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Constructor.
    ///
    ///  @param verbose          Verbose level (0..2).
    ///  @param dbg_image_scale  Scale of debug images.
    ///
    ///////////////////////////////////////////////////////////////////////
    metamorphs2(int verbose, const double& dbg_image_scale) :
        m_verbose(verbose), m_debug_image_scale(dbg_image_scale),
        m_para_a(1.0), m_para_b(1.0), m_para_c(1.0), m_para_d(1.0),
        m_pimage(0), m_plevelset(0)
    {
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Creates the model.
    ///
    ///  @param[in]      image          The input image.
    ///  @param[in]      edge_mask      The edge mask.
    ///  @param[in,out]  levelset       The level set function.
    ///  @param[in]      para_a         Weight of interior shape cost term.
    ///  @param[in]      para_b         Weight of boundary shape cost term.
    ///  @param[in]      para_c         Weight of ROI intensity cost term.
    ///  @param[in]      para_d         Weight of maximum-likelihood
    ///                                 intensity cost term.
    ///
    ///////////////////////////////////////////////////////////////////////
    void create(const image_base_type&      image,
                const mask_image_base_type& edge_mask,
                real_image_type&            levelset,
                const double&               para_a,
                const double&               para_b,
                const double&               para_c,
                const double&               para_d
                );

    ///////////////////////////////////////////////////////////////////////
    ///  Creates the model.
    ///
    ///  @param[in]      task_name      Task name, used as the prefix to
    ///                                 the name of output debug images.
    ///  @param[in]      image          The input image.
    ///  @param[in]      edge_mask      The edge mask.
    ///  @param[in,out]  levelset       The level set function.
    ///  @param[in]      para_a         Weight of interior shape cost term.
    ///  @param[in]      para_b         Weight of boundary shape cost term.
    ///  @param[in]      para_c         Weight of ROI intensity cost term.
    ///  @param[in]      para_d         Weight of maximum-likelihood
    ///                                 intensity cost term.
    ///  @param[in]      verbose        Verbose level.
    ///
    ///////////////////////////////////////////////////////////////////////
    void create(const std::string&          task_name,
                const image_base_type&      image,
                const mask_image_base_type& edge_mask,
                real_image_type&            levelset,
                const double&               para_a,
                const double&               para_b,
                const double&               para_c,
                const double&               para_d,
                int                         verbose
                )
    {
        m_verbose = verbose;
        m_task_name = task_name;
        create(image, edge_mask, levelset, para_a, para_b, para_c, para_d);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Model evolution.
    ///
    ///  @param lambda  Step size.
    ///  @param ctrl_0  X-size of the grid of control points.
    ///  @param ctrl_1  Y-size of the grid of control points.
    ///  @param niters  Maximum number of iterations (-1 meaning no limit).
    ///  @param stop_d  Stopping tolerance.
    ///
    ///////////////////////////////////////////////////////////////////////
    void solve(const double& lambda,
               int           ctrl_0,
               int           ctrl_1,
               int           niters = -1,
               const double& stop_d = .1
               );


private:
    
    typedef fast_marching2<real_value_type,
                           real_value_type>    fmm_type;

    void create_shape_image(const mask_image_base_type& edge_mask);
    void create_roi_shape_image(const mask_image_base_type& roi_edge_mask);
    void create_image_gradients(const image_base_type& image);

    /////////////////////////////////////////////////////////////////////// 
    int                     m_verbose;
    double                  m_debug_image_scale;
    double                  m_para_a, m_para_b, m_para_c, m_para_d;
    const image_base_type*  m_pimage;
    real_image_type*        m_plevelset;
    real_image_type         m_shape_image;
    real_image_type         m_dx_image, m_dy_image;
    real_image_type         m_dx_shape_image, m_dy_shape_image;
    real_image_type         m_dx_roi_shape_image, m_dy_roi_shape_image;
    std::string             m_task_name;
    debug                   m_dbg;

};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_alpha/metamorphs2.cxx>
#   endif

#endif // ___METAMORPHS2_HXX___
