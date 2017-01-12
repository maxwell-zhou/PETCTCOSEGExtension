/*
 ==========================================================================
 |   
 |   $Id: metamorphs2.cxx 167 2005-02-11 20:21:55Z kangli $
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

#ifndef ___METAMORPHS2_CXX___
#   define ___METAMORPHS2_CXX___

#   include <optnet/define.h>
#   include <optnet/_alpha/metamorphs2.hxx>
#   include <limits>

#   ifdef max
#       undef max
#   endif

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Real>
void metamorphs2<_Ty, _Real>::create(const image_base_type&      image,
                                     const mask_image_base_type& edge_mask,
                                     real_image_type&            levelset,
                                     const double&               para_a,
                                     const double&               para_b,
                                     const double&               para_c,
                                     const double&               para_d
                                     )
{
    size_type s0 = image.size_0();
    size_type s1 = image.size_1();

    if (s0 != edge_mask.size_0() || s1 != edge_mask.size_1()) {
        throw_exception(std::invalid_argument(
            "metamorphs2::create: "
            "The size of the image must match the size of the edge mask."
            ));
    }

    if (s0 != levelset.size_0() || s1 != levelset.size_1()) {
        throw_exception(std::invalid_argument(
            "metamorphs2::create: "
            "The size of the image must match the size of the level set."
            ));
    }

    m_pimage = &image;
    m_para_a = para_a;
    m_para_b = para_b;
    m_para_c = para_c;
    m_para_d = para_d;

    // Create "shape image".
    if (para_a != 0 || para_b != 0)
        create_shape_image(edge_mask);
    // TODO: Add roi edge shape image term.
    if (para_d != 0) create_image_gradients(image);

    // Save pointer to the level set.
    m_plevelset = &levelset;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Real>
void metamorphs2<_Ty, _Real>::solve(const double& lambda,
                                    int           ctrl_0,
                                    int           ctrl_1,
                                    int           niters,
                                    const double& stop_d
                                    )
{
    using namespace optnet::utils;

    typedef cubic_spline_ffd2<real_value_type> ffd_type;

    int             iter, s0, s1, i0, i1, c0, c1;
    real_value_type iter_d, rayl_b, rayl_b2, inv_rayl_b2 = 0.0;
    
    real_image_type&        levelset = *m_plevelset;
    const image_base_type&  image    = *m_pimage;

    s0 = (int)levelset.size_0();
    s1 = (int)levelset.size_1();

    if (ctrl_0 < 0) ctrl_0 = s0 / 4;
    if (ctrl_1 < 0) ctrl_1 = s1 / 4;
    if (ctrl_0 < 2) ctrl_0 = 2;
    if (ctrl_1 < 2) ctrl_1 = 2;

    if (niters < 0) {
        niters = std::numeric_limits<int>::max();
    }

    m_dbg.log_begin(NULL);

    // Compute initial region statistics.
    if (m_para_d != 0) {
        if (m_verbose > 0)
            m_dbg.log_printf("Calculating region statistics...\n");

        real_value_type         mean;
        std::vector<pixel_type> init_pixels;
        int                     num_pixels;

        mean = 0;
        num_pixels = 0;
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
                
                if (levelset(i0, i1) > 0.0) {
                    const pixel_type& pixel = image(i0, i1);
                    init_pixels.push_back(pixel);
                    mean = mean + pixel;
                    ++num_pixels;
                }

            } // i0
        } // i1

        mean    = mean / num_pixels;
        rayl_b  = mean / (M_SQRTPI_2);
        rayl_b2 = rayl_b * rayl_b;

        if (rayl_b2 != 0) {
            inv_rayl_b2 = (real_value_type)(1.0 / rayl_b2);
        }
        else {
            inv_rayl_b2 = 0.0;
        }

        init_pixels.clear();
    }

    // Initialize Free-Form Deformation (FFD) engine.
    typename ffd_type::array_type q  (ctrl_0, ctrl_1);
    typename ffd_type::array_type dis(s0, s1);

    ffd_type ffd(q, s0, s1);

    if (m_verbose > 0)
        m_dbg.log_printf("Precomputing...\n");
    ffd.precompute();

    for (iter = 0; iter < niters; ++iter) {
    
        // STEP 1: Initialize q to q0.
        q.clear();
        
        // STEP 2/3:
        #ifdef __OPTNET_PRAGMA_OMP__
        #   pragma omp parallel num_threads(__OPTNET_OMP_NUM_THREADS__) \
                private(c0, c1, i0, i1) \
                shared(ctrl_0, ctrl_1, lambda, q)
            {
        #   pragma omp for
        #endif
        for (c1 = 0; c1 < ctrl_1; ++c1) {

            int begin1 = (int)ffd.get_bound_begin_1(c1);
            int end1 = (int)ffd.get_bound_end_1(c1);

            for (c0 = 0; c0 < ctrl_0; ++c0) {

                int begin0 = (int)ffd.get_bound_begin_0(c0);
                int end0 = (int)ffd.get_bound_end_0(c0);

                real_value_type isdt_x = 0, isdt_y = 0;
                real_value_type bsdt_x = 0, bsdt_y = 0;
                real_value_type ridt_x = 0, ridt_y = 0;
                real_value_type midt_x = 0, midt_y = 0;
                int             vrm = 0, vbm = 0;

                for (i1 = begin1; i1 < end1; ++i1) {
                    
                    size_type j = ffd.get_index_1(i1);
                    size_type m = c1 - j;

                    for (i0 = begin0; i0 < end0; ++i0) {
                        
                        real_value_type tmp1, tmp2, dev, dev2;

                        size_type i = ffd.get_index_0(i0);
                        size_type l = c0 - i;

                        if (levelset(i0, i1) > 0.0) {
                            dev  = ffd.get_derivative(l, m, i0, i1);
                            dev2 = dev * 2;

                            // -- interior shape data term:
                            if (m_para_a != 0) {
                                tmp1 = (m_shape_image(i0, i1) - 
                                        levelset(i0, i1)) * dev2;
                                isdt_x += tmp1 * m_dx_shape_image(i0, i1);
                                isdt_y += tmp1 * m_dy_shape_image(i0, i1);
                            }

                            // -- region-of-interest intensity data term
                            if (m_para_c != 0) {
                                //
                                //TODO: Implement this.
                                //
                            }

                            ++vrm;
                        }
                        else if (levelset(i0, i1) > -1.0) {
                            dev  = ffd.get_derivative(l, m, i0, i1);
                            dev2 = dev * 2;

                            // -- boundary shape data term:
                            if (m_para_b != 0) {
                                tmp1 = m_shape_image(i0, i1) * dev2;
                                bsdt_x += tmp1 * m_dx_shape_image(i0, i1);
                                bsdt_y += tmp1 * m_dy_shape_image(i0, i1);
                            }

                            // -- maximum likelihood intensity data term
                            if (m_para_d != 0) {
                                real_value_type inv_voxel, 
                                    voxel = (real_value_type)image(i0, i1);
                                inv_voxel = (voxel == 0) ? 0.0 : 1.0 / voxel;
                                tmp2 = (inv_voxel - inv_rayl_b2) * dev;
                                midt_x += tmp2 * m_dx_image(i0, i1);
                                midt_y += tmp2 * m_dy_image(i0, i1);
                            }

                            ++vbm;
                        }

                    } // for i0
                } // for i1

                if (vrm > 0) {
                    real_value_type inv_vrm = 1.0 / vrm;
                    if (m_para_a != 0) {
                        isdt_x *= inv_vrm;
                        isdt_y *= inv_vrm;
                    }
                    if (m_para_c != 0) {
                        ridt_x *= inv_vrm;
                        ridt_y *= inv_vrm;
                    }
                }

                if (vbm > 0) {
                    real_value_type inv_vbm = 1.0 / vbm;
                    if (m_para_b != 0) {
                        bsdt_x *= inv_vbm;
                        bsdt_y *= inv_vbm;
                    }
                    if (m_para_d != 0) {
                        midt_x *= inv_vbm;
                        midt_y *= inv_vbm;
                    }
                }

                q(c0, c1)[0] = (
                    m_para_a * isdt_x + 
                    m_para_b * bsdt_x +
                    m_para_c * ridt_x +
                    m_para_d * midt_x) * lambda;
                q(c0, c1)[1] = (
                    m_para_a * isdt_y + 
                    m_para_b * bsdt_y +
                    m_para_c * ridt_y +
                    m_para_d * midt_y) * lambda;

            } // c0
        } // c1
        #ifdef __OPTNET_PRAGMA_OMP__
            }
        #endif

        ffd.compute(dis);
        interp_bspline2<real_value_type> interp;
        interp.create(levelset, 1);

        #ifdef __OPTNET_PRAGMA_OMP__
        #   pragma omp parallel num_threads(__OPTNET_OMP_NUM_THREADS__) \
                private(i0, i1) \
                shared(s0, s1, interp, dis)
            {
        #   pragma omp for
        #endif
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
                levelset(i0, i1) = interp.interp(i0 + dis(i0, i1)[0],
                                                 i1 + dis(i0, i1)[1]
                                                 );
            } // i0
        } // i1
        #ifdef __OPTNET_PRAGMA_OMP__
            }
        #endif

        // Check termination condition.
        iter_d = 0.0;
        for (i0 = 0; i0 < (int)dis.size(); ++i0) {
            real_value_type diff = dis[i0][0] * dis[i0][0] +
                                   dis[i0][1] * dis[i0][1];
            if (sqrt(diff) > iter_d)
                iter_d = diff;
        }

        if (m_verbose > 0) {
            m_dbg.log_printf("%05d - %e\n", iter + 1, iter_d);
        }

        if (m_verbose > 1) {
            char extbuf[64];
            secure_sprintf(extbuf, sizeof(extbuf), ".%05d.png", iter + 1);
            contour_png<real_value_type, pixel_type> contour;
            contour.find(std::string(m_task_name + extbuf).c_str(),
                    levelset, *m_pimage, false, 0,
                    m_debug_image_scale
                );
        }

        if (iter_d < stop_d)
            break; // terminate

    } // main loop

    m_dbg.log_end();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Real>
void metamorphs2<_Ty, _Real>::create_shape_image(
    const mask_image_base_type& edge_mask
    )
{
    fmm_type::node_type node;
    fmm_type            fmm;

    size_type s0 = edge_mask.size_0();
    size_type s1 = edge_mask.size_1();
    size_type i0, i1;

    // ================================================================
    // Create shape image. =======================================
    m_shape_image.create(s0, s1);

    fmm.create(m_shape_image, NULL);
    
    for (i1 = 0; i1 < s1; ++i1) {
        for (i0 = 0; i0 < s0; ++i0) {
            if (edge_mask(i0, i1) != 0) {
                node.index[0] = i0;
                node.index[1] = i1;
                node.value    = 0;
                fmm.trial_nodes().push_back(node);
            }
        } // i0
    } // i1

    // Compute Euclidean distance transformation
    // by fast marching.
    fmm.solve();

    if (m_verbose > 1) {
        // Save the edge shape image to file (nasty!).
        using namespace optnet::io;
        array<real_value_type> shape_copy(m_shape_image);
        array<unsigned char> png_copy;
        shape_copy.scale_to_range(0, 255);
        shape_copy.copy_to(png_copy);
        std::string png_name = m_task_name + ".shape.png";
        png_save(png_copy, png_name.c_str());
    }

    // Compute derivatives of the shape image.
    m_dx_shape_image.create(s0, s1);
    m_dy_shape_image.create(s0, s1);

    #ifdef __OPTNET_PRAGMA_OMP__
    #   pragma omp parallel num_threads(__OPTNET_OMP_NUM_THREADS__) \
            private(i0, i1) shared(s0, s1)
        {
    #   pragma omp for
    #endif
    for (i1 = 1; i1 + 1 < s1; ++i1) {
        for (i0 = 1; i0 + 1 < s0; ++i0) {
            
            // Compute derivatives by central difference.
            m_dx_shape_image(i0, i1) = (m_shape_image(i0 + 1, i1) -
                                        m_shape_image(i0 - 1, i1)) * 0.5;
            m_dy_shape_image(i0, i1) = (m_shape_image(i0, i1 + 1) -
                                        m_shape_image(i0, i1 - 1)) * 0.5;

        } // i0
    } // i1
    #ifdef __OPTNET_PRAGMA_OMP__
        }
    #endif
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Real>
void metamorphs2<_Ty, _Real>::create_roi_shape_image(
    const mask_image_base_type& roi_edge_mask
    )
{
    // Not implemented.
    OPTNET_UNUSED(roi_edge_mask);
}

///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Real>
void metamorphs2<_Ty, _Real>::create_image_gradients(
    const image_base_type& image
    )
{
    int i0, i1;
    int s0 = (int)image.size_0();
    int s1 = (int)image.size_1();

    // ================================================================
    // Compute image gradients. =======================================

    m_dx_image.create(s0, s1);
    m_dy_image.create(s0, s1);

    #ifdef __OPTNET_PRAGMA_OMP__
    #   pragma omp parallel num_threads(__OPTNET_OMP_NUM_THREADS__) \
            private(i0, i1) shared(s0, s1, image)
        {
    #   pragma omp for
    #endif
    for (i1 = 1; i1 < s1 - 1; ++i1) {
        for (i0 = 1; i0 < s0 - 1; ++i0) {
            
            // Compute derivatives by central difference.
            m_dx_image(i0, i1)
                = ((real_value_type)image(i0 + 1, i1    ) -
                   (real_value_type)image(i0 - 1, i1    )) * 0.5;
            m_dy_image(i0, i1)
                = ((real_value_type)image(i0,     i1 + 1) -
                   (real_value_type)image(i0,     i1 - 1)) * 0.5;

        } // i0
    } // i1
    #ifdef __OPTNET_PRAGMA_OMP__
        }
    #endif
}

} // namespace

#endif
