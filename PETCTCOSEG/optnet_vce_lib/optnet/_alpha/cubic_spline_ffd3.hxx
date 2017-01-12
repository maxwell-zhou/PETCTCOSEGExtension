/*
 ==========================================================================
 |   
 |   $Id: cubic_spline_ffd3.hxx 186 2005-02-27 02:04:22Z kangli $
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
  Cubic-spline-based 3-D Free-Form Deformation.
 ==========================================================================
 */

#ifndef ___CUBIC_SPLINE_FFD3_HXX___
#   define ___CUBIC_SPLINE_FFD3_HXX___

#   include <optnet/config.h>
#   include <optnet/_alpha/cubic_spline.hxx>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_base/array2.hxx>
#   include <optnet/_base/point3.hxx>
#   include <vector>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class cubic_spline_ffd3
///  @brief Cubic-spline-based 3-D Free-Form Deformation.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Tg = net_f_xy>
class cubic_spline_ffd3
{
public:
    
    typedef _Ty                         value_type;
    typedef value_type&                 reference;
    typedef const value_type&           const_reference;

    typedef point3<_Ty>                 point_type;
    typedef point_type&                 point_reference;
    typedef const point_type&           point_const_reference;

    typedef array_base<value_type>      value_array_base_type;
    typedef array_ref<value_type>       value_array_ref_type;
    typedef array<value_type>           value_array_type;

    typedef array_base<point_type, _Tg> array_base_type;
    typedef array_ref<point_type, _Tg>  array_ref_type;
    typedef array<point_type, _Tg>      array_type;

    typedef size_t                      size_type;

    ///////////////////////////////////////////////////////////////////////
    ///  Constructor.
    ///
    ///  @param control_points  Array of control points.
    ///  @param image_size_x    X-size of the image.
    ///  @param image_size_y    Y-size of the image.
    ///  @param image_size_z    Z-size of the image.
    ///
    ///////////////////////////////////////////////////////////////////////
    cubic_spline_ffd3(const array_base_type& control_points,
                      size_type              image_size_x,
                      size_type              image_size_y,
                      size_type              image_size_z
                      )
    {
        m_pq = &control_points;
        m_qx = control_points.size_0() - 1;
        m_qy = control_points.size_1() - 1;
        m_qz = control_points.size_2() - 1;
        m_dx = 1.0 / image_size_x;
        m_dy = 1.0 / image_size_y;
        m_dz = 1.0 / image_size_z;
        m_sx = image_size_x;
        m_sy = image_size_y;
        m_sz = image_size_z;
        m_precomputed = false;
    }
    
    ///////////////////////////////////////////////////////////////////////
    ///  Precompute B-spline coefficients.
    ///////////////////////////////////////////////////////////////////////
    void precompute()
    {
        size_type x, y, z, i, j, k;
        double    di, dj, dk;

        m_coe0.create(4,  m_sx);
        m_coe1.create(4,  m_sy);
        m_coe2.create(4,  m_sz);
        
        m_idx0.resize(m_sx    );
        m_id00.resize(m_qx + 1);
        m_id01.resize(m_qx + 1);
        
        m_idx1.resize(m_sy    );
        m_id10.resize(m_qy + 1);
        m_id11.resize(m_qy + 1);
        
        m_idx2.resize(m_sz    );
        m_id20.resize(m_qz + 1);
        m_id21.resize(m_qz + 1);

        #ifdef __OPTNET_PRAGMA_OMP__
        #   pragma omp parallel sections \
                num_threads(__OPTNET_OMP_NUM_THREADS__)
            {
        #endif

        #ifdef __OPTNET_PRAGMA_OMP__
        #   pragma omp section
        #endif
        for (z = 0; z < m_sz; ++z) {
            dk = (z + 1) * m_dz * m_qz;
            k  = static_cast<size_type>(dk);
            m_cz.evaluate(dk - k);
            m_idx2[z]    = k;
            m_coe2(0, z) = m_cz.b[0];
            m_coe2(1, z) = m_cz.b[1];
            m_coe2(2, z) = m_cz.b[2];
            m_coe2(3, z) = m_cz.b[3];
        } // z

        #ifdef __OPTNET_PRAGMA_OMP__
        #   pragma omp section
        #endif
        for (y = 0; y < m_sy; ++y) {
            dj = (y + 1) * m_dy * m_qy;
            j  = static_cast<size_type>(dj);
            m_cy.evaluate(dj - j);
            m_idx1[y]    = j;
            m_coe1(0, y) = m_cy.b[0];
            m_coe1(1, y) = m_cy.b[1];
            m_coe1(2, y) = m_cy.b[2];
            m_coe1(3, y) = m_cy.b[3];
        } // y

        #ifdef __OPTNET_PRAGMA_OMP__
        #   pragma omp section
        #endif
        for (x = 0; x < m_sx; ++x) {
            di = (x + 1) * m_dx * m_qx;
            i = static_cast<size_type>(di);
            m_cx.evaluate(di - i);
            m_idx0[x]    = i;
            m_coe0(0, x) = m_cx.b[0];
            m_coe0(1, x) = m_cx.b[1];
            m_coe0(2, x) = m_cx.b[2];
            m_coe0(3, x) = m_cx.b[3];
        } // x

        #ifdef __OPTNET_PRAGMA_OMP__
        #   pragma omp section
        #endif
        {
            di = m_sx / (double)m_qx;
            for (i = 0; i <= m_qx; ++i) {
                size_type ii;
                ii = i > 3 ? i - 3 : 0;
                m_id00[i] = (size_type)(ii * di);
                ii = i < m_qx ? i + 1 : m_qx;
                m_id01[i] = (size_type)(ii * di);
            } // i
        }

        #ifdef __OPTNET_PRAGMA_OMP__
        #   pragma omp section
        #endif
        {
            dj = m_sy / (double)m_qy;
            for (j = 0; j <= m_qy; ++j) {
                size_type jj;
                jj = j > 3 ? j - 3 : 0;
                m_id10[j] = (size_type)(jj * dj);
                jj = j < m_qy ? j + 1 : m_qy;
                m_id11[j] = (size_type)(jj * dj);
            } // j
        }

        #ifdef __OPTNET_PRAGMA_OMP__
        #   pragma omp section
        #endif
        {
            dk = m_sz / (double)m_qz;
            for (k = 0; k <= m_qz; ++k) {
                size_type kk;
                kk = k > 3 ? k - 3 : 0;
                m_id20[k] = (size_type)(kk * dk);
                kk = k < m_qz ? k + 1 : m_qz;
                m_id21[k] = (size_type)(kk * dk);
            } // k
        }

        #ifdef __OPTNET_PRAGMA_OMP__
            }
        #endif
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Computes the displacement at image voxels.
    ///
    ///  @param[out] output The output array containing the displacements at
    ///                     each image voxel.
    ///
    ///////////////////////////////////////////////////////////////////////
    void compute(array_base_type& output)
    {
        int x, y, z;

        if (output.size_0() != m_sx ||
            output.size_1() != m_sy ||
            output.size_2() != m_sz
            )
        {
            throw_exception(std::invalid_argument(
                "cubic_spline_ffd3::compute: "
                "The size of the output array does not match the value that was given."
                )
            );
        }

        if (!m_precomputed) precompute();

        #ifdef __OPTNET_PRAGMA_OMP__
        #   pragma omp parallel \
                num_threads(__OPTNET_OMP_NUM_THREADS__) private(x, y, z)
            {
        #   pragma omp for
        #endif
        for (z = 0; z < (int)m_sz; ++z) {
            size_type k = m_idx2[z];
            for (y = 0; y < (int)m_sy; ++y) {
                size_type j = m_idx1[y];
                for (x = 0; x < (int)m_sx; ++x) {
                    size_type i = m_idx0[x];
                    double dx = 0, dy = 0, dz = 0;
                    for (size_type n = 0; n < 4; ++n) {
                        size_type kn = k + n;
                        if (kn > m_qz) kn = m_qz;
                        for (size_type m = 0; m < 4; ++m) {
                            size_type jm = j + m;
                            if (jm > m_qy) jm = m_qy;
                            for (size_type l = 0; l < 4; ++l) {
                                size_type il = i + l;
                                if (il > m_qx) il = m_qx;
                                double tmp = m_coe0(l, x) * m_coe1(m, y) * m_coe2(n, z);
                                dx += (*m_pq)(il, jm, kn).v[0] * tmp;
                                dy += (*m_pq)(il, jm, kn).v[1] * tmp;
                                dz += (*m_pq)(il, jm, kn).v[2] * tmp;
                            } // k
                        } // l
                    } // m
                    output(x, y, z).v[0] = (value_type)(dx);
                    output(x, y, z).v[1] = (value_type)(dy);
                    output(x, y, z).v[2] = (value_type)(dz);
                } // x
            } // y
        } // z
        #ifdef __OPTNET_PRAGMA_OMP__
            }
        #endif
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Compute the displacement at a given point.
    ///
    ///  @param[in]  x             The x-coordinate of the point.
    ///  @param[in]  y             The y-coordinate of the point.
    ///  @param[in]  z             The z-coordinate of the point.
    ///  @param[out] displacement  The output displacement of the point.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void compute_at(double          x,
                           double          y,
                           double          z,
                           point_reference displacement
                           )
    {
        double    dx = 0, dy = 0, dz = 0;
        double    di = (x + 1) * m_dx * m_qx;
        double    dj = (y + 1) * m_dy * m_qy;
        double    dk = (z + 1) * m_dz * m_qz;
        size_type i = static_cast<size_type>(di);
        size_type j = static_cast<size_type>(dj);
        size_type k = static_cast<size_type>(dk);
        
        m_cx.evaluate(di - i);
        m_cy.evaluate(dj - j);
        m_cz.evaluate(dk - k);

        for (size_type n = 0; n < 4; ++n) {
            size_type kn = k + n;
            if (kn > m_qz) kn = m_qz;
            for (size_type m = 0; m < 4; ++m) {
                size_type jm = j + m;
                if (jm > m_qy) jm = m_qy;
                double tmp1 = m_cy.b[m] * m_cz.b[n];
                for (size_type l = 0; l < 4; ++l) {
                    size_type il = i + l;
                    if (il > m_qx) il = m_qx;
                    double tmp2 = m_cx.b[l] * tmp1;
                    dx += (*m_pq)(il, jm, kn).v[0] * tmp2;
                    dy += (*m_pq)(il, jm, kn).v[1] * tmp2;
                    dz += (*m_pq)(il, jm, kn).v[2] * tmp2;

                } // k
            } // l
        } // m

        displacement.v[0] = (value_type)(dx);
        displacement.v[1] = (value_type)(dy);
        displacement.v[2] = (value_type)(dz);
    }

    ///////////////////////////////////////////////////////////////////////
    inline const size_type& get_index_0(size_type i0) const
                                { return m_idx0[i0]; }
    inline const size_type& get_bound_begin_0(size_type c0) const
                                { return m_id00[c0]; }
    inline const size_type& get_bound_end_0(size_type c0) const
                                { return m_id01[c0]; }

    inline const size_type& get_index_1(size_type i1) const
                                { return m_idx1[i1]; }
    inline const size_type& get_bound_begin_1(size_type c1) const
                                { return m_id10[c1]; }
    inline const size_type& get_bound_end_1(size_type c1) const
                                { return m_id11[c1]; }

    inline const size_type& get_index_2(size_type i2) const
                                { return m_idx2[i2]; }
    inline const size_type& get_bound_begin_2(size_type c2) const
                                { return m_id20[c2]; }
    inline const size_type& get_bound_end_2(size_type c2) const
                                { return m_id21[c2]; }

    inline double get_derivative(size_type l,
                                 size_type m,
                                 size_type n,
                                 size_type x,
                                 size_type y,
                                 size_type z
                                 ) const
    {
        return m_coe0(l, x) * m_coe1(m, y) * m_coe2(n, z);
    }


private:
    typedef array2<double>          coeff_array_type;
    typedef std::vector<size_type>  index_array_type;

    coeff_array_type                m_coe0, m_coe1, m_coe2;
    index_array_type                m_idx0, m_id00, m_id01;
    index_array_type                m_idx1, m_id10, m_id11;
    index_array_type                m_idx2, m_id20, m_id21;

    const array_base_type*          m_pq;
    cubic_spline<double>            m_cx, m_cy, m_cz;
    size_type                       m_qx, m_qy, m_qz, m_sx, m_sy, m_sz;
    double                          m_dx, m_dy, m_dz;
    bool                            m_precomputed;
};

} // namespace

#endif // ___CUBIC_SPLINE_FFD3_HXX___
