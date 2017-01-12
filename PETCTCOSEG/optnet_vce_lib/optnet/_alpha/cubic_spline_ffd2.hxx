/*
 ==========================================================================
 |   
 |   $Id: cubic_spline_ffd2.hxx 186 2005-02-27 02:04:22Z kangli $
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
  Cubic-spline-based 2-D Free-Form Deformation.
 ==========================================================================
 */

#ifndef ___CUBIC_SPLINE_FFD2_HXX___
#   define ___CUBIC_SPLINE_FFD2_HXX___

#   include <optnet/_alpha/cubic_spline.hxx>
#   include <optnet/_base/array2.hxx>
#   include <optnet/_base/array2_ref.hxx>
#   include <optnet/_base/point2.hxx>
#   include <vector>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class cubic_spline_ffd2
///  @brief Cubic-spline-based 2-D Free-Form Deformation.
///////////////////////////////////////////////////////////////////////////
template <typename _Ty>
class cubic_spline_ffd2
{
public:
    
    typedef _Ty                         value_type;
    typedef value_type&                 reference;
    typedef const value_type&           const_reference;

    typedef point2<_Ty>                 point_type;
    typedef point_type&                 point_reference;
    typedef const point_type&           point_const_reference;

    typedef array2_base<value_type>     value_array_base_type;
    typedef array2_ref<value_type>      value_array_ref_type;
    typedef array2<value_type>          value_array_type;

    typedef array2_base<point_type>     array_base_type;
    typedef array2_ref<point_type>      array_ref_type;
    typedef array2<point_type>          array_type;

    typedef size_t                      size_type;

    ///////////////////////////////////////////////////////////////////////
    ///  Constructor.
    ///
    ///  @param control_points  Array of control points.
    ///  @param image_size_x    X-size of the image.
    ///  @param image_size_y    Y-size of the image.
    ///
    ///////////////////////////////////////////////////////////////////////
    cubic_spline_ffd2(const array_base_type& control_points,
                      size_type              image_size_x,
                      size_type              image_size_y
                      )
    {
        m_pq = &control_points;
        m_qx = control_points.size_0() - 1;
        m_qy = control_points.size_1() - 1;
        m_dx = 1.0 / image_size_x;
        m_dy = 1.0 / image_size_y;
        m_sx = image_size_x;
        m_sy = image_size_y;
        m_precomputed = false;
    }
    
    ///////////////////////////////////////////////////////////////////////
    ///  Precompute B-spline coefficients.
    ///////////////////////////////////////////////////////////////////////
    void precompute()
    {
        size_type x, y, i, j;
        double    di, dj;

        m_coe0.create(4,  m_sx);
        m_coe1.create(4,  m_sy);
        
        m_idx0.resize(m_sx    );
        m_id00.resize(m_qx + 1);
        m_id01.resize(m_qx + 1);

        m_idx1.resize(m_sy    );
        m_id10.resize(m_qy + 1);
        m_id11.resize(m_qy + 1);

        for (y = 0; y < m_sy; ++y) {
            dj = (y + 1) * m_dy * m_qy;
            j = (size_type)(dj);
            m_cy.evaluate(dj - j);
            m_idx1[y] = j;
            m_coe1(0, y) = m_cy.b[0];
            m_coe1(1, y) = m_cy.b[1];
            m_coe1(2, y) = m_cy.b[2];
            m_coe1(3, y) = m_cy.b[3];
        } // y

        for (x = 0; x < m_sx; ++x) {
            di = (x + 1) * m_dx * m_qx;
            i = (size_type)(di);
            m_cx.evaluate(di - i);
            m_idx0[x] = i;
            m_coe0(0, x) = m_cx.b[0];
            m_coe0(1, x) = m_cx.b[1];
            m_coe0(2, x) = m_cx.b[2];
            m_coe0(3, x) = m_cx.b[3];
        } // x

        di = m_sx / (double)m_qx;
        for (i = 0; i <= m_qx; ++i) {
            size_type ii;
            ii = i > 3 ? i - 3 : 0;
            m_id00[i] = (size_type)(ii * di);
            ii = i < m_qx ? i + 1 : m_qx;
            m_id01[i] = (size_type)(ii * di);
        } // i

        dj = m_sy / (double)m_qy;
        for (j = 0; j <= m_qy; ++j) {
            size_type jj;
            jj = j > 3 ? j - 3 : 0;
            m_id10[j] = (size_type)(jj * dj);
            jj = j < m_qy ? j + 1 : m_qy;
            m_id11[j] = (size_type)(jj * dj);
        } // j

        m_precomputed = true;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Computes the displacement at image pixels.
    ///
    ///  @param[out] output The output array containing the displacements at
    ///                     each image pixel.
    ///
    ///////////////////////////////////////////////////////////////////////
    void compute(array_base_type& output)
    {
        double    dx, dy;
        size_type i, j, ik, jl;

        if (output.size_0() != m_sx ||
            output.size_1() != m_sy
            )
        {
            throw_exception(std::invalid_argument(
                "cubic_spline_ffd2::compute: "
                "The size of the output array does not match the given value."
                )
            );
        }

        if (!m_precomputed) precompute();

        for (size_type y = 0; y < m_sy; ++y) {
            j = m_idx1[y];

            for (size_type x = 0; x < m_sx; ++x) {
                i = m_idx0[x];

                dx = dy = 0;
                for (size_type l = 0; l < 4; ++l) {
                    jl = j + l;
                    if (jl > m_qy) jl = m_qy;
                    for (size_type k = 0; k < 4; ++k) {
                        ik = i + k;
                        if (ik > m_qx) ik = m_qx;
                        double tmp = m_coe0(k, x) * m_coe1(l, y);
                        dx += (*m_pq)(ik, jl).v[0] * tmp;
                        dy += (*m_pq)(ik, jl).v[1] * tmp;
                    } // k
                } // l

                output(x, y).v[0] = (value_type)(dx);
                output(x, y).v[1] = (value_type)(dy);

            } // x
        } // y
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Compute the displacement at a given point.
    ///
    ///  @param[in]  x             The x-coordinate of the point.
    ///  @param[in]  y             The y-coordinate of the point.
    ///  @param[out] displacement  The output displacement of the point.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void compute_at(double          x,
                           double          y,
                           point_reference displacement
                           )
    {
        double    dx = 0, dy = 0;
        double    di = (x + 1) * m_dx * m_qx;
        double    dj = (y + 1) * m_dy * m_qy;
        size_type ni = static_cast<size_type>(di);
        size_type nj = static_cast<size_type>(dj);
        
        m_cx.evaluate(di - ni);
        m_cy.evaluate(dj - nj);

        for (size_type l = 0; l < 4; ++l) {
            size_type jl = nj + l;
            if (jl > m_qy) jl = m_qy;

            for (size_type k = 0; k < 4; ++k) {
                size_type ik = ni + k;
                if (ik > m_qx) ik = m_qx;
                double tmp = m_cx.b[k] * m_cy.b[l];

                dx += (*m_pq)(ik, jl).v[0] * tmp;
                dy += (*m_pq)(ik, jl).v[1] * tmp;

            } // k
        } // l

        displacement.v[0] = (value_type)(dx);
        displacement.v[1] = (value_type)(dy);
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

    inline double get_derivative(size_type k,
                                 size_type l,
                                 size_type x,
                                 size_type y
                                 ) const
    {
        return m_coe0(k, x) * m_coe1(l, y);
    }


private:
    typedef array2<double>          coeff_array_type;
    typedef std::vector<size_type>  index_array_type;

    coeff_array_type                m_coe0, m_coe1;
    index_array_type                m_idx0, m_id00, m_id01;
    index_array_type                m_idx1, m_id10, m_id11;

    const array_base_type*          m_pq;
    cubic_spline<double>            m_cx, m_cy;
    size_type                       m_qx, m_qy, m_sx, m_sy;
    double                          m_dx, m_dy;
    bool                            m_precomputed;

};

} // namespace

#endif // ___CUBIC_SPLINE_FFD2_HXX___
