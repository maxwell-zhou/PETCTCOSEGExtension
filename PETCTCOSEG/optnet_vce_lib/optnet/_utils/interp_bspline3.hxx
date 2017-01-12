/* 
 ==========================================================================
 |   
 |   $Id: interp_bspline3.hxx 140 2005-02-07 04:23:55Z kangli $
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

#ifndef ___INTERP_BSPLINE3_HXX___
#   define ___INTERP_BSPLINE3_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/except.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_base/array.hxx>
#   include <limits>
#   include <cmath>

/// @namespace optnet
namespace optnet {
    /// @namespace optnet::utils
    namespace utils {

template <typename _Tv,
          typename _Real = double,
          typename _Tg = net_f_xy>
class interp_bspline3
{
public:
    
    typedef _Tv                             value_type;
    typedef _Real                           real_value_type;

    // array types
    typedef optnet::array_base<_Tv, _Tg>    value_array_base_type;
    typedef optnet::array_ref<_Tv, _Tg>     value_array_ref_type;
    typedef optnet::array<_Tv, _Tg>         value_array_type;

    typedef optnet::array_base<_Real>       real_value_array_base_type; 
    typedef optnet::array_ref<_Real>        real_value_array_ref_type;  
    typedef optnet::array<_Real>            real_value_array_type;  

    typedef size_t                          size_type;


    ///////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    ///////////////////////////////////////////////////////////////////////
    interp_bspline3() :
        m_order(0)
    {
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Constructor.
    ///
    ///  @param  a      Input array.
    ///  @param  order  B-spline order.
    ///
    ///////////////////////////////////////////////////////////////////////
    interp_bspline3(const value_array_base_type& a, int order)
    {
        create(a, order);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Creates the coefficient lookup table, initializes the interpolator.
    ///
    ///  @param  a      Input array.
    ///  @param  order  B-spline order.
    ///
    ///////////////////////////////////////////////////////////////////////
    void create(const value_array_base_type& a, int order)
    {
        int                 npole = 0;
        double              poles[2];
        size_type           i0, i1, i2;
        real_value_type*    pline;
        
        // recover the poles from a lookup table
        switch (order) {
        case 0:
        case 1:
            npole = 0;
            break;
        case 2:
            npole = 1;
            poles[0] = sqrt(8.0) - 3.0;
            break;
        case 3:
            npole = 1;
            poles[0] = sqrt(3.0) - 2.0;
            break;
        case 4:
            npole = 2;
            poles[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
            poles[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
            break;
        case 5:
            npole = 2;
            poles[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0))
                     + sqrt(105.0 / 4.0) - 13.0 / 2.0;
            poles[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0))
                     - sqrt(105.0 / 4.0) - 13.0 / 2.0;
            break;
        default:
            throw_exception(
                std::invalid_argument(
                    "interp_bspline3::create: Invalid order.")
                );
        }
        
        // allocate coefficient array
        m_coeff.create(a.size_0(), a.size_1(), a.size_2());
        
        // convert the image samples into interpolation coefficients
        // -- i2
        pline = new real_value_type[a.size_2()];
        for (i1 = 0; i1 < a.size_1(); ++i1) {
            for (i0 = 0; i0 < a.size_0(); ++i0) {
                for (i2 = 0; i2 < a.size_2(); ++i2)
                    pline[i2] = a(i0, i1, i2);
                convert_to_coeff(pline, (int)a.size_2(), poles, npole);
                for (i2 = 0; i2 < a.size_2(); ++i2)
                    m_coeff(i0, i1, i2) = pline[i2];
            }
        }
        delete [] pline;

       // -- i1
        pline = new real_value_type[a.size_1()];
        for (i2 = 0; i2 < a.size_2(); ++i2) {
            for (i0 = 0; i0 < a.size_0(); ++i0) {
                for (i1 = 0; i1 < a.size_1(); ++i1)
                    pline[i1] = m_coeff(i0, i1, i2);
                convert_to_coeff(pline, (int)a.size_1(), poles, npole);
                for (i1 = 0; i1 < a.size_1(); ++i1)
                    m_coeff(i0, i1, i2) = pline[i1];
            }
        }
        delete [] pline;

        // --i0
        for (i2 = 0; i2 < a.size_2(); ++i2) {
            for (i1 = 0; i1 < a.size_1(); ++i1) {
                convert_to_coeff(
                    &(m_coeff(0, i1, i2)), 
                    (int)a.size_0(), 
                    poles, npole
                );
            }
        }

        m_order = order;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Performs the interpolation.
    ///
    ///  @param  i0  The first  coordinate of the point to be interpolated.
    ///  @param  i1  The second coordinate of the point to be interpolated.
    ///  @param  i2  The third  coordinate of the point to be interpolated.
    ///
    ///  @return The interpolated value.
    ///////////////////////////////////////////////////////////////////////
    real_value_type interp(const real_value_type& i0,
                           const real_value_type& i1,
                           const real_value_type& i2
                           )
    {
        double  w, w2, w4, t, t0, t1;
        double  weights[3][6];
        int     indices[3][6];
        int     ss0, ss1, ss2, s0, s1, s2;
        int     n, i, j, k;
        double  ans;

        s0  = (int)m_coeff.size_0();
        s1  = (int)m_coeff.size_1();
        s2  = (int)m_coeff.size_2();
        ss1 = s1 * 2 - 2;
        ss0 = s0 * 2 - 2;
        ss2 = s2 * 2 - 2;

        // compute the interpolation indexes
        if (m_order & 1)  { // odd
            i = (int)floor(i0) - m_order / 2;
            j = (int)floor(i1) - m_order / 2;
            k = (int)floor(i2) - m_order / 2;
            for (n = 0; n <= m_order; ++n) {
                indices[0][n] = i++;
                indices[1][n] = j++;
                indices[2][n] = k++;
            }
        }
        else {              // even
            i = (int)floor(i0 + 0.5) - m_order / 2;
            j = (int)floor(i1 + 0.5) - m_order / 2;
            k = (int)floor(i2 + 0.5) - m_order / 2;
            for (n = 0; n <= m_order; ++n) {
                indices[0][n] = i++;
                indices[1][n] = j++;
                indices[2][n] = k++;
            }
        }
        switch (m_order) {
        case 3:
            w = i0 - (double) indices[0][1];
            weights[0][3] = (1.0 / 6.0) * w * w * w;
            weights[0][0] = (1.0 / 6.0) + 0.5 * w * (w - 1.0) - weights[0][3];
            weights[0][2] = w + weights[0][0] - 2.0 * weights[0][3];
            weights[0][1] = 1.0 - weights[0][0] - weights[0][2] - weights[0][3];

            w = i1 - (double) indices[1][1];
            weights[1][3] = (1.0 / 6.0) * w * w * w;
            weights[1][0] = (1.0 / 6.0) + 0.5 * w * (w - 1.0) - weights[1][3];
            weights[1][2] = w + weights[1][0] - 2.0 * weights[1][3];
            weights[1][1] = 1.0 - weights[1][0] - weights[1][2] - weights[1][3];

            w = i2 - (double) indices[2][1];
            weights[2][3] = (1.0 / 6.0) * w * w * w;
            weights[2][0] = (1.0 / 6.0) + 0.5 * w * (w - 1.0) - weights[2][3];
            weights[2][2] = w + weights[2][0] - 2.0 * weights[2][3];
            weights[2][1] = 1.0 - weights[2][0] - weights[2][2] - weights[2][3];
            break;

        case 0:
            // implements nearest neighbor
            weights[0][0] = 1;
            weights[1][0] = 1;
            weights[2][0] = 1;
            break;

        case 1:
            w = i0 - (double) indices[0][0];
            weights[0][1] = w;
            weights[0][0] = 1.0 - w;

            w = i1 - (double) indices[1][0];
            weights[1][1] = w;
            weights[1][0] = 1.0 - w;

            w = i2 - (double) indices[2][0];
            weights[2][1] = w;
            weights[2][0] = 1.0 - w;
            break;

        case 2:
            w = i0 - (double) indices[0][1];
            weights[0][1] = 0.75 - w * w;
            weights[0][2] = 0.5 * (w - weights[0][1] + 1.0);
            weights[0][0] = 1.0 - weights[0][1] - weights[0][2];

            w = i1 - (double) indices[1][1];
            weights[1][1] = 0.75 - w * w;
            weights[1][2] = 0.5 * (w - weights[1][1] + 1.0);
            weights[1][0] = 1.0 - weights[1][1] - weights[1][2];

            w = i2 - (double) indices[2][1];
            weights[2][1] = 0.75 - w * w;
            weights[2][2] = 0.5 * (w - weights[2][1] + 1.0);
            weights[2][0] = 1.0 - weights[2][1] - weights[2][2];
            break;

        case 4:
            w = i0 - (double) indices[0][2];
            w2 = w * w;
            t = (1.0 / 6.0) * w2;
            weights[0][0] = 0.5 - w;
            weights[0][0] *= weights[0][0];
            weights[0][0] *= (1.0 / 24.0) * weights[0][0];
            t0 = w * (t - 11.0 / 24.0);
            t1 = 19.0 / 96.0 + w2 * (0.25 - t);
            weights[0][1] = t1 + t0;
            weights[0][3] = t1 - t0;
            weights[0][4] = weights[0][0] + t0 + 0.5 * w;
            weights[0][2] = 1.0 -
                weights[0][0] -
                weights[0][1] -
                weights[0][3] -
                weights[0][4];

            w = i1 - (double) indices[1][2];
            w2 = w * w;
            t = (1.0 / 6.0) * w2;
            weights[1][0] = 0.5 - w;
            weights[1][0] *= weights[1][0];
            weights[1][0] *= (1.0 / 24.0) * weights[1][0];
            t0 = w * (t - 11.0 / 24.0);
            t1 = 19.0 / 96.0 + w2 * (0.25 - t);
            weights[1][1] = t1 + t0;
            weights[1][3] = t1 - t0;
            weights[1][4] = weights[1][0] + t0 + 0.5 * w;
            weights[1][2] = 1.0 -
                weights[1][0] -
                weights[1][1] -
                weights[1][3] -
                weights[1][4];

            w = i2 - (double) indices[2][2];
            w2 = w * w;
            t = (1.0 / 6.0) * w2;
            weights[2][0] = 0.5 - w;
            weights[2][0] *= weights[2][0];
            weights[2][0] *= (1.0 / 24.0) * weights[2][0];
            t0 = w * (t - 11.0 / 24.0);
            t1 = 19.0 / 96.0 + w2 * (0.25 - t);
            weights[2][1] = t1 + t0;
            weights[2][3] = t1 - t0;
            weights[2][4] = weights[2][0] + t0 + 0.5 * w;
            weights[2][2] = 1.0 -
                weights[2][0] -
                weights[2][1] -
                weights[2][3] -
                weights[2][4];
            break;

        case 5:
            w = i0 - (double) indices[0][2];
            w2 = w * w;
            weights[0][5] = (1.0 / 120.0) * w * w2 * w2;
            w2 -= w;
            w4 = w2 * w2;
            w -= 0.5;
            t = w2 * (w2 - 3.0);
            weights[0][0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - weights[0][5];
            t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
            t1 = (-1.0 / 12.0) * w * (t + 4.0);
            weights[0][2] = t0 + t1;
            weights[0][3] = t0 - t1;
            t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
            t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
            weights[0][1] = t0 + t1;
            weights[0][4] = t0 - t1;

            w = i1 - (double) indices[1][2];
            w2 = w * w;
            weights[1][5] = (1.0 / 120.0) * w * w2 * w2;
            w2 -= w;
            w4 = w2 * w2;
            w -= 0.5;
            t = w2 * (w2 - 3.0);
            weights[1][0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - weights[1][5];
            t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
            t1 = (-1.0 / 12.0) * w * (t + 4.0);
            weights[1][2] = t0 + t1;
            weights[1][3] = t0 - t1;
            t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
            t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
            weights[1][1] = t0 + t1;
            weights[1][4] = t0 - t1;

            w = i2 - (double) indices[2][2];
            w2 = w * w;
            weights[2][5] = (1.0 / 120.0) * w * w2 * w2;
            w2 -= w;
            w4 = w2 * w2;
            w -= 0.5;
            t = w2 * (w2 - 3.0);
            weights[2][0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - weights[2][5];
            t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
            t1 = (-1.0 / 12.0) * w * (t + 4.0);
            weights[2][2] = t0 + t1;
            weights[2][3] = t0 - t1;
            t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
            t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
            weights[2][1] = t0 + t1;
            weights[2][4] = t0 - t1;
            break;

        default:
            break;
        }

        // apply the mirror boundary conditions
        for (k = 0; k <= m_order; ++k) {
            indices[0][k] = (s0 == 1) ? (0) :
                (
                    (indices[0][k] < 0) ?
                        (-indices[0][k] - ss0 * ((-indices[0][k]) / ss0)) :
                            (indices[0][k] - ss0 * (indices[0][k] / ss0))
                );
            if (s0 <= indices[0][k]) indices[0][k] = ss0 - indices[0][k];

            indices[1][k] = (s1 == 1) ? (0) :
                (
                    (indices[1][k] < 0) ?
                        (-indices[1][k] - ss1 * ((-indices[1][k]) / ss1)) :
                            (indices[1][k] - ss1 * (indices[1][k] / ss1))
                );
            if (s1 <= indices[1][k]) indices[1][k] = ss1 - indices[1][k];

            indices[2][k] = (s2 == 1) ? (0) :
                (
                    (indices[2][k] < 0) ?
                        (-indices[2][k] - ss2 * ((-indices[2][k]) / ss2)) :
                            (indices[2][k] - ss2 * (indices[2][k] / ss2))
                );
            if (s2 <= indices[2][k]) indices[2][k] = ss2 - indices[2][k];
        }

        // perform interpolation
        ans = 0.0;

        for (k = 0; k <= m_order; ++k) {
            w2 = 0.0;
            for (j = 0; j <= m_order; ++j) {
                w = 0.0;
                for (i = 0; i <= m_order; ++i)
                    w += weights[0][i] * m_coeff(indices[0][i],
                                                 indices[1][j],
                                                 indices[2][k]
                                                 );
                w2 += weights[1][j] * w;
            }
            ans += weights[2][k] * w2;
        }
        return (real_value_type)ans;
    }

    //////////////////////////////////////////////////////////////////////////
    ///  Free the lookup table and set the order to zero.
    //////////////////////////////////////////////////////////////////////////
    void release()
    {
        m_coeff.release();
        m_order = 0;
    }

protected:

    //////////////////////////////////////////////////////////////////////////
    void convert_to_coeff(real_value_type* pline,
                          int              count,
                          double*          poles,
                          int              npole) const
    {
        double  lambda = 1.0;
        int     n, k;

        // special case required by mirror boundaries
        if (count == 1) return;

        if (npole > 0) {
            for (k = 0; k < npole; ++k) {
                // compute the overall gain
                lambda = lambda * (1.0 - poles[k]) * (1.0 - 1.0 / poles[k]);
            }

            // apply the gain
            for (n = 0; n < count; ++n) pline[n] *= (real_value_type)lambda;
        }

        // loop over all poles
        for (k = 0; k < npole; ++k) {
            // causal initialization & recursion
            pline[0] = (real_value_type)
                get_initial_causal_coeff(pline, count, poles[k]);
            for (n = 1; n < count; ++n)
                pline[n] += (real_value_type)poles[k] * pline[n - 1];

            // anticausal initialization & recursion
            pline[count - 1] = (real_value_type)
                get_initial_anticausal_coeff(pline, count, poles[k]);
            for (n = count - 2; 0 <= n; --n)
                pline[n] = (real_value_type)poles[k] * (pline[n + 1] - pline[n]);
        }
    }

    //////////////////////////////////////////////////////////////////////////
    double get_initial_causal_coeff(real_value_type* pline,
                                    int              count,
                                    double           dpole) const
    {
        double const  TOL = std::numeric_limits<double>::epsilon();
        double        sum, zn, z2n, iz;
        int           horizon;
        int           n;

        // this initialization corresponds to mirror boundaries
        horizon = count;

        if (TOL > 0.0)
            horizon = (int)ceil(log(TOL) / log(fabs(dpole)));

        if (horizon < count) {
            // accelerated loop
            zn = dpole;
            sum = pline[0];
            for (n = 1; n < horizon; ++n) {
                sum += zn * pline[n];
                zn *= dpole;
            }
            return(sum);
        }
        else {
            // full loop
            zn = dpole;
            iz = 1.0 / dpole;
            z2n = pow(dpole, (double)(count - 1));
            sum = pline[0] + z2n * pline[count - 1];
            z2n *= z2n * iz;
            for (n = 1; n <= count - 2; ++n) {
                sum += (zn + z2n) * pline[n];
                zn *= dpole;
                z2n *= iz;
            }
            return (sum / (1.0 - zn * zn));
        }
    }

    //////////////////////////////////////////////////////////////////////////
    double get_initial_anticausal_coeff(real_value_type* pline,
                                        int              count,
                                        double           dpole) const
    {
        return ((dpole / (dpole * dpole - 1.0)) *
                (dpole * pline[count - 2] + pline[count - 1])
                );
    }

    //////////////////////////////////////////////////////////////////////////
    int                     m_order;
    real_value_array_type   m_coeff;


private:

    // not implemented
    interp_bspline3(const interp_bspline3&);
    interp_bspline3& operator=(const interp_bspline3&);

};

    } // namespace
} // namespace

#endif
