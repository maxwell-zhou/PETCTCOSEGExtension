/* 
 ==========================================================================
 |   
 |   $Id: hessian.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___HESSIAN_HXX___
#   define ___HESSIAN_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#      pragma once
#      pragma warning(disable: 4786)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/except.hxx>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_utils/gaussian.hxx>
#   include <optnet/_utils/index.hxx>
#   include <cassert>

/// @namespace optnet
namespace optnet {
    /// @namespace optnet::utils
    namespace utils {

#   define _HESSIAN_SIZE 4

///////////////////////////////////////////////////////////////////////////
/// @class hessian
/// @brief A class that evaluate the hessian matrix and related statistics
///        of an image.
///////////////////////////////////////////////////////////////////////////
template <typename _Tx,
          typename _Ty = float,
          typename _Tg = net_f_xy>
class hessian
{
public:
    
    typedef optnet::array_base<_Tx, _Tg>    array_base_type;
    typedef optnet::array_ref <_Tx, _Tg>    array_ref_type;
    typedef optnet::array     <_Tx, _Tg>    array_type;

    typedef _Tx                             input_value_type;
    typedef _Ty                             output_value_type;

    /// Flags selecting which output values be computed.
    enum   output_flag
    {
        COMPUTE_NONE      = 0,
        COMPUTE_GRADIENTS = 0x01,
        COMPUTE_HESSIAN   = 0x02,
        COMPUTE_EIGEN     = 0x04,
        COMPUTE_ALL       = 0xFF
    };

    /// Structure containing output image statistics.
    struct output_info
    {
        /// Flags that indicate which values should be computed.
        output_flag       flags;
        /// The gradients at the current point.
        output_value_type gradients     [_HESSIAN_SIZE];
        /// The eigen values at the current point.
        output_value_type eigen_values  [_HESSIAN_SIZE];
        /// The eigen vectors at the current point.
        output_value_type eigen_vectors [_HESSIAN_SIZE][_HESSIAN_SIZE];
        /// The hessian matrix at the current point.
        output_value_type hessian_matrix[_HESSIAN_SIZE][_HESSIAN_SIZE];
    };


    typedef size_t                          size_type;


    //////////////////////////////////////////////////////////////////////////
    /// Constructor
    //////////////////////////////////////////////////////////////////////////
    hessian(const output_value_type& sigma,
            const output_value_type& scale
            ) :
        m_input(0), m_alpha(3), m_sigma(sigma), m_scale(scale),
        m_plane(0), m_ghalf(0), m_gsize(0)
    {
        m_range = static_cast<int>(m_sigma * m_scale + 0.5);
        
        for (int i = 0; i < 3; ++i) {
            m_gauss[i] = 0;
        }
        
        init(); // Initialize the gaussian kernels.
    }

    //////////////////////////////////////////////////////////////////////////
    ///
    //////////////////////////////////////////////////////////////////////////
    hessian(const array_base_type&   input, 
            const output_value_type& sigma,
            const output_value_type& scale
            ) :
        m_input(&input), m_alpha(3), m_sigma(sigma), m_scale(scale),
        m_plane(0), m_ghalf(0), m_gsize(0)
    {
        m_range = static_cast<int>(m_sigma * m_scale + 0.5);
        
        for (int i = 0; i < 3; ++i) {
            m_gauss[i] = 0;
        }
        
        init(); // Initialize the gaussian kernels.
    }

    //////////////////////////////////////////////////////////////////////////
    ///
    //////////////////////////////////////////////////////////////////////////
    ~hessian()
    {
        // Anything to clean up here?
    }

    //////////////////////////////////////////////////////////////////////////
    ///
    //////////////////////////////////////////////////////////////////////////
    inline void set_input(const array_base_type& input) {
        m_input = &input;
    }

    //////////////////////////////////////////////////////////////////////////
    ///
    //////////////////////////////////////////////////////////////////////////
    void set_sigma(const output_value_type& sigma) {

        assert(0 != m_input);

        if (sigma < 1.0) {
            throw_exception(
                std::invalid_argument(
                "hessian::set_sigma: The standard deviation cannot be less than 1.0."
                )
            );
        }

        if (m_sigma != sigma) {
            m_range = static_cast<int>(sigma * m_scale + 0.5);
            m_sigma = sigma;
            init();
        }
    }

    //////////////////////////////////////////////////////////////////////////
    ///
    //////////////////////////////////////////////////////////////////////////
    void set_scale(const output_value_type& scale) {

        assert(0 != m_input);

        if (scale < 1.0) {
            throw_exception(
                std::invalid_argument(
                    "hessian::set_scale: The scale cannot be less than 1.0."
                )
            );
        }

        if (m_scale != scale) {
            m_range = static_cast<int>(m_sigma * scale + 0.5);
            m_scale = scale;
            init();
        }
    }

    //////////////////////////////////////////////////////////////////////////
    ///
    //////////////////////////////////////////////////////////////////////////
    void compute_at(const  size_type& i0,
                    const  size_type& i1,
                    const  size_type& i2,
                    output_info& info
                    )
    {
        assert(0 != m_input);
        assert(0 != m_gsize);
        assert(0 != m_ghalf);

        if (static_cast<size_type>(m_gsize) > m_input->size_0() ||
            static_cast<size_type>(m_gsize) > m_input->size_1() ||
            static_cast<size_type>(m_gsize) > m_input->size_2())
        {
            throw_exception(std::range_error(
                "hessian::compute_at: Index out of range."
                ));
        }

        int                x, y, z;
        output_flag        flags;
        output_value_type  temp, temp2;

        flags = info.flags; // Save flags.
        memset(&info, 0, sizeof(output_info));
        info.flags = flags; // Restore flags.

        // Computes nothing.
        if (flags == 0) return;

        // Compute first and second derivatives.
        
        // --- Next compute: Dx, Dy, Dxx, Dxy, Dyx, Dyy
        for (y = 0; y < m_gsize; ++y) {
            for (x = 0; x < m_gsize; ++x) {

                for (temp = (output_value_type)0.0, z = 0; z < m_gsize; ++z) {
                    temp += (*m_input)(
                        index_symmetric((int)i0 + x - m_ghalf, 
                                        (int)m_input->size_0()),
                        index_symmetric((int)i1 + y - m_ghalf,
                                        (int)m_input->size_1()),
                        index_symmetric((int)i2 + z - m_ghalf,
                                        (int)m_input->size_2())
                        ) * m_gauss[0][z];  // Gaussian
                }
                
                m_plane[y * m_gsize + x] = temp;
            } // for x
        } // for y

        if (flags & COMPUTE_GRADIENTS) {
            // --- Dx
            for (temp2 = (output_value_type)0.0, y = 0; y < m_gsize; ++y) {

                for (temp = (output_value_type)0.0, x = 0; x < m_gsize; ++x) {
                    temp += m_plane[y * m_gsize + x]
                              * m_gauss[1][x];  // First derivative of Gaussian
                }
                temp2 += temp * m_gauss[0][y];  // Gaussian
            } // for y

            info.gradients[1] = temp2;

            // --- Dy
            for (temp2 = (output_value_type)0.0, y = 0; y < m_gsize; ++y) {

                for (temp = (output_value_type)0.0, x = 0; x < m_gsize; ++x) {
                    temp += m_plane[y * m_gsize + x]
                              * m_gauss[0][x];  // Gaussian
                }
                temp2 += temp * m_gauss[1][y];  // First derivative of Gaussian
            } // for y

            info.gradients[2] = temp2;
        } // if

        if (flags & COMPUTE_HESSIAN) {
            // --- Dxx
            for (temp2 = (output_value_type)0.0, y = 0; y < m_gsize; ++y) {

                for (temp = (output_value_type)0.0, x = 0; x < m_gsize; ++x) {
                    temp += m_plane[y * m_gsize + x]
                              * m_gauss[2][x];  // Second derivative of Gaussian
                }
                temp2 += temp * m_gauss[0][y];  // Gaussian
            } // for y

            info.hessian_matrix[1][1] = temp2;

            // --- Dyy
            for (temp2 = (output_value_type)0.0, y = 0; y < m_gsize; ++y) {

                for (temp = (output_value_type)0.0, x = 0; x < m_gsize; ++x) {
                    temp += m_plane[y * m_gsize + x]
                              * m_gauss[0][x];  // Gaussian
                }
                temp2 += temp * m_gauss[2][y];  // Second derivative of Gaussian
            } // for y

            info.hessian_matrix[2][2] = temp2;

            // --- Dxy
            for (temp2 = (output_value_type)0.0, y = 0; y < m_gsize; ++y) {

                for (temp = (output_value_type)0.0, x = 0; x < m_gsize; ++x) {
                    temp += m_plane[y * m_gsize + x]
                              * m_gauss[1][x];  // First derivative of Gaussian
                }
                temp2 += temp * m_gauss[1][y];  // First derivative of Gaussian
            } // for y

            info.hessian_matrix[1][2] = 
            info.hessian_matrix[2][1] = temp2;
        } // if
        
        // --- Next compute: Dz, Dxz, Dzx, Dyz, Dzy
        for (y = 0; y < m_gsize; ++y) {
            for (x = 0; x < m_gsize; ++x) {

                for (temp = (output_value_type)0.0, z = 0; z < m_gsize; ++z) {
                    temp += (*m_input)(
                        index_symmetric((int)i0 + x - m_ghalf,
                                        (int)m_input->size_0()),
                        index_symmetric((int)i1 + y - m_ghalf,
                                        (int)m_input->size_1()),
                        index_symmetric((int)i2 + z - m_ghalf,
                                        (int)m_input->size_2())
                        ) * m_gauss[1][z];  // First derivative of Gaussian
                }
                
                m_plane[y * m_gsize + x] = temp;
            } // for x
        } // for y

        if (flags & COMPUTE_GRADIENTS) {
            // --- Dz
            for (temp2 = (output_value_type)0.0, y = 0; y < m_gsize; ++y) {

                for (temp = (output_value_type)0.0, x = 0; x < m_gsize; ++x) {
                    temp += m_plane[y * m_gsize + x]
                              * m_gauss[0][x];  // Gaussian
                }
                temp2 += temp * m_gauss[0][y];  // Gaussian
            } // for y

            info.gradients[3] = temp2;
        } // if

        if (flags & COMPUTE_HESSIAN) {

            // --- Dxz
            for (temp2 = (output_value_type)0.0, y = 0; y < m_gsize; ++y) {

                for (temp = (output_value_type)0.0, x = 0; x < m_gsize; ++x) {
                    temp += m_plane[y * m_gsize + x]
                              * m_gauss[1][x];  // First derivative of Gaussian
                }
                temp2 += temp * m_gauss[0][y];  // Gaussian
            } // for y
            
            info.hessian_matrix[1][3] =
            info.hessian_matrix[3][1] = temp2;

            // -- Dyz
            for (temp2 = (output_value_type)0.0, y = 0; y < m_gsize; ++y) {

                for (temp = (output_value_type)0.0, x = 0; x < m_gsize; ++x) {
                    temp += m_plane[y * m_gsize + x]
                              * m_gauss[0][x];  // Gaussian
                }
                temp2 += temp * m_gauss[1][y];  // First derivative of Gaussian
            } // for y
            
            info.hessian_matrix[2][3] =
            info.hessian_matrix[3][2] = temp2;

            // --- Next compute: Dzz
            for (y = 0; y < m_gsize; ++y) {
                for (x = 0; x < m_gsize; ++x) {

                    for (temp = (output_value_type)0.0, z = 0; z < m_gsize; ++z) {
                        temp += (*m_input)(
                            index_symmetric((int)i0 + x - m_ghalf,
                                            (int)m_input->size_0()),
                            index_symmetric((int)i1 + y - m_ghalf,
                                            (int)m_input->size_1()),
                            index_symmetric((int)i2 + z - m_ghalf,
                                            (int)m_input->size_2())
                            ) * m_gauss[2][z];  // Second derivative of Gaussian
                    }
                    
                    m_plane[y * m_gsize + x] = temp;
                } // for x
            } // for y

            for (temp2 = (output_value_type)0.0, y = 0; y < m_gsize; ++y) {

                for (temp = (output_value_type)0.0, x = 0; x < m_gsize; ++x) {
                    temp += m_plane[y * m_gsize + x]
                              * m_gauss[0][x];  // Gaussian
                }
                temp2 += temp * m_gauss[0][y];  // Gaussian
            } // for y

            info.hessian_matrix[3][3] = temp2;

        } // if

        if (flags & COMPUTE_EIGEN) {

            // Compute eigen values/vectors of hessian matrix.
            jacobi(info.hessian_matrix,
                info.eigen_vectors,
                info.eigen_values
                );

            eigsrt(info.eigen_vectors,
                info.eigen_values
                );

        } // if

    }

    //////////////////////////////////////////////////////////////////////////
    ///
    //////////////////////////////////////////////////////////////////////////
    inline const output_value_type& sigma() const { return m_sigma; }

    //////////////////////////////////////////////////////////////////////////
    ///
    //////////////////////////////////////////////////////////////////////////
    inline const output_value_type& scale() const { return m_scale; }
    

private:
    
    //////////////////////////////////////////////////////////////////////////
    void init()
    {
        int half = m_range + m_alpha;
        int size = 2 * half + 1;
        int i;

        assert(half != 0);      

        for (i = 0; i < 3; ++i) {
            if (0 != m_gauss[i]) delete [] m_gauss[i];
            m_gauss[i] = new output_value_type[size];
        }

        for (i = -half; i <= half; ++i) {
            m_gauss[2][i + half] 
                = d2gaussian(static_cast<output_value_type>(i), m_sigma);
            m_gauss[1][i + half] 
                = dgaussian(static_cast<output_value_type>(i), m_sigma);
            m_gauss[0][i + half]
                = gaussian(static_cast<output_value_type>(i), m_sigma);
        }

        m_ghalf = half; m_gsize = size;
        if (0 != m_plane) delete [] m_plane;
        m_plane = new output_value_type[size * size];
    }
    
    //////////////////////////////////////////////////////////////////////////
#   define _HESSIAN_ROTATE(a, i, j, k, l)               \
        g = a[i][j]; h = a[k][l];                       \
        a[i][j] = (output_value_type)(g - s * (h + g * tau)); \
        a[k][l] = (output_value_type)(h + s * (g - h * tau));

    void jacobi(output_value_type a[_HESSIAN_SIZE][_HESSIAN_SIZE],
                output_value_type v[_HESSIAN_SIZE][_HESSIAN_SIZE],
                output_value_type d[_HESSIAN_SIZE]
                )
    {

        int    j, iq, ip, i, n = _HESSIAN_SIZE - 1;
        double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

        b = new double[n + 1];
        z = new double[n + 1];
        
        for (ip = 1; ip <= n; ip++) {
            for (iq = 1; iq <= n; iq++) {
                v[ip][iq] = 0.0;
            }
            v[ip][ip] = 1.0;
        }

        for (ip = 1; ip <= n; ip++) {
            b[ip] = d[ip] = a[ip][ip];
            z[ip] = 0.0;
        }

        for (i = 1; i <= 50; i++) {
            sm = 0.0;
            for (ip = 1; ip <= n - 1; ip++) {
                for (iq = ip + 1; iq <= n; iq++) {
                    sm += fabs(a[ip][iq]);
                }
            }

            if (sm == 0.0) {
                delete [] z;
                delete [] b;
                return;
            }

            if (i < 4) {
                tresh = 0.2f * sm / (n * n);
            } else {
                tresh = 0.0;
            }

            for (ip = 1; ip <= n - 1; ip++) {
                for (iq = ip + 1; iq <= n; iq++) {
                    g = 100.0 * fabs(a[ip][iq]);
                    if ((4 < i)
                    && ((fabs(d[ip]) + g) == fabs(d[ip]))
                    && ((fabs(d[iq]) + g) == fabs(d[iq]))
                    ) {
                        a[ip][iq] = 0.0;
                    }
                    else if (tresh < fabs(a[ip][iq])) {
                        h = d[iq] - d[ip];
                        if ((output_value_type)(fabs(h) + g) == 
                            (output_value_type)(fabs(h))) {
                            t = (a[ip][iq]) / h;
                        }
                        else {
                            theta = 0.5 * h / (a[ip][iq]);
                            t = 1.0 / (fabs(theta) + 
                                sqrt(1.0 + theta * theta));
                            if (theta < 0.0) {
                                t = -t;
                            }
                        }

                        c = 1.0 / sqrt(1 + t * t);
                        s = t * c;
                        tau = s / (1.0 + c);
                        h = t * a[ip][iq];
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= (output_value_type)h;
                        d[iq] += (output_value_type)h;
                        a[ip][iq] = 0.0;
                        for (j = 1; j <= ip - 1; j++) {
                            _HESSIAN_ROTATE(a, j, ip, j, iq);
                        }

                        for (j = ip + 1; j <= iq - 1; j++) {
                            _HESSIAN_ROTATE(a, ip, j, j, iq);
                        }

                        for (j = iq + 1; j <= n; j++) {
                            _HESSIAN_ROTATE(a, ip, j, iq, j);
                        }

                        for (j = 1; j <= n; j++) {
                            _HESSIAN_ROTATE(v, j, ip, j, iq);
                        }

                    }
                }
            }

            for (ip = 1; ip <= n; ip++) {
                b[ip] += z[ip];
                d[ip] = (output_value_type)b[ip];
                z[ip] = 0.0;
            }
        }
    }

#   undef _HESSIAN_ROTATE

    //////////////////////////////////////////////////////////////////////////
    void eigsrt(output_value_type v[_HESSIAN_SIZE][_HESSIAN_SIZE],
                output_value_type d[_HESSIAN_SIZE]
                )
    {

        int    k, j, i, n = _HESSIAN_SIZE - 1;
        double p;

        for (i = 1; i < n; i++) {

            p = d[k = i];
            for (j = i + 1; j <= n; j++) {
                if (p <= d[j]) {
                    p = d[k = j];
                }
            } // for j

            if (k != i) {
                d[k] = d[i];
                d[i] = (output_value_type)p;
                for (j = 1; j <= n; j++) {
                    p = v[j][i];
                    v[j][i] = v[j][k];
                    v[j][k] = (output_value_type)p;
                } // for j
            } // if k

        } // for i

    }

    //////////////////////////////////////////////////////////////////////////
    const array_base_type* m_input;

    int                 m_alpha;
    int                 m_range;
    output_value_type   m_sigma;      // Standard deviation of the Gaussian kernel
    output_value_type   m_scale;      // 
    output_value_type*  m_gauss[3];   // The Gaussian kernels
    output_value_type*  m_plane;
    int                 m_ghalf;
    int                 m_gsize;
};

#   undef _HESSIAN_SIZE

    } // namespace
} // namespace

#endif // ___HESSIAN_HXX___
