/*
 ==========================================================================
 |   
 |   $Id: contour.hxx 186 2005-02-27 02:04:22Z kangli $
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

#ifndef ___CONTOUR_HXX___
#   define ___CONTOUR_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array2.hxx>
#   include <optnet/_base/array2_ref.hxx>
#   include <optnet/_base/io/png.hxx>
#   include <optnet/_utils/xstring.hxx>
#   include <algorithm>
#   include <string>

namespace optnet {
    namespace detail {

#   ifdef min
#       undef min
#   endif
#   ifdef max
#       undef max
#   endif

#   define XSECT(p1, p2) (h[p2] * xh[p1] - h[p1] * xh[p2]) / (h[p2] - h[p1])
#   define YSECT(p1, p2) (h[p2] * yh[p1] - h[p1] * yh[p2]) / (h[p2] - h[p1])

    ///////////////////////////////////////////////////////////////////////////
    ///  @class contour_base
    ///  @brief Provides 2-D contouring routines.
    ///////////////////////////////////////////////////////////////////////////
    template <typename _Ty>
    class contour_base
    {
    public:

        typedef _Ty                         value_type;
        typedef value_type&                 reference;
        typedef const value_type&           const_reference;

        typedef array2_base<value_type>     array_base_type;
        typedef array2_ref<value_type>      array_ref_type;
        typedef array2<value_type>          array_type;

        typedef size_t                      size_type;


        ///////////////////////////////////////////////////////////////////////
        ///  Default constructor.
        ///////////////////////////////////////////////////////////////////////
        contour_base()
        {
        }

        ///////////////////////////////////////////////////////////////////////
        ///  Default destructor.
        ///////////////////////////////////////////////////////////////////////
        virtual ~contour_base()
        {
        }

    protected:

        ///////////////////////////////////////////////////////////////////////
        ///  Find the contours.
        ///
        ///  @param data        The 2-D array in which to find contours.
        ///  @param levels      The contour levels.
        ///  @param num_levels  Number of levels to find.
        ///////////////////////////////////////////////////////////////////////
        void do_find(const array_base_type& data,
                     double*                levels,
                     int                    num_levels
                     )
        {
            static const int OFS0_TBL[4] = {0, 1, 1, 0};
            static const int OFS1_TBL[4] = {0, 0, 1, 1};
            static const int CASE_TBL[3][3][3] =
            {
                {{0, 0, 8}, {0, 2, 5}, {7, 6, 9}},
                {{0, 3, 4}, {1, 3, 1}, {4, 3, 0}},
                {{9, 6, 7}, {5, 2, 0}, {8, 0, 0}}
            };

            int    i0, i1, k, m, m1, m2, m3, case_value, sh[5];
            double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0;
            double temp1, temp2, dmin, dmax;
            double h[5], xh[5], yh[5];
            
            int size0 = (int)data.size_0();
            int size1 = (int)data.size_1();

            for (i1 = (size1 - 2); i1 >= 0; --i1) {
                for (i0 = 0; i0 + 1 < size0; ++i0) {
                
                    temp1 = std::min(data(i0, i1), data(i0, i1 + 1));
                    temp2 = std::min(data(i0 + 1, i1), data(i0 + 1, i1 + 1));
                    dmin  = std::min(temp1, temp2);
                    temp1 = std::max(data(i0, i1), data(i0, i1 + 1));
                    temp2 = std::max(data(i0 + 1, i1), data(i0 + 1, i1 + 1));
                    dmax  = std::max(temp1, temp2);

                    if (dmax >= levels[0] && dmin <= levels[num_levels - 1]) {
                        for (k = 0; k < num_levels; k++) {
                            if (levels[k] >= dmin && levels[k] <= dmax) {
                                for (m = 4; m >= 0; m--) {
                                    if (m > 0) {
                                        h [m] = data(i0 + OFS0_TBL[m - 1],
                                                     i1 + OFS1_TBL[m - 1]) - levels[k];
                                        xh[m] = i0 + OFS0_TBL[m - 1];
                                        yh[m] = i1 + OFS1_TBL[m - 1];
                                    }
                                    else {
                                        h [0] = 0.25 * (h [1] + h[2] + h[3] + h[4]);
                                        xh[0] = 0.50 * (i0 + i0 + 1);
                                        yh[0] = 0.50 * (i1 + i1 + 1);
                                    }
                                    if (h[m] > 0.0) {
                                        sh[m] =  1;
                                    }
                                    else if (h[m] < 0.0) {
                                        sh[m] = -1;
                                    }
                                    else
                                        sh[m] =  0;
                                }
                                //=================================================================
                                //
                                // Note: at this stage the relative heights of the corners and the
                                // centre are in the h data, and the corresponding coordinates are
                                // in the xh and yh arrays. The centre of the box is indexed by 0
                                // and the 4 corners by 1 to 4 as shown below.
                                // Each triangle is then indexed by the parameter m, and the 3
                                // vertices of each triangle are indexed by parameters m1,m2,and
                                // m3.
                                // It is assumed that the centre of the box is always vertex 2
                                // though this isimportant only when all 3 vertices lie exactly on
                                // the same contour level, in which case only the side of the box
                                // is drawn.
                                //
                                //
                                //      vertex 4 +-------------------+ vertex 3
                                //               | \               / |
                                //               |   \    m-3    /   |
                                //               |     \       /     |
                                //               |       \   /       |
                                //               |  m=2    X   m=2   |       the centre is vertex 0
                                //               |       /   \       |
                                //               |     /       \     |
                                //               |   /    m=1    \   |
                                //               | /               \ |
                                //      vertex 1 +-------------------+ vertex 2
                                //
                                //
                                //
                                //               Scan each triangle in the box
                                //
                                //=================================================================
                                for (m = 1; m <= 4; ++m) {
                                    
                                    m1 = m;
                                    m2 = 0;
                                    m3 = (m != 4) ? m + 1 : 1;
                                    
                                    case_value = CASE_TBL[sh[m1] + 1]
                                                         [sh[m2] + 1]
                                                         [sh[m3] + 1];
                                    
                                    if (case_value != 0) {
                                        switch (case_value) {
                                        case 1:
                                            x1 = xh[m1];
                                            y1 = yh[m1];
                                            x2 = xh[m2];
                                            y2 = yh[m2];
                                        break;
                                        case 2:
                                            x1 = xh[m2];
                                            y1 = yh[m2];
                                            x2 = xh[m3];
                                            y2 = yh[m3];
                                        break;
                                        case 3:
                                            x1 = xh[m3];
                                            y1 = yh[m3];
                                            x2 = xh[m1];
                                            y2 = yh[m1];
                                        break;
                                        case 4:
                                            x1 = xh[m1];
                                            y1 = yh[m1];
                                            x2 = XSECT(m2, m3);
                                            y2 = YSECT(m2, m3);
                                        break;
                                        case 5:
                                            x1 = xh[m2];
                                            y1 = yh[m2];
                                            x2 = XSECT(m3, m1);
                                            y2 = YSECT(m3, m1);
                                        break;
                                        case 6:
                                            x1 = xh[m3];
                                            y1 = yh[m3];
                                            x2 = XSECT(m1, m2);
                                            y2 = YSECT(m1, m2);
                                        break;
                                        case 7:
                                            x1 = XSECT(m1, m2);
                                            y1 = YSECT(m1, m2);
                                            x2 = XSECT(m2, m3);
                                            y2 = YSECT(m2, m3);
                                        break;
                                        case 8:
                                            x1 = XSECT(m2, m3);
                                            y1 = YSECT(m2, m3);
                                            x2 = XSECT(m3, m1);
                                            y2 = YSECT(m3, m1);
                                        break;
                                        case 9:
                                            x1 = XSECT(m3, m1);
                                            y1 = YSECT(m3, m1);
                                            x2 = XSECT(m1, m2);
                                            y2 = YSECT(m1, m2);
                                            break;
                                        default:
                                            break;
                                        }

                                        //
                                        // Overridable: Render the line segment.
                                        //
                                        render(x1, y1, x2, y2, levels[k]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////
        ///  Find a contour of a specified level.
        ///////////////////////////////////////////////////////////////////////
        void do_find(const array_base_type& data, double level)
        {
            double levels[1] = { level };
            do_find(data, levels, 1);
        }

        ///////////////////////////////////////////////////////////////////////
        ///  Render a line segment of the contour.
        ///
        ///  @param x1
        ///  @param y1
        ///  @param x2
        ///  @param y2
        ///  @param level
        ///
        ///////////////////////////////////////////////////////////////////////
        virtual void render(double /* x1    */,
                            double /* y1    */,
                            double /* x2    */,
                            double /* y2    */,
                            double /* level */
                            ) = 0;
    };
    
    } // namespace


///////////////////////////////////////////////////////////////////////////
///  @class contour_png
///  @brief
///////////////////////////////////////////////////////////////////////////
template <typename _Ty, typename _Pixel>
class contour_png : 
    public optnet::detail::contour_base<_Ty>
{
    typedef optnet::detail::contour_base<_Ty> _Base;
public:

    typedef typename _Base::value_type          value_type;
    typedef typename _Base::reference           reference;
    typedef typename _Base::const_reference     const_reference;

    typedef typename _Base::array_base_type     array_base_type;
    typedef typename _Base::array_ref_type      array_ref_type;
    typedef typename _Base::array_type          array_type;

    typedef _Pixel                              pixel_type;

    typedef array2_base<pixel_type>             image_base_type;
    typedef array2_ref<pixel_type>              image_ref_type;
    typedef array2<pixel_type>                  image_type;

    typedef typename _Base::size_type           size_type;

    //////////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    //////////////////////////////////////////////////////////////////////////
    contour_png()
    {
    }

    //////////////////////////////////////////////////////////////////////////
    ///  Find a contour and save the result to a file.
    //////////////////////////////////////////////////////////////////////////
    void find(const char*            outfile,
              const array_base_type& data,
              const image_base_type& image,
              bool                   adjust,
              double                 level = 0,
              double                 scale = 1
              )
    {
        using namespace optnet::io;
        using namespace optnet::utils;

        value_type      min_val, max_val;
        int             x, y, xx, yy;
        unsigned char   gray_value;
        double          delta;

        if (NULL == outfile || scale < 0) {
            throw_exception(std::invalid_argument(
                "contour_png::find: Invalid argument."
                ));
        }

        if (image.size_0() != data.size_0() ||
            image.size_1() != data.size_1())
        {
            throw_exception(std::invalid_argument(
                "contour_png::find: Image sizes do not match."
                ));
        }

        m_size0 = (int)(data.size_0() * scale + 0.5);
        m_size1 = (int)(data.size_1() * scale + 0.5);

        if (m_size0 < 2 || m_size1 < 2) {
            m_size0 = m_size1 = 0;
            throw_exception(std::range_error(
                "contour_png::find: The output image is too small."
                )); 
        }
        else if (m_size0 > 10000 || m_size1 > 10000) {
            m_size0 = m_size1 = 0;
            throw_exception(std::range_error(
                "contour_png::find: The output image is too large."
                )); 
        }

        // Save the scale of the output image.
        //   ... will be used by render().
        m_scale     = scale;
        m_scale_inv = 1.0 / scale;

        std::string filename(outfile);
        if (!str_ends_with_no_case(filename, ".png")) {
            filename += ".png";
        }

        // Create the output image on which to overlay the contour.
        m_image.create(m_size0, m_size1, 1);

        data.min_max_element(&min_val, &max_val);
        delta = (max_val == min_val) ? 1.0 : max_val - min_val;

        for (y = 0; y < m_size1; ++y) {
            
            yy = (int)(y * m_scale_inv + 0.5);
            if (yy >= m_size1) yy = m_size1 - 1;
            else if (yy < 0) yy = 0;
            
            for (x = 0; x < m_size0; ++x) {
                
                xx = (int)(x * m_scale_inv + 0.5);
                if (xx >= m_size0) xx = m_size0 - 1;
                else if (xx < 0) xx = 0;
                
                if (adjust) {
                    gray_value = (unsigned char)
                        ((image(xx, yy) - min_val) / delta * 255);
                }
                else {
                    gray_value = (unsigned char)(image(xx, yy));
                }
                
                m_image(x, y, 0)[0] = gray_value;
                m_image(x, y, 0)[1] = gray_value;
                m_image(x, y, 0)[2] = gray_value;
            } // x
        } // y

        do_find(data, level); // find contour!

        png_save(m_image, filename.c_str());
    }


protected:

    array<rgb>  m_image;
    double      m_scale, m_scale_inv;
    int         m_size0, m_size1;

    //////////////////////////////////////////////////////////////////////////
    inline void put_pixel(int x, int y)
    {
        m_image(x, y, 0)[0] = 255;  //
        m_image(x, y, 0)[1] = 0;    // red
        m_image(x, y, 0)[2] = 0;    //
    }

    //////////////////////////////////////////////////////////////////////////
    void put_line(int x1, int y1, int x2, int y2)
    {
        bool is_longer = false;
        int  dec_inc, long_len = x2 - x1, short_len = y2 - y1;

        if (abs(short_len) > abs(long_len)) {
            int swap  = short_len;
            short_len = long_len;
            long_len  = swap;               
            is_longer = true;
        }

        if (long_len == 0) dec_inc = 0;
        else dec_inc = (short_len << 16) / long_len;

        if (is_longer) {
            if (long_len > 0) {
                long_len += y1;
                for (int j = 0x8000 + (x1 << 16); y1 <= long_len; ++y1) {
                    put_pixel(j >> 16, y1);
                    j += dec_inc;
                }
                return;
            }
            long_len += y1;
            for (int j = 0x8000 + (x1 << 16); y1 >= long_len; --y1) {
                put_pixel(j >> 16, y1); 
                j -= dec_inc;
            }
            return; 
        }

        if (long_len>0) {
            long_len += x1;
            for (int j = 0x8000 + (y1 << 16); x1 <= long_len; ++x1) {
                put_pixel(x1, j >> 16);
                j += dec_inc;
            }
            return;
        }
        long_len += x1;
        for (int j = 0x8000 + (y1 << 16); x1 >= long_len; --x1) {
            put_pixel(x1, j >> 16);
            j-=dec_inc;
        }
    }

    //////////////////////////////////////////////////////////////////////////
    virtual void render(double x1,
                        double y1,
                        double x2,
                        double y2,
                        double level
                        )
    {
        int xx1, yy1, xx2, yy2;

        xx1 = (int)(x1 * m_scale + 0.5);
        yy1 = (int)(y1 * m_scale + 0.5);
        xx2 = (int)(x2 * m_scale + 0.5);
        yy2 = (int)(y2 * m_scale + 0.5);

        if (xx1 >= m_size0) xx1 = m_size0 - 1;
        if (yy1 >= m_size1) yy1 = m_size1 - 1;
        if (xx2 >= m_size0) xx2 = m_size0 - 1;
        if (yy2 >= m_size1) yy2 = m_size1 - 1;

        put_line(xx1, yy1, xx2, yy2);
    }

};

} // namespace

#endif // ___CONTOUR_HXX___
