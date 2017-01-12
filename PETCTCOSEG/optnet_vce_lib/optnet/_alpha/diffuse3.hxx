/*
 ==========================================================================
 |   
 |   $Id: diffuse3.hxx 172 2005-02-14 14:18:14Z kangli $
 |
 |   Written by Kang Li <kangl@cmu.edu>
 |   Department of Electrical and Computer Engineering
 |   Carnegie Mellon University
 |   
 ==========================================================================
 |   This file is a part of the OptimalNet library.
 ==========================================================================
 | tmpright (c) 2004-2005 Kang Li <kangl@cmu.edu>. All Rights Reserved.
 | 
 | This software is supplied under the terms of a license agreement or
 | nondisclosure agreement  with the author  and may not be copied  or
 | disclosed except in accordance with the terms of that agreement.
 ==========================================================================
 */

#ifndef ___DIFFUSE3_HXX___
#   define ___DIFFUSE3_HXX___

#   include <optnet/_base/array.hxx>
#   include <optnet/_base/except.hxx>
#   include <limits>
#   include <cmath>
    
#   ifdef max
#       undef max
#   endif

#   ifdef min
#       undef min
#   endif

#   define MIRROR_L(i, n) ((i) < 0) ? -(i) : (i)
#   define MIRROR_R(i, n) ((i) >= (n)) ? ((n) - ((i) - (n)) - 2) : (i)
    
namespace optnet {
    
template <typename _VoxelT, typename _weights, typename _Tg>
bool diffuse3(array_base<_VoxelT, _Tg>&        image,
              const array_base<_weights, _Tg>& weights,
              const double&                    amount,
              const double&                    step_size,
              const double&                    stop_tol,
              int                              max_iter
              )
{
    typedef double real_type;
    
    int     iter, src, dst;
    int     x, y, z, sx, sy, sz, xx, yy, zz;
    double  term1, term2, term3, temp, change;
    double  alpha = step_size / (amount * amount);
    double  denom = 1.0 / (1.0 + alpha);
    bool    converged = false;
    
    sx = (int)image.size_0();
    sy = (int)image.size_1();
    sz = (int)image.size_2();
    
    if (sx != (int)weights.size_0() ||
        sy != (int)weights.size_1() ||
        sz != (int)weights.size_2())
    {
        throw_exception(std::runtime_error("Size mismatch."));
    }
    
    double per = 1.0 / (sx * sy * sz);
    
    array<real_type> tmp(sx, sy, sz, 2);
    
    for (z = 0; z < sz; ++z) {
        for (y = 0; y < sy; ++y) {
            for (x = 0; x < sx; ++x) {
                tmp(x, y, z, 0) = image(x, y, z);
            } // x
        } // y
    } // z
    
    src = 1;
    dst = 0;
    
    for (iter = 0; iter < max_iter; ++iter) {
        
        change = .0;

        src = !src;
        dst = !dst;
            
        // Loop through all voxels.
        for (z = 0; z < sz; ++z) {
            for (y = 0; y < sy; ++y) {
                for (x = 0; x < sx; ++x) {
        
                    term1 = term2 = 0;
                    
                    /* 6-neighbor system */
                    xx = MIRROR_L(x - 1, sx);       
                    temp = sqrt((double)(weights(xx, y, z) * weights(x, y, z)));
                    term1 += temp;
                    term2 += temp * tmp(xx, y, z, src);
                    /*---*/
                    xx = MIRROR_R(x + 1, sx);
                    temp = sqrt((double)(weights(xx, y, z) * weights(x, y, z)));
                    term1 += temp;
                    term2 += temp * tmp(xx, y, z, src);
                    /*---*/
                    yy = MIRROR_L(y - 1, sy);
                    temp = sqrt((double)(weights(x, yy, z) * weights(x, y, z)));
                    term1 += temp;
                    term2 += temp * tmp(x, yy, z, src);
                    /*---*/
                    yy = MIRROR_R(y + 1, sy);
                    temp = sqrt((double)(weights(x, yy, z) * weights(x, y, z)));
                    term1 += temp;
                    term2 += temp * tmp(x, yy, z, src);
                    /*---*/
                    zz = MIRROR_L(z - 1, sz);
                    temp = sqrt((double)(weights(x, y, zz) * weights(x, y, z)));
                    term1 += temp;
                    term2 += temp * tmp(x, y, zz, src);
                    /*---*/
                    zz = MIRROR_R(z + 1, sz);
                    temp = sqrt((double)(weights(x, y, zz) * weights(x, y, z)));
                    term1 += temp;
                    term2 += temp * tmp(x, y, zz, src);
                    
                    term1 = (1 - step_size * term1) * tmp(x, y, z, src);
                    term2 = step_size * term2;
                    term3 = alpha * image(x, y, z);
                    
                    tmp(x, y, z, dst) = (term1 + term2 + term3) * denom;
                
                    change += fabs(tmp(x, y, z, dst) - 
                                   tmp(x, y, z, src));
                                                        
                } // x
            } // y
        } // z

        if ((change * per) <= stop_tol) {
            converged = true;
            break;
        }
    } // iter
    
    _VoxelT vmax = std::numeric_limits<_VoxelT>::max();
    _VoxelT vmin = std::numeric_limits<_VoxelT>::min();
    
    for (z = 0; z < sz; ++z) {
        for (y = 0; y < sy; ++y) {
            for (x = 0; x < sx; ++x) {
                real_type v = tmp(x, y, z, dst);
                if (v > vmax) v = vmax;
                else if (v < vmin) v = vmin;
                image(x, y, z) = static_cast<_VoxelT>(v);
            } // x
        } // y
    } // z
    
    return converged;
}

} // namespace

#   undef MIRROR_L
#   undef MIRROR_R

#endif
