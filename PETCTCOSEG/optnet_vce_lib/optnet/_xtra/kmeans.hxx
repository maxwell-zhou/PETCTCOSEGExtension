/*
 ==========================================================================
 |   
 |   $Id: kmeans.hxx 186 2005-02-27 02:04:22Z kangli $
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
  - Purpose:
     
     This file implements the fast k-means algorithm.
  
  - Reference:
     
     Charles Elkan (ELKAN@CS.UCSD.EDU)
     Using the Triangle Inequality to Accelerate k-Means
     Proceedings of the Twentieth International Conference on Machine
         Learning (ICML-2003), Washington DC, 2003.
 ==========================================================================
 */

#ifndef ___KMEANS_HXX___
#   define ___KMEANS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <cstdlib>
#   include <optnet/_base/iterator.hxx>
#   include <optnet/_base/memory.hxx>
#   include <optnet/_base/deref.hxx>


/// @namespace optnet
namespace optnet {
    /// @namespace xtra
    namespace xtra {

///////////////////////////////////////////////////////////////////////////
///  K-means clustering.
///////////////////////////////////////////////////////////////////////////
template <typename _IIt, typename _OIt, typename _Dist> inline
int kmeans(_IIt     input_begin,
           _IIt     input_end,
           _OIt     means_begin,
           _OIt     means_end,
           int*     cids,
           _Dist    dist
           )
{
    return kmeans(input_begin, input_end, means_begin, means_end, cids,
                  dist, make_deref(input_begin));
} // kmeans


///////////////////////////////////////////////////////////////////////////
///  K-means clustering.
///////////////////////////////////////////////////////////////////////////
template <typename _IIt, typename _OIt, typename _Dist, typename _Deref>
int kmeans(_IIt     input_begin,
           _IIt     input_end,
           _OIt     means_begin,
           _OIt     means_end,
           int*     cids,
           _Dist    dist,
           _Deref   deref
           )
{
    using optnet::xtra::detail;
    typedef typename _OIt::value_type ovalue_type;
    typedef typename _IIt::value_type ivalue_type;

    // number of clusters
    int kused = 0;
    int k = (int)abs(&*means_end - &*means_begin);
    int n = (int)abs(&*input_end - &*input_begin);
    
    assert(k != 0);
    assert(n != 0);

    // inter-centroid distances, lower bounds and upper bounds
    double** cd = new_matrix<double>(k, k, 0.0);
    double** lb = new_matrix<double>(k, n, 0.0);
    double*  ub = new double[n];

    //
    // initialization
    //
    {
        _IIt        iit = input_begin;
        _OIt        oit2, oit = means_begin;
        ivalue_type c, c0;
        int         i, j, k;
        
        c0 = deref(iit); ++iit;
        for (j = 1; iit != input_end; ++iit, ++j)
            c0 += deref(iit); // sum of input data
        c0 /= j; // average

        // use the mean of all input data as the first center
        *oit++ = (ovalue_type)c0;
        for (iit = input_begin, j = 0; iit != input_end; ++iit, ++j) {
            lb[0][j] = dist(deref(iit), c0);
            ub[j]    = lb[j][0];
            cids[j]  = 0;
        }
        
        // find additional centers (furthest first)
        for (i = 1; oit != means_end; ++oit, ++i) {
            double d;
            
            k = 0;
            d = ub[0];
            for (j = 1; j < n; ++j) {
                if (ub[j] > d) {
                    d = ub[j];
                    k = j;
                }
            }

            *oit = c = deref(input_begin + k);
            
            for (oit2 = means_begin, j = 0; oit2 != oit; ++oit2, ++j)
                cd[i][j] = dist(*oit2, c); // inter-center distances

            for (iit  = input_begin, j = 0;
                 iit != input_end;
                 ++iit, ++j) {

                if (2 * ub[j] > cd[i][j]) {
                    lb[i][j] = d = dist(deref(iit), c);
                    if (d < ub[j]) {
                        ub[j]   = d;
                        cids[j] = j;
                    }
                }
            }

        } // for i
    }

    //
    // begin k-means clustering
    //

    // TODO: Implement this.

    // release all resources
    delete_matrix(cd);
    delete_matrix(lb);
    delete[] ub;

    return kused;
} // kmeans


    } // namespace
} // namespace

#endif // ___KMEANS_HXX___
