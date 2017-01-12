/*
 ==========================================================================
 |   
 |   $Id: distance.hxx 159 2005-02-10 02:40:05Z kangli $
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

#ifndef ___DISTANCE_HXX___
#   define ___DISTANCE_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/point2.hxx>
#   include <optnet/_base/point3.hxx>

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
/// Manhattan distance.
#   ifndef __OPTNET_CRAPPY_MSC__
    template <typename _T>
    inline double distance_manhattan(const _T& a, const _T& b)
    {
        return (double)std::abs(a - b);
    }
#   else // VC6
#       define  DISTANCE_MANHATTAN_PROTO(type)                          \
        inline double distance_manhattan(const type& a, const type& b)  \
        {                                                               \
            return (a > b) ? (double)(a - b) : (double)(b - a);         \
        }
    DISTANCE_MANHATTAN_PROTO(float          )
    DISTANCE_MANHATTAN_PROTO(double         )
    DISTANCE_MANHATTAN_PROTO(long double    )
    DISTANCE_MANHATTAN_PROTO(char           )
    DISTANCE_MANHATTAN_PROTO(int            )
    DISTANCE_MANHATTAN_PROTO(long           )
    DISTANCE_MANHATTAN_PROTO(short          )
    DISTANCE_MANHATTAN_PROTO(unsigned char  )
    DISTANCE_MANHATTAN_PROTO(unsigned int   )
    DISTANCE_MANHATTAN_PROTO(unsigned long  )
    DISTANCE_MANHATTAN_PROTO(unsigned short )
#   endif

///////////////////////////////////////////////////////////////////////////
/// Euclidean distance.
template <typename _T, int _Dim>
inline double distance_manhattan(const point<_T, _Dim>& a,
                                 const point<_T, _Dim>& b
                                 )
{
    double sum = 0.0;
    for (int i = 0; i < _Dim; ++i) {
        sum += (a[i] > b[i]) ? 
            (double)(a[i] - b[i]) : (double)(b[i] - a[i]);
    }
    return sum;
}

///////////////////////////////////////////////////////////////////////////
/// Euclidean distance.
#   ifndef __OPTNET_CRAPPY_MSC__
template <typename _T>
    inline double distance_euclidean(const _T& a, const _T& b)
    {
        return sqrt((double)((a - b) * (a - b)));
    }
#   else // VC6
#       define  DISTANCE_EUCLIDEAN_PROTO(type)                          \
        inline double distance_euclidean(const type& a, const type& b)  \
        {                                                               \
            return (double)((a - b) * (a - b));                         \
        }
    DISTANCE_EUCLIDEAN_PROTO(float          )
    DISTANCE_EUCLIDEAN_PROTO(double         )
    DISTANCE_EUCLIDEAN_PROTO(long double    )
    DISTANCE_EUCLIDEAN_PROTO(char           )
    DISTANCE_EUCLIDEAN_PROTO(int            )
    DISTANCE_EUCLIDEAN_PROTO(long           )
    DISTANCE_EUCLIDEAN_PROTO(short          )
    DISTANCE_EUCLIDEAN_PROTO(unsigned char  )
    DISTANCE_EUCLIDEAN_PROTO(unsigned int   )
    DISTANCE_EUCLIDEAN_PROTO(unsigned long  )
    DISTANCE_EUCLIDEAN_PROTO(unsigned short )
#   endif

///////////////////////////////////////////////////////////////////////////
/// Euclidean distance.
template <typename _T, int _Dim>
inline double distance_euclidean(const point<_T, _Dim>& a,
                                 const point<_T, _Dim>& b
                                 )
{
    double sum = 0.0;
    for (int i = 0; i < _Dim; ++i) {
        sum += (double)((a[i] - b[i]) * (a[i] - b[i]));
    }
    return sqrt(sum);
}

///////////////////////////////////////////////////////////////////////////
/// Euclidean distance.
template <typename _T>
inline double distance_euclidean(const point2<_T>& a,
                                 const point2<_T>& b
                                 )
{
    double sum;
    sum = (double)((a[0] - b[0]) * (a[0] - b[0])) +
          (double)((a[1] - b[1]) * (a[1] - b[1]));
    return sqrt(sum);
}

///////////////////////////////////////////////////////////////////////////
/// Euclidean distance.
template <typename _T>
inline double distance_euclidean(const point3<_T>& a,
                                 const point3<_T>& b
                                 )
{
    double sum;
    sum = (double)((a[0] - b[0]) * (a[0] - b[0])) +
          (double)((a[1] - b[1]) * (a[1] - b[1])) +
          (double)((a[2] - b[2]) * (a[2] - b[2]));
    return sqrt(sum);
}

///////////////////////////////////////////////////////////////////////////
/// Squared Euclidean distance.
#   ifndef __OPTNET_CRAPPY_MSC__
    template <typename _T>
    inline double distance_euclidean2(const _T& a, const _T& b)
    {
        return (double)((a - b) * (a - b));
    }
#   else // VC6
#       define  DISTANCE_EUCLIDEAN2_PROTO(type)                         \
        inline double distance_euclidean2(const type& a, const type& b) \
        {                                                               \
            return (double)((a - b) * (a - b));                         \
        }
    DISTANCE_EUCLIDEAN2_PROTO(float         )
    DISTANCE_EUCLIDEAN2_PROTO(double        )
    DISTANCE_EUCLIDEAN2_PROTO(long double   )
    DISTANCE_EUCLIDEAN2_PROTO(char          )
    DISTANCE_EUCLIDEAN2_PROTO(int           )
    DISTANCE_EUCLIDEAN2_PROTO(long          )
    DISTANCE_EUCLIDEAN2_PROTO(short         )
    DISTANCE_EUCLIDEAN2_PROTO(unsigned char )
    DISTANCE_EUCLIDEAN2_PROTO(unsigned int  )
    DISTANCE_EUCLIDEAN2_PROTO(unsigned long )
    DISTANCE_EUCLIDEAN2_PROTO(unsigned short)
#   endif

///////////////////////////////////////////////////////////////////////////
/// Squared Euclidean distance.
template <typename _T, int _Dim>
inline double distance_euclidean2(const point<_T, _Dim>& a,
                                  const point<_T, _Dim>& b
                                  )
{
    double sum = 0.0;
    for (int i = 0; i < _Dim; ++i) {
        sum += (double)((a[i] - b[i]) * (a[i] - b[i]));
    }
    return sum;
}

///////////////////////////////////////////////////////////////////////////
/// Squared Euclidean distance.
template <typename _T>
inline double distance_euclidean2(const point2<_T>& a,
                                  const point2<_T>& b
                                  )
{
    double sum;
    sum = (double)((a[0] - b[0]) * (a[0] - b[0])) +
          (double)((a[1] - b[1]) * (a[1] - b[1]));
    return sum;
}

///////////////////////////////////////////////////////////////////////////
/// Squared Euclidean distance.
template <typename _T>
inline double distance_euclidean2(const point3<_T>& a,
                                  const point3<_T>& b
                                  )
{
    double sum;
    sum = (double)((a[0] - b[0]) * (a[0] - b[0])) +
          (double)((a[1] - b[1]) * (a[1] - b[1])) +
          (double)((a[2] - b[2]) * (a[2] - b[2]));
    return sum;
}


} // namespace

#endif // ___DISTANCE_HXX___
