/*
 ==========================================================================
 |   
 |   $Id: optnet.h 21 2005-01-14 15:52:31Z kangli $
 |
 |   OptimalNet Library Standard C Interface Header
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


#ifndef ___OPTNET_H___
#   define ___OPTNET_H___

/** \file optnet.h
 *  The header file containing the C interface functions and definitions.
 */

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/config.h>   /* Configurations */
#   include <optnet/define.h>   /* Type Definitions */
#   include <optnet/errors.h>   /* Error Codes */

/* ===================================================================== */
/*  Algorithm Options                                                    */
/* ===================================================================== */

#   define OPTNET_M_IA              0   /* Implicit-Arc  */
#   define OPTNET_M_FS              1   /* Forward-Star  */
#   define OPTNET_M_PR              2   /* Push-Relabel  */

#   define OPTNET_F_XY              0   /* Net = f(x, y) */
#   define OPTNET_F_YZ              1   /* Net = f(y, z) */
#   define OPTNET_F_ZX              2   /* Net = f(z, x) */


/* ===================================================================== */
/*  Function Prototypes                                                  */
/* ===================================================================== */

#   ifdef __cplusplus
extern "C" {
#   endif

long OPTNET_API optnet_get_version();

long OPTNET_API optnet_3d_1b(int*                    net,
                             const unsigned char*    cost,
                             unsigned int            s0,
                             unsigned int            s1,
                             unsigned int            s2,
                             int                     m0,
                             int                     m1,
                             int                     c0,
                             int                     c1,
                             int                     orient,
                             int                     method
                             );

long OPTNET_API optnet_3d_1d(int*                    net,
                             const double*           cost,
                             unsigned int            s0,
                             unsigned int            s1,
                             unsigned int            s2,
                             int                     m0,
                             int                     m1,
                             int                     c0,
                             int                     c1,
                             int                     orient,
                             int                     method
                             );

long OPTNET_API optnet_3d_1f(int*                    net,
                             const float*            cost,
                             unsigned int            s0,
                             unsigned int            s1,
                             unsigned int            s2,
                             int                     m0,
                             int                     m1,
                             int                     c0,
                             int                     c1,
                             int                     orient,
                             int                     method
                             );

long OPTNET_API optnet_3d_1s(int*                    net,
                             const short*            cost,
                             unsigned int            s0,
                             unsigned int            s1,
                             unsigned int            s2,
                             int                     m0,
                             int                     m1,
                             int                     c0,
                             int                     c1,
                             int                     orient,
                             int                     method
                             );

#   ifdef __cplusplus
}
#   endif

#   include <optnet/version.h> /* Version */

#endif
