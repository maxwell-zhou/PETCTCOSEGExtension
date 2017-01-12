/*
 ==========================================================================
 |   
 |   $Id: optnet_xtra.h 21 2005-01-14 15:52:31Z kangli $
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


#ifndef ___OPTNET_XTRA_H___
#   define ___OPTNET_XTRA_H___

/** \file optnet_xtra.h
 *  The header file containing the C interface functions and definitions.
 */

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/config.h>   /* Configurations */
#   include <optnet/errors.h>   /* Error Codes */


/* ===================================================================== */
/*  Function Prototypes                                                  */
/* ===================================================================== */

#   ifdef __cplusplus
extern "C" {
#   endif

long OPTNET_API optnet_igc_bk_b  (unsigned char*       mask,
                                  const unsigned char* image,
                                  const unsigned char* trimap,
                                  int                  s0,
                                  int                  s1,
                                  int                  s2,
                                  int                  nh,
                                  double               lambda,
                                  double               beta
                                  );

long OPTNET_API optnet_igc_bk_d  (unsigned char*       mask,
                                  const double*        image,
                                  const unsigned char* trimap,
                                  int                  s0,
                                  int                  s1,
                                  int                  s2,
                                  int                  nh,
                                  double               lambda,
                                  double               beta
                                  );

long OPTNET_API optnet_igc_bk_rgb(unsigned char*       mask,
                                  const unsigned char* image,
                                  const unsigned char* trimap,
                                  int                  s0,
                                  int                  s1,
                                  int                  s2,
                                  int                  nh,
                                  double               lambda,
                                  double               beta
                                  );

#   ifdef __cplusplus
}
#   endif

#endif
