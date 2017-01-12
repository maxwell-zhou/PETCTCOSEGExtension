/*
 ==========================================================================
 |   
 |   $Id: optnet_utils.hxx 126 2005-02-05 18:06:44Z kangli $
 |
 |   OptimalNet Library Standard C++ Interface Header
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


#ifndef ___OPTNET_UTILS_HXX___
#   define ___OPTNET_UTILS_HXX___

/** \file optnet_utils.hxx
 *  The header file containing the C++ interfaces of the utility classes
 *  and functions provided by the library.
 */

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/config.h>  // Configurations

#   include <optnet/_utils/arith.hxx>
#   include <optnet/_utils/endian.hxx>
#   include <optnet/_utils/filename.hxx>
#   include <optnet/_utils/gaussian.hxx>
#   include <optnet/_utils/hessian.hxx>
#   include <optnet/_utils/index.hxx>
#   include <optnet/_utils/interp_bspline2.hxx>
#   include <optnet/_utils/interp_bspline3.hxx>
#   include <optnet/_utils/regiongrow.hxx>
#   include <optnet/_utils/statistics.hxx>
#   include <optnet/_utils/timer.hxx>
#   include <optnet/_utils/xstring.hxx>

#endif
