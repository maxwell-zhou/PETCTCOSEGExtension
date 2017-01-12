/*
 ==========================================================================
 |   
 |   $Id: optnet.hxx 156 2005-02-09 23:35:57Z kangli $
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


#ifndef ___OPTNET_HXX___
#   define ___OPTNET_HXX___

/** \file optnet.hxx
 *  The header file containing the C++ interfaces of the main algorithms
 *  provided by the library.
 */

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/config.h>  // Configurations

#   include <optnet/_base/comp.hxx>
#   include <optnet/_base/debug.hxx>
#   include <optnet/_base/point2.hxx>
#   include <optnet/_base/point3.hxx>

#   ifndef __OPTNET_NOALGO__
#       include <optnet/_ia/optnet_ia_3d.hxx>
#       include <optnet/_ia/optnet_ia_3d_double.hxx>

#       include <optnet/_fs/optnet_fs_3d.hxx>
#       include <optnet/_fs/optnet_fs_3d_multi.hxx>

#       include <optnet/_pr/optnet_pr_3d.hxx>
#       include <optnet/_pr/optnet_pr_3d_multi.hxx>
#   endif // __OPTNET_BASICS__

#endif
