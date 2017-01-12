/*
 ==========================================================================
 |   
 |   $Id: optnet_sys.hxx 146 2005-02-08 05:26:28Z kangli $
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


#ifndef ___OPTNET_SYS_HXX___
#   define ___OPTNET_SYS_HXX___

/** \file optnet_sys.hxx
 *  The header file containing the C++ interfaces of the system-related
 *  classes and functions provided by the library.
 */

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/config.h>  // Configurations

#   include <optnet/_sys/getopt.hxx>
#   include <optnet/_sys/process_info.hxx>
    
#   ifdef __OPTNET_HAS_GLUT__
#       include <optnet/_sys/gui/app_glut.hxx>
#   endif
    
#endif
