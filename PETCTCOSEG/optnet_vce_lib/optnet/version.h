/*
 ==========================================================================
 |   
 |   $Id: version.h 194 2005-02-27 04:04:12Z kangli $
 |
 |   OptimalNet Library Platform-Specific Configurations
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

#ifndef ___VERSION_H___
#   define ___VERSION_H___

/** \file version.h
 *  The header file containing the version number of the library.
 */

/**
 *  \def OPTNET_VERSION
 *  Defines the version number of the library.
 *    - OPTNET_VERSION / 100000      is the major version
 *    - OPTNET_VERSION / 100 \% 1000 is the minor version
 *    - OPTNET_VERSION \% 100        is the sub-minor version
 */

#   define OPTNET_VERSION 100100

#endif
