/*
 ==========================================================================
 |   
 |   $Id: optnet_io.hxx 21 2005-01-14 15:52:31Z kangli $
 |
 |   OptimalNet Library I/O Facility Standard C++ Interface Header
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


#ifndef ___OPTNET_IO_HXX___
#   define ___OPTNET_IO_HXX___

/** \file optnet_io.hxx
 *  The header file containing the C++ interfaces of the data I/O
 *  facilities provided by the library.
 */

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/config.h>  // Configurations

#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>

#   include <optnet/_base/io/analyze.hxx>
#   include <optnet/_base/io/png.hxx>

#endif
