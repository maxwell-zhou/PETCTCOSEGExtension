/*
 ==========================================================================
 |   
 |   $Id: consts.hxx 186 2005-02-27 02:04:22Z kangli $
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

#ifndef ___CONSTS_HXX___
#   define ___CONSTS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

namespace optnet {

static const double INV3 = 1.0 / 3.0; /// Reciprocal of 3.
static const double INV6 = 1.0 / 6.0; /// Reciprocal of 6.
static const double INV7 = 1.0 / 7.0; /// Reciprocal of 7.
static const double INV9 = 1.0 / 9.0; /// Reciprocal of 9.

} // namespace

#endif
