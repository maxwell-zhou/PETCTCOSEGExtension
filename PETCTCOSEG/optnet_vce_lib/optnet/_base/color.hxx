/*
 ==========================================================================
 |   
 |   $Id: color.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___COLOR_HXX___
#   define ___COLOR_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/_base/point.hxx>

/// @namespace optnet
namespace optnet {

// RGB and RGBA color type
typedef point<unsigned char, 3> rgb;
typedef point<unsigned char, 4> rgba;

// Short integer RGB and RGBA color types
typedef point<short, 3>         short_rgb;
typedef point<short, 4>         short_rgba;

// Float RGB and RGBA color types
typedef point<float, 3>         float_rgb;
typedef point<float, 4>         float_rgba;

// Double RGB and RGBA color types
typedef point<double, 3>        double_rgb;
typedef point<double, 4>        double_rgba;

} // namespace

#endif // ___COLOR_HXX___
