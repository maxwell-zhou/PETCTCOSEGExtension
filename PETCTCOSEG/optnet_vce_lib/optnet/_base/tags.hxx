/*
 ==========================================================================
 |   
 |   $Id: tags.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___TAGS_HXX___
#   define ___TAGS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class net_f_xy
///  @brief Surface direction tag.
///
///  This indicates that the surface can be expressed as z = f(x, y).
///////////////////////////////////////////////////////////////////////////
class net_f_xy {};

///////////////////////////////////////////////////////////////////////////
///  @class net_f_yz
///  @brief Surface direction tag.
///
///  This indicates that the surface can be expressed as x = f(y, z).
///////////////////////////////////////////////////////////////////////////
class net_f_yz {};

///////////////////////////////////////////////////////////////////////////
///  @class net_f_zx
///  @brief Surface direction tag.
///
///  This indicates that the surface can be expressed as y = f(z, x).
///////////////////////////////////////////////////////////////////////////
class net_f_zx {};

///////////////////////////////////////////////////////////////////////////
///  @class net_f_xy_zflipped
///  @brief Surface direction tag.
///
///  This indicates that the surface can be written as z = Z - 1 - f(x, y).
///////////////////////////////////////////////////////////////////////////
class net_f_xy_zflipped {};

///////////////////////////////////////////////////////////////////////////
///  @class net_f_yz_xflipped
///  @brief Surface direction tag.
///
///  This indicates that the surface can be written as x = X - 1 - f(y, z).
///////////////////////////////////////////////////////////////////////////
class net_f_yz_xflipped {};

///////////////////////////////////////////////////////////////////////////
///  @class net_f_zx_yflipped
///  @brief Surface direction tag.
///
///  This indicates that the surface can be written as y = Y - 1 - f(z, x).
///////////////////////////////////////////////////////////////////////////
class net_f_zx_yflipped {};

} // namespace


#endif 
