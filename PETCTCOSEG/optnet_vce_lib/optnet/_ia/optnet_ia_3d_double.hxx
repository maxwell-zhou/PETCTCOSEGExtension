/*
 ==========================================================================
 |   
 |   $Id: optnet_ia_3d_double.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___OPTNET_IA_3D_DOUBLE_HXX___
#   define ___OPTNET_IA_3D_DOUBLE_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

namespace optnet {

//TODO: Implement the implicit-arc version of the double-surface
//      Optimal-Net algorithm here.

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_ia/optnet_ia_3d_double.cxx>
#   endif

#endif
