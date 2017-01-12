/*
 ==========================================================================
 |   
 |   $Id: config.h 113 2005-02-03 07:14:47Z kangli $
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

#ifndef ___CONFIG_H___
#   define ___CONFIG_H___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   if (defined(__unix__) || defined(unix)) && !defined(USG)
#       include <sys/param.h>
#   endif

/* operating system recognition */
#   if defined(linux) || defined(__linux) || defined(__linux__)
        /* Linux */
#       define __OPTNET_OS_LINUX__
#   elif !defined(SAG_COM) \
    && (defined(WIN64) || defined(_WIN64) || defined(__WIN64__))
        /* Windows */
#       define __OPTNET_OS_WIN32__
#       define __OPTNET_OS_WIN64__
#       define __OPTNET_OS_WINNT__
#   elif !defined(SAG_COM) \
    && ( \
        defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined \
            (__NT__) \
        )
#       define __OPTNET_OS_WIN32__
#       define __OPTNET_OS_WINNT__
#   elif defined(__hpux)
        /* HP Unix */
#       define __OPTNET_OS_HPUX__
#   elif defined(BSD) || defined(__FreeBSD__)
        /* BSD */
#       define __OPTNET_OS_BSD__
#   elif defined(_MAC) || defined(__APPLE__)
        /* Mac */
#       define __OPTNET_OS_MAC__
#   endif

/* compiler specific configuration */
#   if defined(__GNUC__)
#       define __OPTNET_CC_GNUC__
#       define __OPTNET_CC_GNUC_VER__
#   elif defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ECL)
#       define __OPTNET_CC_INTEL__
#       define __OPTNET_CC_INTEL_VER__
#   elif defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#       define __OPTNET_CC_MSC__
#       define __OPTNET_CC_MSC_VER__   _MSC_VER
#   elif defined(__LCC__)
#       define __OPTNET_CC_LCC__
#       define __OPTNET_CC_LCC_VER__
#   else
#       define __OPTNET_CC_UNKNOWN__
#   endif

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200) // VC 6.0
#       define __OPTNET_CRAPPY_MSC__
#   endif

#   if defined(_MSC_VER) && (_MSC_VER >= 1400) // VC 8.0
#       define __OPTNET_SECURE_STR__
#       if defined(_OPENMP)
#           define __OPTNET_PRAGMA_OMP__
#           define __OPTNET_OMP_NUM_THREADS__ 4
#       endif
#   endif

#   ifndef OPTNET_IMPEXP
#      ifdef OPTNET_EXPORTS
#         define OPTNET_IMPEXP __declspec(dllexport)
#      else
#         define OPTNET_IMPEXP __declspec(dllimport)
#      endif
#   endif

/* API calling convention */
#   if defined(__OPTNET_CC_MSC__) || defined(__OPTNET_CC_INTEL__) \
        || defined(__OPTNET_CC_LCC__)
#       define OPTNET_API /* OPTNET_IMPEXP */
#   else
#       define OPTNET_API
#   endif

/* debugging facilities */
#   if defined(_DEBUG) || defined(DEBUG) && !defined(NDEBUG) \
    || (defined(__LCC__) && __LCCDEBUGLEVEL >= 2)
#       define __OPTNET_DEBUG__ 1
#   endif

/* compiler features */
#   if defined(__OPTNET_CC_MSC__)
#       if (__OPTNET_CC_MSC_VER__ > 1200)
#           define __OPTNET_MEMBER_TEMPLATES__
#       endif
#   else
#       define __OPTNET_MEMBER_TEMPLATES__
#       define __OPTNET_REMOVE_UNUSED_ARG__
#   endif

#   if defined(__OPTNET_REMOVE_UNUSED_ARG__)
#       define OPTNET_UNUSED(arg) /* arg */
#   else  // stupid, broken compiler
#       define OPTNET_UNUSED(arg) arg
#   endif

/* library features */
#   define __OPTNET_SUPPORT_ROI__

#   if defined(__OPTNET_OS_WINNT__) && \
            (defined(__OPTNET_CC_MSC__) || defined(__OPTNET_CC_INTEL__))
#       define __OPTNET_USE_MSXML__ // use MSXML
#   endif

#endif
