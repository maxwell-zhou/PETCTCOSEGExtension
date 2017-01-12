/*
 ==========================================================================
 |   
 |   $Id: errors.h 177 2005-02-17 20:26:01Z kangli $
 |
 |   OptimalNet Library Standard C Interface Header
 |
 |   Written by Kang Li <likang@ieee.org>
 |   Department of Electrical and Computer Engineering
 |   Carnegie Mellon University
 |   
 ==========================================================================
 |   This file is a part of the OptimalNet library.
 ==========================================================================
 | Copyright (c) 2004 Kang Li <likang@ieee.org>. All Rights Reserved.
 | 
 | This software is supplied under the terms of a license agreement or
 | nondisclosure agreement  with the author  and may not be copied  or
 | disclosed except in accordance with the terms of that agreement.
 ==========================================================================
 */


#ifndef ___ERRORS_H___
#   define ___ERRORS_H___

/** \file errors.h
 *  The header file containing error codes for the library.
 */

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

/* ===================================================================== */
/*  Error Codes                                                          */
/* ===================================================================== */

/** The operation completed successfully. */
#   define OPTNET_S_OK              0x00000000L

/** An undetermined error occurred inside the library. */
#   define OPTNET_E_FAILED          0xFFFFFFFFL

/** An invalid parameter was passed to the returning function. */
#   define OPTNET_E_INVALID_ARG     0xF0000001L

/** The funtion is not implemented and the task can not be fulfilled. */
#   define OPTNET_E_UNIMPLEMENTED   0xF0000002L

/** The program could not allocate sufficient memory to complete the call. */
#   define OPTNET_E_OUT_OF_MEMORY   0xF0000003L

/** A runtime error occurred inside the library. */
#   define OPTNET_E_RUNTIME_ERROR   0xF0000004L

/** An I/O error occured inside the library. */
#   define OPTNET_E_IO_ERROR        0xF0000005L


#endif
