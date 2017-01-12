/*
 ==========================================================================
 |   
 |   $Id: datatype.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___DATATYPE_HXX___
#   define ___DATATYPE_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <cstddef>
#   include <optnet/_base/color.hxx>

/// @namespace optnet
namespace optnet {

/// @enum datatype_id
enum datatype_id
{
    OPTNET_DT_UNKNOWN           = 0x00,
    // integer datatypes
    OPTNET_DT_CHAR              = 0x01,
    OPTNET_DT_UNSIGNED_CHAR     = 0x02,
    OPTNET_DT_SHORT             = 0x03,
    OPTNET_DT_UNSIGNED_SHORT    = 0x04,
    OPTNET_DT_INT               = 0x05,
    OPTNET_DT_UNSIGNED_INT      = 0x06,
    OPTNET_DT_LONG              = 0x07,
    OPTNET_DT_UNSIGNED_LONG     = 0x08,
    OPTNET_DT_INT64             = 0x09, //
    OPTNET_DT_UNSIGNED_INT64    = 0x0A, //
    OPTNET_DT_LONGLONG          = 0x0B, // not formally supported
    OPTNET_DT_UNSIGNED_LONGLONG = 0x0C, //
    // floating point datatypes
    OPTNET_DT_FLOAT             = 0x11,
    OPTNET_DT_DOUBLE            = 0x12,
    OPTNET_DT_LONG_DOUBLE       = 0x13,
    // color datatypes
    OPTNET_DT_RGB               = 0x21,
    OPTNET_DT_FLOAT_RGB         = 0x22,
    OPTNET_DT_SHORT_RGB         = 0x23,
    OPTNET_DT_RGBA              = 0x31,
    OPTNET_DT_FLOAT_RGBA        = 0x32,
    OPTNET_DT_SHORT_RGBA        = 0x33
};

///////////////////////////////////////////////////////////////////////////
///  A template class to obtain the numerical id of datatypes.
///////////////////////////////////////////////////////////////////////////
template <typename _T>
struct datatype // This should never be called.
{ static long id() { assert(0); return OPTNET_DT_UNKNOWN; }};

template <> struct datatype<char>
{ static long id() { return OPTNET_DT_CHAR; }};

template <> struct datatype<unsigned char>
{ static long id() { return OPTNET_DT_UNSIGNED_CHAR; }};

template <> struct datatype<short>
{ static long id() { return OPTNET_DT_SHORT; }};

template <> struct datatype<unsigned short>
{ static long id() { return OPTNET_DT_UNSIGNED_SHORT; }};

template <> struct datatype<int>
{ static long id() { return OPTNET_DT_INT; }};

template <> struct datatype<unsigned int>
{ static long id() { return OPTNET_DT_UNSIGNED_INT; }};

template <> struct datatype<long>
{ static long id() { return OPTNET_DT_LONG; }};

template <> struct datatype<unsigned long>
{ static long id() { return OPTNET_DT_UNSIGNED_LONG; }};

template <> struct datatype<float>
{ static long id() { return OPTNET_DT_FLOAT; }};

template <> struct datatype<double>
{ static long id() { return OPTNET_DT_DOUBLE; }};

template <> struct datatype<long double>
{ static long id() { return OPTNET_DT_LONG_DOUBLE; }};

//
// color datatypes
//
template <> struct datatype<rgb>
{ static long id() { return OPTNET_DT_RGB; }};

template <> struct datatype<float_rgb>
{ static long id() { return OPTNET_DT_FLOAT_RGB; }};

template <> struct datatype<short_rgb>
{ static long id() { return OPTNET_DT_SHORT_RGB; }};

template <> struct datatype<rgba>
{ static long id() { return OPTNET_DT_RGBA; }};

template <> struct datatype<float_rgba>
{ static long id() { return OPTNET_DT_FLOAT_RGBA; }};

template <> struct datatype<short_rgba>
{ static long id() { return OPTNET_DT_SHORT_RGBA; }};

} // namespace

#endif 
