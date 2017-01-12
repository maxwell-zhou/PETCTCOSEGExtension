/*
 ==========================================================================
 |   
 |   $Id: io_options.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___IO_OPTIONS__
#   define ___IO_OPTIONS__

#   include <string>

namespace optnet {
    namespace io {

///////////////////////////////////////////////////////////////////////////
///  @class save_options
///  @brief Save options.
///////////////////////////////////////////////////////////////////////////
struct save_options
{
    //  Options:
    bool        single_file;        /// Save the image as a single entity,
                                    ///   if possible.
    int         compression_level;  /// Compression level, 0 if no compression.
    int         slice_to_save;      /// Index of the slice needs to save.
    int         phase_to_save;      /// Index of the phase needs to save.
};

    } // namespace
} // namespace

#endif
