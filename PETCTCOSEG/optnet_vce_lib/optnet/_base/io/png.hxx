/*
 ==========================================================================
 |   
 |   $Id: png.hxx 58 2005-01-21 01:44:10Z kangli $
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

#ifndef ___PNG_HXX__
#   define ___PNG_HXX__

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/_base/except.hxx>
#   include <optnet/_base/io/png_core.hxx>
#   include <optnet/_base/io/io_options.hxx>
#   include <optnet/_utils/filename.hxx>

/// @namespace optnet
namespace optnet {
    
    /// @namespace optnet::io
    namespace io {
    
///////////////////////////////////////////////////////////////////////////
///  Loads a PNG image file into memory.
///
///  @param[out] a      The a for storing the image data.
///  @param[in]  name   The name of the image file.
///
///  @exception  optnet::io_error Failed reading image header or data, or
///                               the data type of a does not match
///                               that used in the file.
///
///////////////////////////////////////////////////////////////////////////
template <typename _Array>
void
png_load(_Array& a, const char* name)
{
    using namespace optnet::io::detail;

    FILE*               fp;
    optnet_png_slice    slice;
    unsigned int        sig_read;
    long                dt = OPTNET_DT_UNKNOWN;
    int                 i;

    // Get a reference to the extension structure of a.
    typename _Array::extension_type& ext = a.extension();

    // Infer datatype of a.
#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200) // VC 6.0
    dt = datatype</*typename*/ _Array::value_type>::id();
#   else
    dt = datatype<  typename   _Array::value_type>::id();
#   endif

    fp = fopen(name, "rb");
    if (NULL == fp) {
        throw_exception(io_error(
            "png_load: Failed opening file"
        ));
    }

    if (!png_check_format(fp, &sig_read)) {
        fclose(fp);
        throw_exception(io_error(
            "png_load: Invalid file format"
        ));
    }

    if (!png_load_info(fp, &slice, sig_read, dt)) {
        fclose(fp);
        std::string msg;
        switch (slice.errcode) {
        case OPTNET_PNG_E_FORMAT:
            msg = "png_load: Invalid file format";
            break;
        case OPTNET_PNG_E_READ_STRUCT:
            msg = "png_load: Failed png_create_read_struct";
            break;
        case OPTNET_PNG_E_INFO_STRUCT:
            msg = "png_load: Failed png_create_info_struct";
            break;
        case OPTNET_PNG_E_SETJMP:
            msg = "png_load: Failed reading image.";
            break;
        case OPTNET_PNG_E_SIZE:
            msg = "png_load: Incorrect image size";
            break;
        default:
            msg = "png_load: Unknown error";
        }
        throw_exception(io_error(msg));
    }

    // Allocate storage for the image data.
    a.create(slice.width, slice.height, 1);
    slice.data = reinterpret_cast<unsigned char*>(a.data());
    if (NULL == slice.data) {
        fclose(fp);
        throw_exception(io_error("png_load: Not enough memory"));
    }

    ext.voxel_size[0] = (slice.x_ppm != 0) ?
        1.0 / slice.x_ppm : 0;
    ext.voxel_size[1] = (slice.y_ppm != 0) ?
        1.0 / slice.y_ppm : 0;
    ext.voxel_size[2] = 0;
    ext.voxel_size[3] = 0;

    // Save textual information if present.
    if (NULL != slice.title)        ext.title       = slice.title;
    if (NULL != slice.author)       ext.author      = slice.author;
    if (NULL != slice.description)  ext.description = slice.description;
    if (NULL != slice.copyright)    ext.copyright   = slice.copyright;

    // Load image data into memory.
    if (!png_load_data(&slice)) {
        fclose(fp);
        throw_exception(io_error(
            "png_load: Failed loading image data"
            ));
    }

    // Retrieve and store palette.
    if ((slice.has_palette) && (slice.num_palette > 0)) {
        // Get a reference to the palette of a.
        std::vector<rgb>& palette = ext.palette;
        // Copy the palette.
        palette.resize(slice.num_palette);
        for (i = 0; i < slice.num_palette; ++i) {
            palette[i][0] = slice.palette[i].red;
            palette[i][1] = slice.palette[i].green;
            palette[i][2] = slice.palette[i].blue;
        } // for
    } // if

    // Cleaning up.
    png_load_free(&slice);
    fclose(fp);

} // png_load

///////////////////////////////////////////////////////////////////////////
///  Saves an image into a PNG format image file.
///
///  @param[in] a        The a for storing the image data.
///  @param[in] name     The name of the image file.
///  @param[in] options  Options.
///
///  @exception optnet::io_error Failed writing image header or data.
///
///////////////////////////////////////////////////////////////////////////
template <typename _Array>
void
png_save(_Array&             a,
         const char*         name,
         const save_options* options = NULL
         )
{
    using namespace optnet::io::detail;
    using namespace optnet::utils;

    if (a.size() == 0) return; // Nothing to write.

    FILE*               fp;
    optnet_png_slice    slice;
    size_t              offset;
    size_t              first_slice;
    size_t              slices_to_save;
    long                dt = OPTNET_DT_UNKNOWN;

    // Get a reference to the extension structure of a.
    const typename _Array::extension_type&  ext   = a.extension();
    const typename _Array::size_type*       sizes = a.sizes();

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200) // VC 6.0
    dt = datatype</*typename*/ _Array::value_type>::id();
#   else
    dt = datatype<  typename   _Array::value_type>::id();
#   endif

    memset(&slice, 0, sizeof(optnet_png_slice));

    // Translate image dimension info.
    slice.width      = (png_uint_32)sizes[0];
    slice.height     = (png_uint_32)sizes[1];
    slice.image_size = slice.width * slice.height;

    // Initialize PNG data structure.
    if (!png_save_init(&slice, dt)) {
        std::string msg;
        switch (slice.errcode) {
        case OPTNET_PNG_E_DATATYPE:
            msg = "png_save: Unsupported datatype";
            break;
        default:
            msg = "png_save: Unknown error";
        }
        throw_exception(io_error(msg));
    }

    // Fill in the additional fields.
    if (ext.voxel_size[0] != 0)
        slice.x_ppm  = (png_uint_32)(1.0 / ext.voxel_size[0]);
    if (ext.voxel_size[1] != 0)
        slice.y_ppm  = (png_uint_32)(1.0 / ext.voxel_size[1]);

    // Check if a palette is present and needed.
    slice.num_palette = (int)ext.palette.size();

    if (slice.num_palette > 0 && (
        (dt == OPTNET_DT_UNSIGNED_CHAR) || 
        (dt == OPTNET_DT_UNSIGNED_SHORT))) {

        slice.palette = new png_color[slice.num_palette];

        // Get a reference to the palette of a.
        const std::vector<rgb>& palette = ext.palette;
        // Copy the palette.
        for (int k = 0; k < slice.num_palette; ++k) {
            slice.palette[k].red   = palette[k][0];
            slice.palette[k].green = palette[k][1];
            slice.palette[k].blue  = palette[k][2];
        } // for

        slice.has_palette = true;
    }

    // Set textual information.
    if (!ext.title.empty()) {
        slice.title              = const_cast<png_charp>(ext.title.c_str());
        slice.title_length       = (png_size_t)ext.title.length();
    }
    if (!ext.author.empty()) {
        slice.author             = const_cast<png_charp>(ext.author.c_str());
        slice.author_length      = (png_size_t)ext.author.length();
    }
    if (!ext.description.empty()) {
        slice.description        = const_cast<png_charp>(ext.description.c_str());
        slice.description_length = (png_size_t)ext.description.length();
    }
    if (!ext.copyright.empty()) {
        slice.copyright          = const_cast<png_charp>(ext.copyright.c_str());
        slice.copyright_length   = (png_size_t)ext.copyright.length();
    }

    // Save all slices (default).
    first_slice    = 0;
    slices_to_save = sizes[2] * sizes[3] * sizes[4];

    // Save all requested slices.
    for (size_t i = first_slice; i < slices_to_save; ++i) {
        
        std::string filename = make_serial_filename(
            name,
            (unsigned int)i + 1,
            (unsigned int)slices_to_save,
            ".png"
            );

        fp = fopen(filename.c_str(), "wb");
        if (NULL == fp) {
            // Cleanup.
            if (NULL != slice.palette) delete [] slice.palette;

            throw_exception(io_error(
                std::string("png_save: Failed opening file: ") +
                filename
            ));
        } // if

        // Save PNG image info.
        if (!png_save_info(fp, &slice)) {
            // Cleanup.
            if (NULL != slice.palette) delete [] slice.palette;
            fclose(fp);

            std::string msg;
            switch (slice.errcode) {
            case OPTNET_PNG_E_WRITE_STRUCT:
                msg = "png_save: Failed png_create_write_struct";
                break;
            case OPTNET_PNG_E_LONG_PALETTE:
                msg = "png_save: The palette is too long.";
                break;
            case OPTNET_PNG_E_SETJMP:
                msg = "png_save: Failed writing image";
                break;
            default:
                msg = "png_save: Unknown error";
            } // switch
            throw_exception(io_error(msg));
        } // if

        // Get pointer to image data.
        offset = i * slice.image_size;
        slice.data = const_cast<unsigned char*>
                                (reinterpret_cast<const unsigned char*>
                                    (a.data() + offset));

        // Save PNG image data.
        if (!png_save_data(&slice)) {
            // Cleanup.
            if (NULL != slice.palette) delete [] slice.palette;
            fclose(fp);

            throw_exception(io_error(
                "png_load: Failed saving image data"
                ));
        }

        // Cleaning up.
        png_save_free(&slice);
        fclose(fp);

    }

    if (NULL != slice.palette)
        delete [] slice.palette;

} // png_save


    } // namespace
} // namespace

#endif
