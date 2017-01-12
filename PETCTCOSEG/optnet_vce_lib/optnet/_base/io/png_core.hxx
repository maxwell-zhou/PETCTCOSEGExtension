/*
 ==========================================================================
 |   
 |   $Id: png_core.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___PNG_CORE_HXX__
#   define ___PNG_CORE_HXX__

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4611)
#       pragma warning(disable: 4786)
#   endif

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       if defined(_DEBUG)
#           if defined(_MT)
#               if defined(_DLL)
#                   pragma comment(lib, "libpngMDd.lib")
#                   pragma comment(lib, "zlibMDd.lib")
#               else
#                   pragma comment(lib, "libpngMTd.lib")
#                   pragma comment(lib, "zlibMTd.lib")
#               endif
#           else
#               pragma comment(lib, "libpngd.lib")
#               pragma comment(lib, "zlibd.lib")
#           endif
#       else // !_DEBUG
#           if defined(_MT)
#               if defined(_DLL)
#                   pragma comment(lib, "libpngMD.lib")
#                   pragma comment(lib, "zlibMD.lib")
#               else
#                   pragma comment(lib, "libpngMT.lib")
#                   pragma comment(lib, "zlibMT.lib")
#               endif
#           else
#               pragma comment(lib, "libpng.lib")
#               pragma comment(lib, "zlib.lib")
#           endif
#       endif
#   endif

#   include <optnet/config.h> // configurations

#   include <cstdio>
#   include <cstdlib>

#   ifdef __OPTNET_USE_SYSTEM_LIBPNG__
#       include <png.h>
#   else
#       include <optnet/_base/io/detail/png.h>
#   endif

#   include <optnet/_base/datatype.hxx>

//
// namespace optnet::io::detail
//
namespace optnet { namespace io { namespace detail {

static const size_t PNG_BYTES_TO_CHECK = 4;

enum optnet_png_errcode
{
    OPTNET_PNG_S_OK             =  0,
    OPTNET_PNG_E_FORMAT         = -1,
    OPTNET_PNG_E_READ_STRUCT    = -2,
    OPTNET_PNG_E_WRITE_STRUCT   = -3,
    OPTNET_PNG_E_INFO_STRUCT    = -4,
    OPTNET_PNG_E_SETJMP         = -5,
    OPTNET_PNG_E_SIZE           = -6,
    OPTNET_PNG_E_DATATYPE       = -7,
    OPTNET_PNG_E_LONG_PALETTE   = -8
};

struct optnet_png_slice
{
    png_structp     png_ptr;
    png_infop       info_ptr;
    png_uint_32     width;
    png_uint_32     height;
    png_uint_32     image_size;
    int             bit_depth;
    int             color_type;
    int             interlace_type;
    png_uint_32     x_ppm;
    png_uint_32     y_ppm;
    png_colorp      palette;
    int             num_palette;
    bool            has_palette;
    png_uint_32     image_bytes;
    png_uint_32     pixel_bytes;
    png_uint_32     row_bytes;
    png_bytep*      row_pointers;
    unsigned char*  data;
    int             errcode;
    // embedded text
    png_charp       title;
    png_size_t      title_length;
    png_charp       author;
    png_size_t      author_length;
    png_charp       description;
    png_size_t      description_length;
    png_charp       copyright;
    png_size_t      copyright_length;
};

///////////////////////////////////////////////////////////////////////////
// Check if the file is in PNG format.
inline bool png_check_format(FILE* fp, unsigned int* psig_read)
{
    char  buf[PNG_BYTES_TO_CHECK];

    // First check the file pointer.
    if (NULL == fp) return false;

    // Read in some of the signature bytes
    if (fread(buf, 1, PNG_BYTES_TO_CHECK, fp) != PNG_BYTES_TO_CHECK)
        return false;

    if (NULL != psig_read)
        *psig_read = PNG_BYTES_TO_CHECK;

    // Compare the first PNG_BYTES_TO_CHECK bytes of the signature.
    //  Return nonzero (true) if they match.
    return (0 == png_sig_cmp((png_bytep)buf,
                             (png_size_t)0,
                             PNG_BYTES_TO_CHECK)
                             );
}

///////////////////////////////////////////////////////////////////////////
inline bool
png_load_info(FILE*             fp,
              optnet_png_slice* slice,
              unsigned int      sig_read,
              long              datatype
              )
{
    int         i;
    png_structp png_ptr;
    png_infop   info_ptr;

    assert(NULL != slice);
    assert(NULL != fp);
    assert(sig_read <= 8);

    // Clear slice structure.
    memset(slice, 0, sizeof(optnet_png_slice));

    /* Create and initialize the png_struct with the desired error handler
     * functions.  If you want to use the default stderr and longjump method,
     * you can supply NULL for the last three parameters.  We also supply the
     * the compiler header file version, so that we know if the application
     * was compiled with a compatible version of the library.  REQUIRED
     */
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
        NULL, NULL, NULL);
    if (NULL == png_ptr) {
        slice->errcode = OPTNET_PNG_E_READ_STRUCT;
        return false;
    }

    /* Allocate/initialize the memory for image information.  REQUIRED. */
    info_ptr = png_create_info_struct(png_ptr);
    if (NULL == info_ptr) {
        png_destroy_read_struct(&png_ptr, png_infopp_NULL, png_infopp_NULL);
        slice->errcode = OPTNET_PNG_E_INFO_STRUCT;
        return false;
    }

    /* Set error handling if you are using the setjmp/longjmp method (this is
     * the normal method of doing things with libpng).  REQUIRED unless you
     * set up your own error handlers in the png_create_read_struct() earlier.
     */
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
        slice->errcode = OPTNET_PNG_E_SETJMP;
        return false;
    }

    /* Set up the input control if you are using standard C streams */
    png_init_io(png_ptr, fp);


    /* If we have already read some of the signature */
    png_set_sig_bytes(png_ptr, sig_read);

    /* The call to png_read_info() gives us all of the information from the
     * PNG file before the first IDAT (image data chunk).  REQUIRED
     */
    png_read_info(png_ptr, info_ptr);

    png_get_IHDR(png_ptr,
                 info_ptr,
                 &(slice->width),
                 &(slice->height),
                 &(slice->bit_depth),
                 &(slice->color_type),
                 &(slice->interlace_type),
                 int_p_NULL,
                 int_p_NULL
                 );
    slice->image_size = slice->width * slice->height;

    /* Extract multiple pixels with bit depths of 1, 2, and 4 from a single
     * byte into separate bytes (useful for paletted and grayscale images).
     */
    png_set_packing(png_ptr);

    /* Expand grayscale images to the full 8 bits from 1, 2, or 4 bits/pixel */
    if (slice->color_type == PNG_COLOR_TYPE_GRAY && 
        slice->bit_depth < 8) {
        png_set_gray_1_2_4_to_8(png_ptr);
        slice->bit_depth = 8; // change bit depth
    }

    png_color_16 my_background = {0, 255, 255, 255, 0};
    png_color_16p image_background;

    if (png_get_bKGD(png_ptr, info_ptr, &image_background))
        png_set_background(png_ptr, image_background,
          PNG_BACKGROUND_GAMMA_FILE, 1, 1.0);
    else
        png_set_background(png_ptr, &my_background,
          PNG_BACKGROUND_GAMMA_SCREEN, 0, 1.0);


    /* Determine if the given datatype matches the image. */
    bool supported = false;
    switch (datatype) {
    case OPTNET_DT_UNSIGNED_CHAR:
        if ((slice->bit_depth == 8) && (
            (slice->color_type == PNG_COLOR_TYPE_GRAY) ||
            (slice->color_type == PNG_COLOR_TYPE_PALETTE) ||
            (slice->color_type == PNG_COLOR_TYPE_GRAY_ALPHA) ||
            (slice->color_type == PNG_COLOR_TYPE_RGB) ||
            (slice->color_type == PNG_COLOR_TYPE_RGB_ALPHA))) {

            if (slice->color_type == PNG_COLOR_TYPE_PALETTE) {
                png_set_palette_to_rgb(png_ptr);
                png_set_rgb_to_gray_fixed(png_ptr, 1, -1, -1);
            }
            else if (slice->color_type == PNG_COLOR_TYPE_RGB ||
                slice->color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
                png_set_rgb_to_gray_fixed(png_ptr, 1, -1, -1);
            }

            /* Strip alpha bytes from the input data without combining with the
             * background (not recommended).
             */
            if (slice->color_type & PNG_COLOR_MASK_ALPHA) {
                png_set_strip_alpha(png_ptr);
            }

            slice->pixel_bytes = 1;
            supported = true;
        }
        break;
    case OPTNET_DT_UNSIGNED_SHORT:
        if ((slice->bit_depth == 16) && (
            (slice->color_type == PNG_COLOR_TYPE_GRAY) ||
            (slice->color_type == PNG_COLOR_TYPE_PALETTE) ||
            (slice->color_type == PNG_COLOR_TYPE_GRAY_ALPHA) ||
            (slice->color_type == PNG_COLOR_TYPE_RGB) ||
            (slice->color_type == PNG_COLOR_TYPE_RGB_ALPHA))) {

            if (slice->color_type == PNG_COLOR_TYPE_RGB ||
                slice->color_type == PNG_COLOR_TYPE_RGB_ALPHA)
                png_set_rgb_to_gray_fixed(png_ptr, 1, -1, -1);

            /* Strip alpha bytes from the input data without combining with the
             * background (not recommended).
             */
            if (slice->color_type & PNG_COLOR_MASK_ALPHA) {
                png_set_strip_alpha(png_ptr);
            }

            slice->pixel_bytes = 2;
            supported = true;
        }
        break;
    case OPTNET_DT_RGBA:
        if ((slice->bit_depth == 8) && (
            (slice->color_type == PNG_COLOR_TYPE_GRAY) ||
            (slice->color_type == PNG_COLOR_TYPE_GRAY_ALPHA) ||
            (slice->color_type == PNG_COLOR_TYPE_PALETTE) ||
            (slice->color_type == PNG_COLOR_TYPE_RGB) ||
            (slice->color_type == PNG_COLOR_TYPE_RGBA))) {

            if (slice->color_type == PNG_COLOR_TYPE_GRAY ||
                slice->color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
                png_set_gray_to_rgb(png_ptr);

            /* Expand paletted colors into true RGB triplets */
            if (slice->color_type == PNG_COLOR_TYPE_PALETTE)
                png_set_palette_to_rgb(png_ptr);

            /* Expand paletted or RGB images with transparency to full alpha channels
             * so the data will be available as RGBA quartets.
             */
            if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
                png_set_tRNS_to_alpha(png_ptr);

            png_set_filler(png_ptr, 0xFF, PNG_FILLER_AFTER);

            slice->pixel_bytes = 4;
            supported = true;
        }
        break;
    case OPTNET_DT_RGB:
        if ((slice->bit_depth == 8) && (
            (slice->color_type == PNG_COLOR_TYPE_GRAY) ||
            (slice->color_type == PNG_COLOR_TYPE_GRAY_ALPHA) ||
            (slice->color_type == PNG_COLOR_TYPE_PALETTE) ||
            (slice->color_type == PNG_COLOR_TYPE_RGB) ||
            (slice->color_type == PNG_COLOR_TYPE_RGBA))) {

            if (slice->color_type == PNG_COLOR_TYPE_GRAY ||
                slice->color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
                png_set_gray_to_rgb(png_ptr);

            /* Expand paletted colors into true RGB triplets */
            if (slice->color_type == PNG_COLOR_TYPE_PALETTE)
                png_set_palette_to_rgb(png_ptr);

            /* Strip alpha bytes from the input data without combining with the
             * background (not recommended).
             */
            if (slice->color_type & PNG_COLOR_MASK_ALPHA) {
                png_set_strip_alpha(png_ptr);
            }

            slice->pixel_bytes = 3;
            supported = true;
        }
        break;
    default:
        ;
    }

    if (!supported) {
        png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
        slice->errcode = OPTNET_PNG_E_FORMAT;
        return false;
    }

    /* Optional call to gamma correct and add the background to the palette
     * and update info structure.  REQUIRED if you are expecting libpng to
     * update the palette for you (ie you selected such a transform above).
     */
    png_read_update_info(png_ptr, info_ptr);

    /* Allocate the memory to hold the image using the fields of info_ptr. */
    slice->row_bytes = png_get_rowbytes(png_ptr, info_ptr);

    if (png_get_PLTE(png_ptr,
                     info_ptr,
                     &slice->palette,
                     &slice->num_palette
                     )) {
        slice->has_palette = true;
    }

    slice->image_bytes = slice->row_bytes * slice->height;

    if ((slice->pixel_bytes * slice->image_size)
            != slice->image_bytes) {
        // This should not happen.
        png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
        slice->errcode = OPTNET_PNG_E_SIZE;
        return false;
    }

    /* Get image metrics. */
    slice->x_ppm = png_get_x_pixels_per_meter(png_ptr, info_ptr);
    slice->y_ppm = png_get_y_pixels_per_meter(png_ptr, info_ptr);

    /* Get embedded textual information. */
    png_textp text_ptr;
    int       num_text;
    png_get_text(png_ptr, info_ptr, &text_ptr, &num_text);

    for (i = 0; i < num_text; ++i) {
        if (0 == strcmp(text_ptr[i].key, "Title")) {
            slice->title = text_ptr[i].text;
            slice->title_length = text_ptr[i].text_length;
        }
        else if (0 == strcmp(text_ptr[i].key, "Author")) {
            slice->author = text_ptr[i].text;
            slice->author_length = text_ptr[i].text_length;
        }
        else if (0 == strcmp(text_ptr[i].key, "Description")) {
            slice->description = text_ptr[i].text;
            slice->description_length = text_ptr[i].text_length;
        }
        else if (0 == strcmp(text_ptr[i].key, "Copyright")) {
            slice->copyright = text_ptr[i].text;
            slice->copyright_length = text_ptr[i].text_length;
        }
    }

    /* Save png data structure handles. */
    slice->info_ptr = info_ptr;
    slice->png_ptr = png_ptr;

    return true;
}


///////////////////////////////////////////////////////////////////////////
inline bool
png_load_data(optnet_png_slice* slice)
{
    assert(NULL != slice);
    assert(NULL != slice->png_ptr);
    assert(NULL != slice->info_ptr);

    // Reset error code (no error).
    slice->errcode = OPTNET_PNG_S_OK;

    slice->row_pointers = (png_bytep*)png_malloc(
        slice->png_ptr,
        slice->height * png_sizeof(png_bytep)
        );

    for (png_uint_32 row = 0; row < slice->height; ++row) {
        /* Initialize row pointers. */
        slice->row_pointers[row] = slice->data + row * slice->row_bytes;
    }

    /* Now it's time to read the image. */
    png_read_image(slice->png_ptr, slice->row_pointers);

    /* Read rest of file, and get additional chunks in info_ptr - REQUIRED */
    png_read_end(slice->png_ptr, slice->info_ptr);

    return true;
}

///////////////////////////////////////////////////////////////////////////
inline void
png_load_free(optnet_png_slice* slice)
{
    assert(NULL != slice);
    assert(NULL != slice->png_ptr);
    assert(NULL != slice->info_ptr);

    /* At this point you have read the entire image */
    /* clean up after the read, and free any memory allocated - REQUIRED */
    png_destroy_read_struct(&(slice->png_ptr),
                            &(slice->info_ptr),
                            png_infopp_NULL
                            );
}

///////////////////////////////////////////////////////////////////////////
inline bool
png_save_init(optnet_png_slice* slice,
              long              datatype
              )
{
    assert(NULL != slice);

    // Reset error code (no error).
    slice->errcode = OPTNET_PNG_S_OK;

    // Let's derive the color_type and bit_depth of the image.
    switch (datatype) {
    case OPTNET_DT_UNSIGNED_CHAR:   // unsigned char
        slice->bit_depth   = 8;
        slice->color_type  = PNG_COLOR_TYPE_GRAY;
        slice->pixel_bytes = 1;
        break;
    case OPTNET_DT_UNSIGNED_SHORT:  // unsigned short
        slice->bit_depth   = 16;
        slice->color_type  = PNG_COLOR_TYPE_GRAY;
        slice->pixel_bytes = 2;
        break;
    case OPTNET_DT_RGB:             // rgb
        slice->bit_depth   = 8;
        slice->color_type  = PNG_COLOR_TYPE_RGB;
        slice->pixel_bytes = 3;
        break;
    case OPTNET_DT_RGBA:            // rgba
        slice->bit_depth   = 8;
        slice->color_type  = PNG_COLOR_TYPE_RGBA;
        slice->pixel_bytes = 4;
        break;
    default:
        // unsupported datatype
        slice->errcode = OPTNET_PNG_E_DATATYPE;
        return false;
    }

    slice->row_bytes   = slice->pixel_bytes * slice->width;
    slice->image_bytes = slice->row_bytes * slice->height;

    return true;
}

///////////////////////////////////////////////////////////////////////////
inline bool
png_save_info(FILE*             fp,
              optnet_png_slice* slice
              )
{
    png_structp png_ptr;
    png_infop   info_ptr;

    assert(NULL != fp);
    assert(NULL != slice);

    // Reset error code (no error).
    slice->errcode = OPTNET_PNG_S_OK;

    /* Create and initialize the png_struct with the desired error handler
     * functions.  If you want to use the default stderr and longjump method,
     * you can supply NULL for the last three parameters.  We also check that
     * the library version is compatible with the one used at compile time,
     * in case we are using dynamically linked libraries.  REQUIRED.
     */
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
        (png_voidp)NULL, NULL, NULL);
    if (NULL == png_ptr) {
        slice->errcode = OPTNET_PNG_E_WRITE_STRUCT;
        return false;
    }

    /* Allocate/initialize the image information data.  REQUIRED */
    info_ptr = png_create_info_struct(png_ptr);
    if (NULL == info_ptr) {
        png_destroy_write_struct(&png_ptr,  png_infopp_NULL);
        return false;
    }

    /* Set error handling.  REQUIRED if you aren't supplying your own
     * error handling functions in the png_create_write_struct() call.
     */
    if (setjmp(png_jmpbuf(png_ptr))) {
        /* If we get here, we had a problem reading the file */
        png_destroy_write_struct(&png_ptr, &info_ptr);
        slice->errcode = OPTNET_PNG_E_SETJMP;
        return false;
    }

    /* set up the output control with standard C streams */
    png_init_io(png_ptr, fp);

    /* Set the image information here.  Width and height are up to 2^31,
     * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
     * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
     * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
     * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
     * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
     * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
     */
    png_set_IHDR(png_ptr,
                 info_ptr,
                 slice->width,
                 slice->height,
                 slice->bit_depth,
                 slice->color_type,
                 PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT
                 );

    if (slice->has_palette && (slice->num_palette > 0)) {
        if (slice->num_palette > PNG_MAX_PALETTE_LENGTH) {
            // Palette too long.
            png_destroy_write_struct(&png_ptr, &info_ptr);
            slice->errcode = OPTNET_PNG_E_LONG_PALETTE;
            return false;
        }
        /* set the palette if there is one.  REQUIRED for indexed-color images */
        png_set_PLTE(png_ptr, info_ptr, slice->palette, slice->num_palette);
    }

    /* Optionally write comments into the image */
    png_text text_ptr [5];
    int      num_text = 1;

    static const png_charp SOFTWARE = "OPTNET Library PNG Writer";
    text_ptr[0].key         = "Software";
    text_ptr[0].text        = SOFTWARE;
    text_ptr[0].text_length = strlen(SOFTWARE);
    text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;

    if ((NULL != slice->title) && (0 != slice->title_length)) {
        text_ptr[num_text].key         = "Title";
        text_ptr[num_text].text        = slice->title;
        text_ptr[num_text].text_length = slice->title_length;
        text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
        ++num_text;
    }
    if ((NULL != slice->author) && (0 != slice->author_length)) {
        text_ptr[num_text].key         = "Author";
        text_ptr[num_text].text        = slice->author;
        text_ptr[num_text].text_length = slice->author_length;
        text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
        ++num_text;
    }
    if ((NULL != slice->description) && (0 != slice->description_length)) {
        text_ptr[num_text].key         = "Description";
        text_ptr[num_text].text        = slice->description;
        text_ptr[num_text].text_length = slice->description_length;
        text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_zTXt;
        ++num_text;
    }
    if ((NULL != slice->copyright) && (0 != slice->copyright_length)) {
        text_ptr[num_text].key         = "Copyright";
        text_ptr[num_text].text        = slice->copyright;
        text_ptr[num_text].text_length = slice->copyright_length;
        text_ptr[num_text].compression = PNG_TEXT_COMPRESSION_NONE;
        ++num_text;
    }
    png_set_text(png_ptr, info_ptr, text_ptr, num_text);

    /* Set pixels/unit physical resolutions. */
    png_set_pHYs(png_ptr, info_ptr, slice->x_ppm, slice->y_ppm, 
        PNG_RESOLUTION_METER);

    /* Write the file header information.  REQUIRED */
    png_write_info(png_ptr, info_ptr);
    
    /* Save png data structure handles. */
    slice->info_ptr = info_ptr;
    slice->png_ptr = png_ptr;

    return true;
}

///////////////////////////////////////////////////////////////////////////
inline bool
png_save_data(optnet_png_slice* slice)
{
    assert(NULL != slice);

    // Now writing the actual image pixels.
    slice->row_pointers = (png_bytep*)png_malloc(
        slice->png_ptr,
        slice->height * png_sizeof(png_bytep)
        );

    for (png_uint_32 k = 0; k < slice->height; ++k) {
        // Initialize row pointers.
        slice->row_pointers[k] = slice->data + k * slice->row_bytes;
    }

    /* write out the entire image data in one call */
    png_write_image(slice->png_ptr, slice->row_pointers);

    /* It is REQUIRED to call this to finish writing the rest of the file */
    png_write_end(slice->png_ptr, slice->info_ptr);

    return true;
}

///////////////////////////////////////////////////////////////////////////
inline void
png_save_free(optnet_png_slice* slice)
{
    assert(NULL != slice);
    assert(NULL != slice->png_ptr);
    assert(NULL != slice->info_ptr);

    /* clean up after the write, and free any memory allocated */
    png_destroy_write_struct(&(slice->png_ptr), &(slice->info_ptr));
}

} } } // namespace

#endif
