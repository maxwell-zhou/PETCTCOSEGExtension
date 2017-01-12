/*
 ==========================================================================
 |   
 |   $Id: filename.hxx 80 2005-01-22 07:49:12Z Administrator $
 |
 |   Written by Kang Li <kangl@cmu.edu>
 |   College of Engineering Imaging Group
 |   The University of Iowa
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

#ifndef ___FILENAME_HXX__
#   define ___FILENAME_HXX__

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#       pragma warning(disable: 4786)
#   endif

#   include <string>
#   include <optnet/config.h>
#   include <optnet/_base/secure_s.hxx>
#   include <optnet/_utils/xstring.hxx>

/// @namespace optnet
namespace optnet { 
    /// @namespace optnet::utils
    namespace utils {

///////////////////////////////////////////////////////////////////////////
///  Extract file name without directory.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_name(const _CharType* path)
{
    return extract_file_name(std::basic_string<_CharType>(path));
}

///////////////////////////////////////////////////////////////////////////
///  Extract file name without directory.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_name(const std::basic_string<_CharType>& path)
{
    typename std::basic_string<_CharType>::size_type pos1, pos2;

#   ifdef __OPTNET_OS_WINNT__
    pos1 = path.rfind((_CharType)('\\'));
    pos2 = path.rfind((_CharType)('/' ));
    if (
      (std::basic_string<_CharType>::npos != pos2)
    && ((pos2 > pos1) || (std::basic_string<_CharType>::npos == pos1)))
        pos1 = pos2;
    pos2 = path.rfind((_CharType)(':' ));
    if (
      (std::basic_string<_CharType>::npos != pos2)
    && ((pos2 > pos1) || (std::basic_string<_CharType>::npos == pos1)))
        pos1 = pos2;
#   else // Unixes
    pos1 = path.rfind((_CharType)('/' ));
    pos2 = 0; // unused
#   endif
    if (std::basic_string<_CharType>::npos != pos1)
        return path.substr(pos1 + 1);
    // else
    return path;
}

///////////////////////////////////////////////////////////////////////////
///  Extract file basename, including directory
///    (e.g., ./myfile.txt becomes ./myfile).
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_base(const _CharType* path)
{
    return extract_file_base(std::basic_string<_CharType>(path));
}

///////////////////////////////////////////////////////////////////////////
///  Extract file basename, including directory
///    (e.g., ./myfile.txt becomes ./myfile).
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_base(const std::basic_string<_CharType>& path)
{
    typename std::basic_string<_CharType>::size_type pos, len;

#   ifdef __OPTNET_OS_WINNT__
    len = path.length();
    if (len == 0 || path[len-1] == '\\' ||
        path[len-1] == '/' || path[len-1] == ':')
        return std::basic_string<_CharType>();
#   else // Unixes
    len = path.length();
    if (len == 0 || path[len-1] == '/')
        return std::basic_string<_CharType>();
#   endif

    pos = path.rfind((_CharType)('.'));
    if (std::basic_string<_CharType>::npos != pos) {
        return path.substr(0, pos);
    }
    // else
    return path;
}

///////////////////////////////////////////////////////////////////////////
///  Extract file directory.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_dir (const _CharType* path)
{
    return extract_file_dir (std::basic_string<_CharType>(path));
}

///////////////////////////////////////////////////////////////////////////
///  Extract file directory.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_dir (const std::basic_string<_CharType>& path)
{
    typename std::basic_string<_CharType>::size_type pos1, pos2;

#   ifdef __OPTNET_OS_WINNT__
    pos1 = path.rfind((_CharType)('\\'));
    pos2 = path.rfind((_CharType)('/' ));
    if (
      (std::basic_string<_CharType>::npos != pos2)
    && ((pos2 > pos1) || (std::basic_string<_CharType>::npos == pos1)))
        pos1 = pos2;
    pos2 = path.rfind((_CharType)(':' ));
    if (
      (std::basic_string<_CharType>::npos != pos2)
    && ((pos2 > pos1) || (std::basic_string<_CharType>::npos == pos1)))
        pos1 = pos2;
#   else // Unixes
    pos1 = path.rfind((_CharType)('/' ));
    pos2 = 0; // unused
#   endif
    if (std::basic_string<_CharType>::npos != pos1)
        return path.substr(0, pos1 + 1);
    // else
    return std::basic_string<_CharType>();
}

///////////////////////////////////////////////////////////////////////////
///  Extract file extension.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_ext (const _CharType* path)
{
    return extract_file_ext (std::basic_string<_CharType>(path));
}

///////////////////////////////////////////////////////////////////////////
///  Extract file extension.
///////////////////////////////////////////////////////////////////////////
template <typename _CharType>
inline std::basic_string<_CharType>
extract_file_ext (const std::basic_string<_CharType>& path)
{
    typename std::basic_string<_CharType>::size_type pos;

    pos = path.rfind((_CharType)('.'));
    if (std::basic_string<_CharType>::npos != pos)
        return path.substr(pos + 1);
    // else
    return std::basic_string<_CharType>();
}

///////////////////////////////////////////////////////////////////////////
/// Generate a serially numbered filename, e.g. pix001.png
///////////////////////////////////////////////////////////////////////////
inline std::string
make_serial_filename(const std::string& filename,
                     unsigned int       index,
                     unsigned int       total,
                     const char*        ext
                     )
{
    std::string name;

    if ((total > 0) && (index <= total)) {
        char tmp1[32], tmp2[32];

        if (str_ends_with_no_case(filename, ext))
            name = extract_file_base(filename);
        else
            name = filename;

        if (total > 1) {
            secure_sprintf(tmp1, sizeof(tmp1), "%d", total);
            secure_sprintf(tmp1, sizeof(tmp1), "%%0%dd", strlen(tmp1));
            secure_sprintf(tmp2, sizeof(tmp2), tmp1, index); 
            name += tmp2;
        }

        name += ext;
    }

    return name;
}


    } // namespace
} // namespace

#endif
