#ifndef ___SECURE_S_HXX___
#   define ___SECURE_S_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/except.hxx>
#   include <stdarg.h>
#   include <stdlib.h>
#   include <errno.h>
#   include <stdio.h>

/// @namespace optnet
namespace optnet {

//
// Macros
//
#   if defined(__OPTNET_SECURE_STR__) // secure version
    /// When available, uses more secure version of sscanf.
#       define secure_sscanf sscanf_s
#   else // insecure version
    /// When available, uses more secure version of sscanf.
#       define secure_sscanf sscanf
#   endif

//
// Functions
//

///////////////////////////////////////////////////////////////////////////
///  More secure version of fopen.
///////////////////////////////////////////////////////////////////////////
inline int secure_fopen(FILE**      ppfile,
                        const char* filename,
                        const char* mode
                        )
{
#   if defined(__OPTNET_SECURE_STR__) // secure version
    return fopen_s(ppfile, filename, mode);
#   else // insecure version
    if (NULL == ppfile) return EINVAL;
    *ppfile = fopen(filename, mode);
    return (NULL == *ppfile) ? 
        EINVAL : 0;
#   endif
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of strtok.
///////////////////////////////////////////////////////////////////////////
inline char* secure_strtok(char*       token,
                           const char* delimit,
                           char**      context
                           )
{
#   if defined(__OPTNET_SECURE_STR__) // secure version
    return strtok_s(token, delimit, context);
#   else // insecure version
    OPTNET_UNUSED(context); // unused
    return strtok(token, delimit);
#   endif
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of vsprintf.
///////////////////////////////////////////////////////////////////////////
inline int secure_vsprintf(char*       buffer,
                           size_t      size_in_bytes,
                           const char* format,
                           va_list     argptr
                           )
{
#   if defined(__OPTNET_SECURE_STR__) // secure version
    return vsprintf_s(buffer, size_in_bytes, format, argptr);
#   elif defined(__OPTNET_HAS_VSNPRINTF__) // secure version
    return vsnprintf(buffer, size_in_bytes, format, argptr);
#   else // insecure version
    OPTNET_UNUSED(size_in_bytes);
    return vsprintf(buffer, format, argptr);
#   endif
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of vsprintf.
///////////////////////////////////////////////////////////////////////////
inline int secure_vswprintf(wchar_t*       buffer,
                            size_t         size_in_bytes,
                            const wchar_t* format,
                            va_list        argptr
                            )
{
#   if defined(__OPTNET_SECURE_STR__) // secure version
    return vswprintf_s(buffer, size_in_bytes, format, argptr);
#   elif defined(__OPTNET_HAS_VSNPRINTF__) // secure version
    return vsnwprintf(buffer, size_in_bytes, format, argptr);
/*
	#   else // insecure version
    OPTNET_UNUSED(size_in_bytes);
    return vswprintf(buffer, format, argptr);
	*/
#	else               //Debug
	return 0;
#   endif
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of sprintf.
///////////////////////////////////////////////////////////////////////////
inline int secure_sprintf(char*       buffer,
                          size_t      size_in_bytes,
                          const char* format,
                          ...
                          )
{
    va_list argptr;
    va_start(argptr, format);
    int ret = secure_vsprintf(buffer, size_in_bytes, format, argptr);
    va_end(argptr);
    return ret;
}

///////////////////////////////////////////////////////////////////////////
///  More secure version of wcstombs.
///////////////////////////////////////////////////////////////////////////
inline size_t secure_wcstombs(char*          mbstr,
                              const wchar_t* wcstr,
                              size_t         count
                              )
{
#   if defined(__OPTNET_SECURE_STR__) // secure version
    size_t converted = 0;
    errcode err = wcstombs_s(&converted, mbstr, count, wcstr, count);
    if (err != 0) {
        throw_exception(
            std::invalid_argument("Call to wcstombs_s() failed."));
    }
    return converted;
#   else // insecure version
    return wcstombs(mbstr, wcstr, count);
#   endif
}

} // namespace

#endif // ___SECURE_S_HXX___
