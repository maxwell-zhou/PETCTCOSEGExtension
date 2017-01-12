/*
 ==========================================================================
 |   
 |   $Id: debug.hxx 186 2005-02-27 02:04:22Z kangli $
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

#ifndef ___DEBUG_HXX___
#   define ___DEBUG_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/secure_s.hxx>

#   include <stdio.h>
#   ifdef __OPTNET_OS_WINNT__
#       include <windows.h>
#   endif

#   ifdef _UNICODE
#       include <stdlib.h>
#       include <string.h>
#   endif
    
// namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class debug
///  @brief Provide convenient debugging functions.
///////////////////////////////////////////////////////////////////////////
struct debug
{

#   if defined(_DEBUG) || defined(DEBUG) && !defined(NDEBUG)

    ///////////////////////////////////////////////////////////////////////
    static inline void outputv(const char* format, va_list argptr)
    {
        static const size_t _BUFSIZE = 4096;

#       if defined(__OPTNET_OS_WINNT__)

#           if defined(_UNICODE) || defined(UNICODE)

            // Unicode only
            wchar_t buffer[_BUFSIZE];
            size_t count = mbstowcs(NULL, format, 0) + 1;
            wchar_t *wformat = new wchar_t[count];
            mbstowcs(wformat, format, count);
            secure_vswprintf(buffer, sizeof(buffer), wformat, argptr);
            delete [] wformat;
            // End Unicode only

#           else

            // MBCS only
            char buffer[_BUFSIZE];
            secure_vsprintf(buffer, sizeof(buffer), format, argptr);
            // End MBCS only

#           endif

            // Send output to the debugger, if available.
            OutputDebugString(buffer);

#       else // Other operating systems.

        vfprintf(stdout, format, argptr);

#       endif
    }

    ///////////////////////////////////////////////////////////////////////
    static inline void output(const char* format, ...)
    {
        va_list argptr;
        va_start(argptr, format);
        outputv(format, argptr);
        va_end(argptr);
    }

#   else    // !_DEBUG

    static inline void outputv(const char*, va_list) {}
    static inline void output(const char*, ...) {}

#   endif   // _DEBUG

    ///////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    ///////////////////////////////////////////////////////////////////////
    debug() : m_plog(0) { }

    ///////////////////////////////////////////////////////////////////////
    ///  Destructor.
    ///////////////////////////////////////////////////////////////////////
    ~debug()            { log_end(); }

    ///////////////////////////////////////////////////////////////////////
    ///  Begin message logging.
    ///
    ///  @param filename  The filename of the log file, or NULL to use
    ///                   stdout.
    ///////////////////////////////////////////////////////////////////////
    void log_begin(const char* filename = NULL)
    {
        log_end(); // close up any files in use.
        if (NULL != filename)
            secure_fopen(&m_plog, filename, "at");
        if (NULL == m_plog) m_plog = stdout;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Prints a formatted message to the log file.
    ///////////////////////////////////////////////////////////////////////
    void log_vprintf(const char* format, va_list argptr)
    {
        if (m_plog == NULL) log_begin(NULL);
        vfprintf(m_plog, format, argptr);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Prints a formatted message to the log file.
    ///////////////////////////////////////////////////////////////////////
    void log_printf(const char* format, ...)
    {
        va_list argptr;
        va_start(argptr, format);
        log_vprintf(format, argptr);
        va_end(argptr);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Finish logging.
    ///////////////////////////////////////////////////////////////////////
    void log_end()
    {
        if (NULL != m_plog) {
            if (m_plog != stdout)
                fclose(m_plog);
            m_plog = NULL;
        }
    }


private:
    FILE* m_plog;
};

} // namespace

#   ifdef __OPTNET_ENABLE_TRACE__
#       define OPTNET_TRACE printf
#   else
        inline void fake_printf(const char*, ...) {}
#       define OPTNET_TRACE fake_printf
#   endif

#endif
