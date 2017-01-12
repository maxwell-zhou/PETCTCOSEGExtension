/*
 ==========================================================================
 |   
 |   $Id: except.hxx 156 2005-02-09 23:35:57Z kangli $
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

#ifndef ___EXCEPT_HXX___
#   define ___EXCEPT_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <optnet/config.h>

# include <stdio.h>

#   ifndef ___OPTNET_NO_EXCEPTIONS__
#       include <stdexcept>
#       include <string>
#   endif

///  @namespace optnet
namespace optnet {

#   ifdef ___OPTNET_NO_EXCEPTIONS__

    template<class _E>
        inline void throw_exception(const _E& e) { /*dummy*/}

#   else // !___OPTNET_NO_EXCEPTIONS__

    template<class _E>
        inline void throw_exception(const _E& e) { throw e; }

#   endif

#   ifndef ___OPTNET_NO_EXCEPTIONS__
    /// @namespace optnet::io
    namespace io {

        ///////////////////////////////////////////////////////////////////
        ///  @class io_error
        ///  @brief The class serves as the base class for all exceptions
        ///         thrown to report an I/O error.
        ///////////////////////////////////////////////////////////////////
        class io_error: public std::runtime_error
        {
        public:
            io_error(const std::string& message) :
                std::runtime_error(message)
            {
            }
        };
    } // namespace
#   endif

inline void
runtime_assert(const char* expr, const char* file, unsigned int line)
{
    std::string errmsg("Assertion failed: ");
    errmsg.append(expr);
    if (NULL != file) {
        const char* format = "\nFile: %s\nLine: %d";
        char  buffer[4096];
#   if defined(__OPTNET_SECURE_STR__) // secure version
        sprintf_s(buffer, sizeof(buffer), format, file, line);
#   else
        sprintf(buffer, format, file, line);
#   endif
        errmsg.append(buffer);
    }
#   ifndef ___OPTNET_NO_EXCEPTIONS__
    throw_exception(std::runtime_error(errmsg.c_str()));
#   else
    //FIXME
#   endif
}

} // namespace

/* ===================================================================== */
/*   EXCEPTION HANDLING MACROS                                           */
/* ===================================================================== */

#   ifdef RUNTIME_ASSERT
#       undef RUNTIME_ASSERT
#   endif

#   ifdef  NDEBUG
#       define RUNTIME_ASSERT(expr) \
        (void)((expr) || (runtime_assert(#expr, NULL, 0), 0))
#   else
#       define RUNTIME_ASSERT(expr) \
        (void)((expr) || (runtime_assert(#expr, __FILE__, __LINE__), 0))
#   endif

#   ifndef ___OPTNET_NO_EXCEPTIONS__
#       define __OPTNET_TRY try
#       define __OPTNET_CATCH_AND_RETURN            \
            catch (optnet::io::io_error&) {         \
                return OPTNET_E_IO_ERROR;           \
            }                                       \
            catch (std::invalid_argument&) {        \
                return OPTNET_E_INVALID_ARG;        \
            }                                       \
            catch (std::runtime_error&) {           \
                return OPTNET_E_RUNTIME_ERROR;      \
            }                                       \
            catch (std::exception&) {               \
                return OPTNET_E_FAILED;             \
            }
#   else
#       define __OPTNET_TRY
#       define __OPTNET_CATCH_AND_RETURN
#   endif

#endif 
