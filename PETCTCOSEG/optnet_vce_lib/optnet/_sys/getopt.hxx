/*
 ==========================================================================
 |   
 |   $Id$
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

#ifndef ___GETOPT_HXX___
#   define ___GETOPT_HXX___

#   include <wchar.h>
#   include <optnet/_base/except.hxx>
#   include <optnet/_utils/xstring.hxx>

/// @namespace optnet
namespace optnet {
    /// @namespace system
    namespace system {

template <typename _CharT>
class getopt
{
public:

    typedef _CharT  char_type;

    static const char_type* EMPTY_STRING;

    ///////////////////////////////////////////////////////////////////////
    ///
    ///////////////////////////////////////////////////////////////////////
    getopt(int argc, char_type** argv) :
        m_is_option(false), m_index(1), m_count(argc),
        m_item(NULL), m_option(NULL), m_value(NULL),
        m_argv(argv)
    {
    }

    ///////////////////////////////////////////////////////////////////////
    ///
    ///////////////////////////////////////////////////////////////////////
    bool next()
    {
        if (m_index >= m_count) {
            m_option = NULL;
            return false;
        }

        m_item = m_argv[m_index++];

        if (m_item[0] == (char_type)('-')) {

            m_is_option = true;

            // Short option begins with '-', long option begins with '--'.
            m_option = m_item + 1;
            if (*m_option == (char_type)('-')) ++m_option;

            m_value = m_option;
            while (*m_value != (char_type)('\0') &&
                *m_value != (char_type)('=')) {
                ++m_value;
            }

            if (*m_value == (char_type)('=')) {
                *m_value++ = (char_type)('\0');
                if (*m_value == (char_type)('\"')) {
                    // Resolve quoted strings.
                    char_type* p = ++m_value;
                    while (*p != (char_type)('\0') &&
                        *p != (char_type)('\"')) {
                        ++p;
                    }
                    *p = (char_type)('\0');
                }
            }
            else {
                m_value = const_cast<char_type*>(EMPTY_STRING);
            }
        }
        else {
            m_is_option = false;
            m_option    = NULL;
            m_value     = NULL;
        }

        return true;
    }

    ///////////////////////////////////////////////////////////////////////
    ///
    ///////////////////////////////////////////////////////////////////////
    bool match(const char_type* short_option,
               const char_type* long_option,
               bool             has_value = false
               ) const
    {
        return samestr(short_option, m_option) ||
               samestr(long_option, m_option) && 
               (!has_value || 0 != m_value);
    }

    ///////////////////////////////////////////////////////////////////////
    ///
    ///////////////////////////////////////////////////////////////////////
    inline const char_type* option()  const { return m_option;  }
    inline const char_type* value()   const { return m_value;   }
    inline const char_type* item()    const { return m_item;    }

    ///////////////////////////////////////////////////////////////////////
    const double value_as_double() const
    {
        double v;
        try {
            v = str_to_double(m_value);
        } catch(std::exception&) {
            char errmsg[4096];
            secure_sprintf(errmsg, 4096,
                "Could not convert '%s' to double in argument '%s'.",
                m_value, m_item);
            throw_exception(std::invalid_argument(errmsg));
        }
        return v;
    }

    ///////////////////////////////////////////////////////////////////////
    const long value_as_long() const
    {
        long v;
        try {
            v = str_to_long(m_value);
        } catch(std::exception&) {
            char errmsg[4096];
            secure_sprintf(errmsg, 4096,
                "Could not convert '%s' to long integer in argument '%s'.",
                m_value, m_item);
            throw_exception(std::invalid_argument(errmsg));
        }
        return v;
    }

    ///////////////////////////////////////////////////////////////////////
    const unsigned long value_as_unsigned_long() const
    {
        unsigned long v;
        try {
            v = str_to_unsigned_long(m_value);
        } catch(std::exception&) {
            char errmsg[4096];
            secure_sprintf(errmsg, 4096,
                "Could not convert '%s' to unsigned long integer in argument '%s'.",
                m_value, m_item);
            throw_exception(std::invalid_argument(errmsg));
        }
        return v;
    }

    ///////////////////////////////////////////////////////////////////////
    ///
    ///////////////////////////////////////////////////////////////////////
    inline bool is_option() const         { return m_is_option; }


private:
    
    ///////////////////////////////////////////////////////////////////////
    bool samestr(const char_type* str1,
                 const char_type* str2
                 ) const
    {
        if (str1 == NULL || str2 == NULL)
            return false;

        while (*str1 == *str2 &&
               *str2 != (char_type)('\0')) {
           ++str1; ++str2;
        }

        return (*str1 == *str2);
    }

    ///////////////////////////////////////////////////////////////////////
    bool        m_is_option;
    int         m_index, m_count;
    char_type   *m_item, *m_option, *m_value, **m_argv;
};

template <>
const getopt<char>::char_type* getopt<char>::EMPTY_STRING = "";

template <>
const getopt<wchar_t>::char_type* getopt<wchar_t>::EMPTY_STRING = L"";

    } // namespace
} // namespace

#endif // ___GETOPT_HXX___
