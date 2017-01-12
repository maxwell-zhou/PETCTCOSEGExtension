/*
 ==========================================================================
 |   
 |   $Id: timer.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___TIMER_HXX___
#   define ___TIMER_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/config.h>

#   if defined(__OPTNET_OS_WINNT__)
#       include <windows.h>
#   else
#       include <ctime>
#   endif


/// @namespace optnet
namespace optnet {
    
    /// @namespace optnet::utils
    /// @brief The namespace that contains utility classes and funtions.
    namespace utils {


#   if defined(__OPTNET_OS_WINNT__)

///////////////////////////////////////////////////////////////////////////
///  @class timer
///  @brief Simple timer class for measuring elapsed time.
///
///  Under Windows, this class uses the Win32 high-resolution timing APIs:
///  QueryPerformanceFrequency() and QueryPerformanceCounter(). Under the
///  other platforms, this class use the C Standard Library clock()
///  function.
///////////////////////////////////////////////////////////////////////////
class timer
{
public:

    ///////////////////////////////////////////////////////////////////////
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    timer()
    {
        LARGE_INTEGER cfreq;

        // Set the clock frequency if it is not yet set.
        QueryPerformanceFrequency(&cfreq);
        m_cfreq = (double)cfreq.QuadPart;
        
        // Get the current value of the high-res performance counter.
        QueryPerformanceCounter(&m_start);
    }
    
    ///////////////////////////////////////////////////////////////////////
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    void restart()
    {
        // Get the current value of the high-res performance counter.
        QueryPerformanceCounter(&m_start);
    }
  
    ///////////////////////////////////////////////////////////////////////
    ///  Returns elapsed time since the timer is created or restarted.
    ///
    ///  @return (double) Elapsed time in seconds.
    ///
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    double elapsed() const
    {
        LARGE_INTEGER   now;
        QueryPerformanceCounter(&now);
        return double(now.QuadPart - m_start.QuadPart) / m_cfreq;
    }


private:
    
    // Frequency setting is based on the hardware clock that does not
    // change between calling, so set this one only once.
    DOUBLE          m_cfreq;

    // The starting time stored as a LARGE_INTEGER.
    LARGE_INTEGER   m_start;

};

#   else // !__OPTNET_OS_WINNT__

///////////////////////////////////////////////////////////////////////////
///  @class timer
///  @brief Simple timer class for measuring elapsed time.
///
///  Under Windows, this class uses the Win32 high-resolution timing APIs:
///  QueryPerformanceFrequency() and QueryPerformanceCounter(). Under the
///  other platforms, this class use the C Standard Library clock()
///  function.
///////////////////////////////////////////////////////////////////////////
class timer
{
public:

    ///////////////////////////////////////////////////////////////////////
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    timer()             { m_start = clock(); }
    
    ///////////////////////////////////////////////////////////////////////
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    void restart()      { m_start = clock(); }
  
    ///////////////////////////////////////////////////////////////////////
    ///  Returns elapsed time since the timer is created or restarted.
    ///
    ///  @return (double) Elapsed time in seconds.
    ///
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    double elapsed() const
    {
        return double(clock() - m_start) / CLOCKS_PER_SEC;
    }

private:

    // The starting time stored as a clock_t.
    clock_t m_start;
};

#   endif // __OPTNET_OS_WINNT__

    } // namespace
} // namespace

#endif
