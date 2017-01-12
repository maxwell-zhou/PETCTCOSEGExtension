/*
 ==========================================================================
 |   
 |   $Id: process_info.hxx 151 2005-02-08 22:53:02Z kangli $
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

#ifndef ___PROCESS_INFO_HXX___
#   define ___PROCESS_INFO_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>
#   include <stddef.h>
#   if defined(__OPTNET_OS_WINNT__)
#       include <windows.h>
#       include <psapi.h>
#       pragma comment(lib, "psapi.lib")
#   elif defined(__OPTNET_OS_LINUX__)
#       include <stdio.h>
#   elif defined(__OPTNET_OS_BSD__)
#       include <sys/time.h>
#       include <sys/resource.h>
#       include <unistd.h>
#   endif

/// @namespace optnet
namespace optnet {
    /// @namespace system
    namespace system {

class process_info
{
public:
    ///////////////////////////////////////////////////////////////////////
    /// Returns the memory usage of the current process.
    ///
    /// @return The memory used in bytes.
    ///////////////////////////////////////////////////////////////////////
    inline size_t get_memory_usage()
    {
#   ifdef __OPTNET_OS_WINNT__

        //
        // Windows
        //
        DWORD dwProcID;
        HANDLE hProcess;
        PROCESS_MEMORY_COUNTERS pmc;

        dwProcID = ::GetCurrentProcessId();

        // Print information about the memory usage of the process.
        hProcess = ::OpenProcess(PROCESS_QUERY_INFORMATION |
            PROCESS_VM_READ, FALSE, dwProcID);

        if (NULL == hProcess)
            return 0;
        
        if (::GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)) == FALSE)
            return 0;

        return (size_t)pmc.WorkingSetSize;

#   elif defined(__OPTNET_OS_LINUX__)

        //
        // Linux
        //
        FILE* fp;
        int size = 0, resident = 0, shared = 0, total = 0;
    
        fp = fopen("/proc/self/statm", "r");
        fscanf(fp, "%d %d %d",
            &size, &resident, &shared);
        fclose(fp);
        
        total = (size + shared)
                    * getpagesize(); // page size: 4kb
        
        return (size_t)(total);
        
#   else

        return 0; // Not implemented.

#   endif
    }
    
    
    ///////////////////////////////////////////////////////////////////////
    /// Returns the peak memory usage of the current process.
    ///
    /// @return The peak memory used in bytes.
    ///////////////////////////////////////////////////////////////////////
    inline size_t get_peak_memory_usage()
    {
#   ifdef __OPTNET_OS_WINNT__

        //
        // Windows
        //
        DWORD dwProcID;
        HANDLE hProcess;
        PROCESS_MEMORY_COUNTERS pmc;

        dwProcID = ::GetCurrentProcessId();

        // Print information about the memory usage of the process.
        hProcess = ::OpenProcess(PROCESS_QUERY_INFORMATION |
            PROCESS_VM_READ, FALSE, dwProcID);

        if (NULL == hProcess)
            return 0;
        
        if (::GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)) == FALSE)
            return 0;

        return (size_t)pmc.PeakWorkingSetSize;

#   elif defined(__OPTNET_OS_BSD__)

        //
        // BSD-like Unix
        //
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage); // does not work under Linux.
        
        return (size_t)(usage.ru_maxrss + usage.ru_ixrss);

#   else

        return 0; // Not implemented.

#   endif
    }
    
};

    } // namespace
} // namespace

#endif // ___PROCESS_INFO_HXX___
