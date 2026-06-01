//
// Created by sunsh on 4/26/2022.
//

#include "mem_usage.h"

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <unistd.h>
#endif
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

//////////////////////////////////////////////////////////////////////////////
// Get from stackoverflow: https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0


void mem_usage::process_mem_usage(double& vm_usage_KB, double& resident_set_KB)
{
    vm_usage_KB     = 0.0;
    resident_set_KB = 0.0;

#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (!GetProcessMemoryInfo(GetCurrentProcess(),
                              reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc),
                              sizeof(pmc))) {
        return;
    }

    vm_usage_KB = static_cast<double>(pmc.PrivateUsage) / 1024.0;
    resident_set_KB = static_cast<double>(pmc.WorkingSetSize) / 1024.0;
    return;
#else
    using std::ios_base;
    using std::ifstream;
    using std::string;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat",ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage_KB     = vsize / 1024.0;
    resident_set_KB = rss * page_size_kb;
#endif
}