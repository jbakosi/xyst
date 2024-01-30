// *****************************************************************************
/*!
  \file      src/Base/Benchmark.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Macro definitions for fine-grained benchmarking
*/
// *****************************************************************************
#pragma once

#include "sys/time.h"

//! Start-time macro for fine-grained profiling. Put this in the beginning of
//! the section of code to be profiled.
#define STARTTIME \
struct timeval START_TIME, END_TIME; \
long total_usecs; \
gettimeofday( &START_TIME, NULL );

//! End-time macro for fine-grained profiling. Put this at the end of the
//! section of code to be profiled.
#define ENDTIME \
gettimeofday( &END_TIME, NULL ); \
total_usecs = (END_TIME.tv_sec-START_TIME.tv_sec) * 1000000 + (END_TIME.tv_usec-START_TIME.tv_usec); \
printf("Total time was %ld uSec.\n", total_usecs);
