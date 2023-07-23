// *****************************************************************************
/*!
  \file      src/Main/Init.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Common initialization routines for main() functions for multiple
     exectuables
  \details   Common initialization routines for main() functions for multiple
     exectuables. The functions in this file are used by multiple execitables
     to ensure code-reuse and a uniform screen-output.
*/
// *****************************************************************************
#pragma once

#include <string>
#include <vector>

#include "Print.hpp"
#include "ProcessException.hpp"

#include "NoWarning/charm++.hpp"

namespace tk {

//! Executable types for which an ascii logo is available in tk::Print
enum class HeaderType : uint8_t { INCITER
                                , UNITTEST
                                , MESHCONV
                                };

//! Wrapper for the standard C library's gettimeofday() from
std::string curtime();

//! Echo program header
void echoHeader( HeaderType header );

//! Echo build environment
void echoBuildEnv( const std::string& executable );

//! Echo runtime environment
void echoRunEnv( int argc, char** argv, int quiescence );

//! Finalize function for different executables
void finalize( const std::vector< tk::Timer >& timer,
               std::vector< std::pair< std::string,
                                       tk::Timer::Watch > >& timestamp,
               bool clean = true );

} // tk::
