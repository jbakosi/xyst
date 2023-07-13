// *****************************************************************************
/*!
  \file      src/Main/Init.cpp
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

#include <ctime>
#include <unistd.h>

#include "XystConfig.hpp"
#include "Exception.hpp"
#include "Tags.hpp"
#include "Keywords.hpp"
#include "Init.hpp"

namespace tk {

static std::string workdir()
// *****************************************************************************
// Wrapper for POSIX API's getcwd() from unistd.h
//! \return A stirng containing the current working directory
// *****************************************************************************
{
  char cwd[1024];

  if ( getcwd(cwd, sizeof(cwd)) != nullptr )
    return std::string( cwd );
  else
    Throw( "Error from POSIX API's getcwd()" );
}

std::string curtime()
// *****************************************************************************
//  Wrapper for the standard C library's gettimeofday() from
//! \return A stirng containing the current date and time
// *****************************************************************************
{
  time_t current_time;
  char* c_time_string;

  // Obtain current time as seconds elapsed since the Epoch
  current_time = time( nullptr );

  if (current_time == static_cast<time_t>(-1))
    Throw( "Failure to compute the current time" );

  // Convert to local time format
  c_time_string = ctime(&current_time);

  if (c_time_string == nullptr)
    Throw( "Failure to convert the current time" );

  // Convert to std::string and remove trailing newline
  std::string str( c_time_string );
  str.erase( std::remove(str.begin(), str.end(), '\n'), str.end() );

  return str;
}

void echoHeader( const Print& print, HeaderType header )
// *****************************************************************************
//  Echo program header
//! \param[in] print Pretty printer
//! \param[in] header Header type enum indicating which header to print
// *****************************************************************************
{
  if ( header == HeaderType::INCITER )
    print.headerInciter();
  else if ( header == HeaderType::UNITTEST )
    print.headerUnitTest();
  else if ( header == HeaderType::MESHCONV )
    print.headerMeshConv();
  else
    Throw( "Header not available" );
}

void echoBuildEnv( const Print& print, const std::string& executable )
// *****************************************************************************
//  Echo build environment
//! \details Echo information read from build_dir/Base/Config.h filled by
//!    CMake based on src/Main/Config.h.in.
//! \param[in] print Pretty printer
//! \param[in] executable Name of the executable
// *****************************************************************************
{
  print.section( "Build environment" );
  print.item( "Hostname", build_hostname() );
  print.item( "Executable", executable );
  if (!git_commit().empty()) print.item( "Git sha1", git_commit() );
  print.item( "Build type", build_type() );

#ifdef NDEBUG
  print.item( "Asserts", "off (turn on: CMAKE_BUILD_TYPE=DEBUG)" );
#else
  print.item( "Asserts", "on (turn off: CMAKE_BUILD_TYPE=RELEASE)" );
#endif

  print.item( "C++ compiler", compiler() );
  print.item( "Build date", build_date() );
}

void echoRunEnv( const Print& print, int argc, char** argv,
                 bool quiescence, bool trace,
                 const std::string& input_log )
// *****************************************************************************
//  Echo runtime environment
//! \param[in] print Pretty printer
//! \param[in] argc Number of command-line arguments to executable
//! \param[in] argv C-style string array to command-line arguments to executable
//! \param[in] quiescence True if quiescence detection is enabled
//! \param[in] trace True if call and stack trace output is enabled
//! \param[in] input_log Input log file name
// *****************************************************************************
{
  print.section( "Run-time environment" );

  print.item( "Date, time", curtime() );
  print.item( "Work directory", workdir() );
  print.item( "Executable (relative to work dir)", argv[0] );

  print << "Command line arguments : \'";
  if (argc>1) {
    for (auto i=1; i<argc-1; ++i) print << argv[i] << ' ';
    print << argv[argc-1];
  }
  print << "'\n";

  print.item( "Input log file", input_log );
  print.item( "Number of processing elements",
              std::to_string( CkNumPes() ) + " (" +
              std::to_string( CkNumNodes() ) + 'x' +
              std::to_string( CkNumPes()/CkNumNodes() ) + ')' );
  print.item( "Quiescence detection, -" + *kw::quiescence::alias(),
              quiescence ? "on" : "off" );
  print.item( "Call and stack trace, -" + *kw::trace::alias(),
              trace ? "on" : "off" );
}

// *****************************************************************************
//  Finalize function for different executables
//! \param[in] timer Vector of timers, held by the main chare
//! \param[in,out] timestamp Vector of time stamps in h:m:s with labels
//! \param[in] clean True if we should exit with a zero exit code, false to
//!   exit with a nonzero exit code
// *****************************************************************************
void finalize( const std::vector< tk::Timer >& timer,
               std::vector< std::pair< std::string,
                                       tk::Timer::Watch > >& timestamp,
               bool clean )
{
  try {

    if (!timer.empty()) {
      timestamp.emplace_back( "Total runtime", timer[0].hms() );
      tk::Print().time( "Timers (h:m:s)", timestamp );
    }

    if (clean) CkExit(); else CkAbort("Failed");

  } catch (...) { tk::processExceptionCharm(); }
}

} // tk::
