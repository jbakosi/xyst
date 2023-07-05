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
#ifndef Init_h
#define Init_h

#include <string>
#include <unordered_map>

#include "NoWarning/charm++.hpp"

#include "XystConfig.hpp"
#include "Exception.hpp"
#include "Print.hpp"
#include "ProcessException.hpp"
#include "Tags.hpp"

namespace tk {

//! Executable types for which an ascii logo is available in tk::Print
enum class HeaderType : uint8_t { INCITER=0,
                                  UNITTEST,
                                  MESHCONV,
                                };

//! Wrapper for the standard C library's gettimeofday() from
std::string curtime();

//! Echo program header
void echoHeader( const Print& print, HeaderType header );

//! Echo build environment
void echoBuildEnv( const Print& print, const std::string& executable );

//! Echo runtime environment
void echoRunEnv( const Print& print, int argc, char** argv,
                 bool verbose, bool quiescence, bool trace,
                 const std::string& screen_log, const std::string& input_log );

//! \brief Generic Main() used for all executables for code-reuse and a uniform
//!    output
//! \details The template arguments configure this Main class that is
//!   practically used instead of the usual main(). This allows code-reuse and a
//!   unfirom screen-output. The template arguments are:
//!   - Driver, specializaing the driver type to be created, see tk::Driver
//!   - Printer, specializaing the pretty printer type to use, see tk::Print
//!   - CmdLine, specializing the command line object storing data parsed from
//!     the command line
//! \param[in] argc Number of command-line arguments to executable
//! \param[in] argv C-style string array to command-line arguments to executable
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \param[in] header Header type enum indicating which executable header to
//!   print
//! \param[in] executable Name of the executable
//! \param[in] def Default log file name
//! \param[in] nrestart Number of times restarted
//! \return Instantiated driver object which can then be used to execute()
//!   whatever it is intended to drive
template< class Driver, class CmdLine >
Driver Main( int argc, char* argv[],
             const CmdLine& cmdline,
             HeaderType header,
             const std::string& executable,
             const std::string& def,
             int nrestart )
{
  // Create pretty printer
  tk::Print
    print( cmdline.logname( def, nrestart ),
           cmdline.template get< tag::verbose >() ? std::cout : std::clog );

  // Echo program header
  echoHeader( print, header );

  // Echo environment
  print.part( "Environment" );
  // Build environment
  echoBuildEnv( print, executable );
  // Runtime environment
  echoRunEnv( print, argc, argv, cmdline.template get< tag::verbose >(),
              cmdline.template get< tag::quiescence >(),
              cmdline.template get< tag::trace >(),
              cmdline.logname( def, nrestart ),
              executable + "_input.log" );

  // Create and return driver
  return Driver( cmdline, nrestart );
}

//! Generic Main Charm++ module constructor for all executables
//! \tparam MainProxy Main Charm++ chare proxy for the executable
//! \tparam CmdLine Executable-specific tagged tuple storing the rusult of the
//!    command line parser
//! \param[in,out] mp MainProxy to set for the main chare
//! \param[in] thisProxy 'thisProxy' to set as MainProxy
//! \param[in,out] timer Vector of timers, held by the main chare, in which to
//!   start the first timer, measuring the migration of global-scope data
//! \param[in] cmdline Command line grammar stack for the executable (assumed
//!   already parsed)
//! \param[in] quiescenceTarget Pre-created Charm++ callback to use as the
//!   target function to call if quiescence is detected
template< class MainProxy, class CmdLine >
void MainCtor( MainProxy& mp,
               const MainProxy& thisProxy,
               std::vector< tk::Timer >& timer,
               const CmdLine& cmdline,
               const CkCallback& quiescenceTarget )
{
  // Set Charm++ main proxy
  mp = thisProxy;

  // Optionally enable quiscence detection
  if (cmdline.template get< tag::quiescence >()) CkStartQD( quiescenceTarget );

  // Start new timer measuring the migration of global-scope data
  timer.emplace_back();
}

//! Generic finalize function for different executables
//! \param[in] cmdline Command line grammar stack for the executable
//! \param[in] timer Vector of timers, held by the main chare
//! \param[in,out] timestamp Vector of time stamps in h:m:s with labels
//! \param[in] def Default log file name
//! \param[in] nrestart Number of times restarted
//! \param[in] clean True if we should exit with a zero exit code, false to
//!   exit with a nonzero exit code
template< class CmdLine >
void finalize( const CmdLine& cmdline,
               const std::vector< tk::Timer >& timer,
               std::vector< std::pair< std::string,
                                       tk::Timer::Watch > >& timestamp,
               const std::string& def,
               int nrestart,
               bool clean = true )
{
  try {

    if (!timer.empty()) {
      timestamp.emplace_back( "Total runtime", timer[0].hms() );
       tk::Print print( cmdline.logname( def, nrestart ),
         cmdline.template get< tag::verbose >() ? std::cout : std::clog,
         std::ios_base::app );
      print.time( "Timers (h:m:s)", timestamp );
      print.endpart();
    }

    if (clean) CkExit(); else CkAbort("Failed");

  } catch (...) { tk::processExceptionCharm(); }
}

} // tk::

#endif // Init_h
