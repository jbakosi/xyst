// *****************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Parser.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     UnitTest's command line parser
  \details   This file defines the command-line argument parser for the unit
     test suite, UnitTest.
*/
// *****************************************************************************

#include "NoWarning/pegtl.hpp"

#include "Print.hpp"
#include "UnitTest/Types.hpp"
#include "UnitTest/CmdLine/Parser.hpp"
#include "UnitTest/CmdLine/Grammar.hpp"
#include "UnitTest/CmdLine/CmdLine.hpp"

using unittest::CmdLineParser;

CmdLineParser::CmdLineParser( int argc,
                              char** argv,
                              ctr::CmdLine& cmdline,
                              bool& helped ) :
  StringParser( argc, argv )
// *****************************************************************************
//  Contructor: parse the command line for UnitTest
//! \param[in] argc Number of C-style character arrays in argv
//! \param[in] argv C-style character array of character arrays
//! \param[inout] cmdline Command-line stack where data is stored from parsing
//! \param[inout] helped Boolean indicating if command-line help was requested
// *****************************************************************************
{
  // Create CmdLine (a tagged tuple) to store parsed input
  ctr::CmdLine cmd;

  // Parse command line string by populating the underlying tagged tuple
  tao::pegtl::memory_input<> in( m_string, "command line" );
  tao::pegtl::parse< cmd::read_string, tk::grm::action >( in, cmd );

  tk::Print print;

  // Echo errors and warnings accumulated during parsing
  diagnostics( print, cmd.get< tag::error >() );

  // Strip command line (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  cmdline = std::move( cmd );

  // Print out help on all command-line arguments if requested
  const auto helpcmd = cmdline.get< tag::help >();
  if (helpcmd) {
    print.help( tk::unittest_executable(),
                cmdline.get< tag::cmdinfo >(),
                "Command-line Parameters:", "-" );
   print.mandatory( "None of the arguments are mandatory." );
   print.usage( "charmrun +p4 " + tk::unittest_executable(),
                "will execute all unit tests on 4 CPUs" );
  }

  // Print out version information if it was requested
  const auto version = cmdline.get< tag::version >();
  if (version) {
    print.version( tk::unittest_executable(), tk::git_commit() );
  }

  // Will exit in main chare constructor if any help was output
  if (cmdline.get< tag::help >() || version) {
    helped = true;
  } else {
    helped = false;
  }
}
