// *****************************************************************************
/*!
  \file      src/Control/Inciter/CmdLine/Parser.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's command line parser
  \details   This file defines the command-line argument parser for the
     computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************

#include "NoWarning/pegtl.hpp"
#include "NoWarning/charm.hpp"

#include "XystConfig.hpp"
#include "Exception.hpp"
#include "Print.hpp"
#include "Keywords.hpp"
#include "Inciter/Types.hpp"
#include "Inciter/CmdLine/Parser.hpp"
#include "Inciter/CmdLine/Grammar.hpp"
#include "Inciter/CmdLine/CmdLine.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::CmdLineParser;

CmdLineParser::CmdLineParser( int argc, char** argv, ctr::CmdLine& cmdline ) :
  StringParser( argc, argv )
// *****************************************************************************
//  Contructor: parse the command line for Inciter
//! \param[in] argc Number of C-style character arrays in argv
//! \param[in] argv C-style character array of character arrays
//! \param[in,out] cmdline Command-line stack where data is stored from parsing
// *****************************************************************************
{
  // Create CmdLine (a tagged tuple) to store parsed input
  ctr::CmdLine cmd( g_inputdeck.get< tag::cmd, tag::ctrinfo >() );

  // Parse command line string by populating the underlying tagged tuple:
  tao::pegtl::memory_input<> in( m_string, "command line" );
  tao::pegtl::parse< cmd::read_string, tk::grm::action >( in, cmd );

  tk::Print print;

  // Echo errors and warnings accumulated during parsing
  diagnostics( print, cmd.get< tag::error >() );

  // Strip command line (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  cmdline = std::move( cmd );

  // Print out help on all command-line arguments if the executable was invoked
  // without arguments or the help was requested
  const auto helpcmd = cmdline.get< tag::help >();
  if (argc == 1 || helpcmd) {
    print.help( tk::inciter_executable(),
                 cmdline.get< tag::cmdinfo >(),
                 "Command-line Parameters:", "-" );
    print.mandatory(
     "The '--" + kw::input().string() + " <filename>' and the "
     "'--" + kw::control().string() + " <filename>' arguments are mandatory." );
    print.usage(
      "charmrun +p4 " + tk::inciter_executable() + " -" + *kw::control().alias()
      + " vort.q -" + *kw::input().alias() + " unitcube.exo",
      "will execute the simulation configured in the control file 'vort.q' "
      "using the mesh in 'unitcube.exo' on 4 CPUs" );
  }

  // Print out help on all control file keywords if they were requested
  const auto helpctr = cmdline.get< tag::helpctr >();
  if (helpctr)
    print.help( tk::inciter_executable(),
                cmdline.get< tag::ctrinfo >(),
                "Control File Keywords:" );

  // Print out verbose help for a single keyword if requested
  const auto helpkw = cmdline.get< tag::helpkw >();
  if (!helpkw.keyword.empty())
    print.helpkw( tk::inciter_executable(), helpkw );

  // Print out version information if it was requested
  const auto version = cmdline.get< tag::version >();
  if (version) {
    print.version( tk::inciter_executable(), tk::git_commit() );
  }

  // Immediately exit if any help was output or was called without any argument
  // or version info was requested with zero exit code
  if (argc == 1 || helpcmd || helpctr || !helpkw.keyword.empty() || version) {
    CkExit();
  }

  // Make sure mandatory arguments are set
  auto ctralias = kw::control().alias();
  ErrChk( !(cmdline.get< tag::io, tag::control >().empty()),
          "Mandatory control file not specified. "
          "Use '--" + kw::control().string() + " <filename>'" +
          ( ctralias ? " or '-" + *ctralias + " <filename>'" : "" ) + '.' );
}
