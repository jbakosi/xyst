// *****************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/Parser.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     MeshConv's command line parser
  \details   This file defines the command-line argument parser for the mesh
     file converter, MeshConv.
*/
// *****************************************************************************

#include "NoWarning/pegtl.hpp"
#include "NoWarning/charm.hpp"

#include "XystConfig.hpp"
#include "Exception.hpp"
#include "Print.hpp"
#include "Keywords.hpp"
#include "MeshConv/Types.hpp"
#include "MeshConv/CmdLine/Parser.hpp"
#include "MeshConv/CmdLine/Grammar.hpp"

namespace tk {
namespace grm {

tk::Print g_print;

} // grm::
} // tk::

using meshconv::CmdLineParser;

CmdLineParser::CmdLineParser( int argc,
                              char** argv,
                              const tk::Print& print,
                              ctr::CmdLine& cmdline ) :
  StringParser( argc, argv )
// *****************************************************************************
//  Contructor: parse the command line for MeshConv
//! \param[in] argc Number of C-style character arrays in argv
//! \param[in] argv C-style character array of character arrays
//! \param[in] print Pretty printer
//! \param[inout] cmdline Command-line stack where data is stored from parsing
// *****************************************************************************
{
  // Create CmdLine (a tagged tuple) to store parsed input
  ctr::CmdLine cmd;

  // Reset parser's output stream to that of print's. This is so that mild
  // warnings emitted during parsing can be output using the pretty printer.
  // Usually, errors and warnings are simply accumulated during parsing and
  // printed during diagnostics after the parser has finished. However, in some
  // special cases we can provide a more user-friendly message right during
  // parsing since there is more information available to construct a more
  // sensible message. This is done in e.g., tk::grm::store_option. Resetting
  // the global g_print, to that of passed in as the constructor argument allows
  // not to have to create a new pretty printer, but use the existing one.
  tk::grm::g_print.reset( print.save() );

  // Parse command line string by populating the underlying tagged tuple
  tao::pegtl::memory_input<> in( m_string, "command line" );
  tao::pegtl::parse< cmd::read_string, tk::grm::action >( in, cmd );

  // Echo errors and warnings accumulated during parsing
  diagnostics( print, cmd.get< tag::error >() );

  // Strip command line (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  cmdline = std::move( cmd );

  // If we got here, the parser has succeeded
  print.item("Parsed command line", "success");

  // Print out help on all command-line arguments if the executable was invoked
  // without arguments or the help was requested
  const auto helpcmd = cmdline.get< tag::help >();
  if (argc == 1 || helpcmd) {
    print.help< tk::QUIET >( tk::meshconv_executable(),
                             cmdline.get< tag::cmdinfo >(),
                             "Command-line Parameters:", "-" );
    print.mandatory< tk::QUIET >(
     "The '--" + kw::input().string() + " <filename>' and the "
     "'--" + kw::output().string() + " <filename>' arguments are mandatory." );
    print.usage< tk::QUIET >(
      tk::meshconv_executable(),
      tk::meshconv_executable() + " -" + *kw::input().alias() + " in.msh -" +
        *kw::output().alias() + " out.exo",
      "will read data from 'in.msh' (in Gmsh format) and output it to "
      "out.exo' (in ExodusII format)" );
  }

  // Print out verbose help for a single keyword if requested
  const auto helpkw = cmdline.get< tag::helpkw >();
  if (!helpkw.keyword.empty())
    print.helpkw< tk::QUIET >( tk::meshconv_executable(), helpkw );

  // Print out version information if it was requested
  const auto version = cmdline.get< tag::version >();
  if (version)
    print.version< tk::QUIET >( tk::meshconv_executable(),
                                tk::git_commit() );

  // Immediately exit if any help was output or was called without any argument
  // or version info was requested with zero exit code
  if (argc == 1 || helpcmd || !helpkw.keyword.empty() || version) CkExit();

  // Make sure mandatory arguments are set
  auto ialias = kw::input().alias();
  auto oalias = kw::output().alias();
  ErrChk( !(cmdline.get< tag::io, tag::input >().empty()),
          "Mandatory input file not specified. "
          "Use '--" + kw::input().string() + " <filename>'" +
          ( ialias ? " or '-" + *ialias + " <filename>'" : "" ) + '.' );
  ErrChk( !(cmdline.get< tag::io, tag::output >().empty()),
          "Mandatory output file not specified. "
          "Use '--" + kw::output().string() + " <filename>'" +
          ( oalias ? " or '-" + *oalias + " <filename>'" : "" ) + '.' );
}
