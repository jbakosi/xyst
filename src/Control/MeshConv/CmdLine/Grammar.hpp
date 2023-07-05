// *****************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/Grammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     MeshConv's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef MeshConvCmdLineGrammar_h
#define MeshConvCmdLineGrammar_h

#include "CommonGrammar.hpp"
#include "Keywords.hpp"

namespace meshconv {
//! Mesh converter command line grammar definition
namespace cmd {

  using namespace tao;

  //! Specialization of tk::grm::use for MeshConv's command line parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords::set >;

  // MeshConv's CmdLine state

  // MeshConv's CmdLine grammar

  //! Match and set verbose switch (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< use, kw::verbose, tag::verbose > {};

  //! Match and set reorder switch (i.e., reorder mesh nodes or not)
  struct reorder :
         tk::grm::process_cmd_switch< use, kw::reorder_cmd, tag::reorder > {};

  //! Match and set io parameter
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< use, keyword,
                               tk::grm::Store< tag::io, io_tag >,
                               pegtl::any,
                               tag::io, io_tag > {};

  //! Match help on command-line parameters
  struct help :
         tk::grm::process_cmd_switch< use, kw::help, tag::help > {};

  //! Match help on a single command-line or control file keyword
  struct helpkw :
         tk::grm::process_cmd< use, kw::helpkw,
                               tk::grm::helpkw,
                               pegtl::alnum > {};

  //! Match raise_signal
  struct raise_signal :
         tk::grm::process_cmd< use, kw::raise_signal,
                               tk::grm::Store< tag::signal >,
                               tk::grm::number > {};

  //! Match help on control file keywords
  struct quiescence :
         tk::grm::process_cmd_switch< use, kw::quiescence,
                                      tag::quiescence > {};

  //! Match switch on trace output
  struct trace :
         tk::grm::process_cmd_switch< use, kw::trace,
                                      tag::trace > {};

  //! Match switch on version output
  struct version :
         tk::grm::process_cmd_switch< use, kw::version,
                                      tag::version > {};

  //! Match all command line keywords
  struct keywords :
         pegtl::sor< verbose,
                     reorder,
                     help,
                     helpkw,
                     quiescence,
                     trace,
                     version,
                     raise_signal,
                     io< kw::input, tag::input >,
                     io< kw::screen, tag::screen >,
                     io< kw::output, tag::output > > {};

  //! Grammar entry point: parse keywords until end of string
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // meshconv::

#endif // MeshConvCmdLineGrammar_h
