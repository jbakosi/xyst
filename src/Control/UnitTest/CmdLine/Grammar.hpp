// *****************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Grammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     UnitTest's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef UnitTestCmdLineGrammar_h
#define UnitTestCmdLineGrammar_h

#include "CommonGrammar.hpp"
#include "Keywords.hpp"

namespace unittest {
//! UnitTest command line grammar definition
namespace cmd {

  using namespace tao;

  //! \brief Specialization of tk::grm::use for UnitTest's command line parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords::set >;

  //! \brief Match help on command-line parameters
  struct help :
         tk::grm::process_cmd_switch< use, kw::help,
                                      tag::help > {};

  //! \brief Match test group name(s) and only run those
  struct group :
         tk::grm::process_cmd< use, kw::group,
                               tk::grm::Store< tag::group >,
                               pegtl::any,
                               tag::group > {};

  //! Match switch on quiescence
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

  //! \brief Match and set io parameter
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< use, keyword,
                               tk::grm::Store< tag::io, io_tag >,
                               pegtl::any,
                               tag::io, io_tag > {};

  //! \brief Match all command line keywords
  struct keywords :
         pegtl::sor< help
                   , group
                   , quiescence
                   , trace
                   , version
                   > {};

  //! \brief Grammar entry point: parse keywords until end of string
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // unittest::

#endif // UnitTestCmdLineGrammar_h
