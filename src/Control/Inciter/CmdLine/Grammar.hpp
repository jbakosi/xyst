// *****************************************************************************
/*!
  \file      src/Control/Inciter/CmdLine/Grammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef InciterCmdLineGrammar_h
#define InciterCmdLineGrammar_h

#include "CommonGrammar.hpp"
#include "Keywords.hpp"

namespace inciter {
//! Inciter command line grammar definition
namespace cmd {

  using namespace tao;

  //! Specialization of tk::grm::use for Inciter's command line parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::CmdLine::keywords::set >;

  // Inciter's CmdLine grammar

  //! Match and set non-blocking (migration) switch
  struct nonblocking :
         tk::grm::process_cmd_switch< use, kw::nonblocking,
                                      tag::nonblocking > {};


  //! Match and set benchmark switch (i.e., benchmark mode)
  struct benchmark :
         tk::grm::process_cmd_switch< use, kw::benchmark,
                                      tag::benchmark > {};

  //! Match and set feedback switch (i.e., feedback mode)
  struct feedback :
         tk::grm::process_cmd_switch< use, kw::feedback,
                                      tag::feedback > {};

  //! Match and set virtualization parameter
  struct virtualization :
         tk::grm::process_cmd< use, kw::virtualization,
                               tk::grm::Store< tag::virtualization >,
                               tk::grm::number,
                               tag::virtualization > {};

  //! Match and set io parameter
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< use, keyword,
                               tk::grm::Store< tag::io, io_tag >,
                               pegtl::any,
                               tag::io, io_tag > {};

  //! Match help on command-line parameters
  struct help :
         tk::grm::process_cmd_switch< use, kw::help,
                                      tag::help > {};


  //! Match on quiescence switch
  struct quiescence :
         tk::grm::process_cmd_switch< use, kw::quiescence,
                                      tag::quiescence > {};

  //! Match and set load-balancing frequency
  struct lbfreq :
         tk::grm::process_cmd< use, kw::lbfreq,
                               tk::grm::Store< tag::lbfreq >,
                               tk::grm::number,
                               tag::lbfreq > {};

  //! Match and set checkpoint/restartfrequency
  struct rsfreq :
         tk::grm::process_cmd< use, kw::rsfreq,
                               tk::grm::Store< tag::rsfreq >,
                               tk::grm::number,
                               tag::rsfreq > {};

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
         pegtl::sor< nonblocking,
                     benchmark,
                     feedback,
                     virtualization,
                     help,
                     quiescence,
                     lbfreq,
                     rsfreq,
                     trace,
                     version,
                     io< kw::control, tag::control >,
                     io< kw::input, tag::input >,
                     io< kw::output, tag::output >,
                     io< kw::diagnostics_cmd, tag::diag >,
                     io< kw::restart, tag::restart > > {};

  //! Grammar entry point: parse keywords until end of string
  struct read_string :
         tk::grm::read_string< keywords > {};

} // cmd::
} // inciter::

#endif // InciterCmdLineGrammar_h
