// *****************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/CmdLine.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
#            2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     UnitTest's command line
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the command-line parsing of the unit test
     suite, UnitTest.
*/
// *****************************************************************************
#ifndef UnitTestCmdLine_h
#define UnitTestCmdLine_h

#include <string>

#include <brigand/algorithms/for_each.hpp>

#include "Macro.hpp"
#include "XystConfig.hpp"
#include "TaggedTuple.hpp"
#include "HelpFactory.hpp"
#include "Keywords.hpp"
#include "UnitTest/Types.hpp"
#include "PrintUtil.hpp"

namespace unittest {
//! UnitTest control facilitating user input to internal data transfer
namespace ctr {

//! Member data for tagged tuple
using CmdLineMembers = brigand::list<
    tag::io,         ios
  , tag::verbose,    bool
  , tag::chare,      bool
  , tag::help,       bool
  , tag::quiescence, bool
  , tag::trace,      bool
  , tag::version,    bool
  , tag::license,    bool
  , tag::cmdinfo,    tk::ctr::HelpFactory
  , tag::ctrinfo,    tk::ctr::HelpFactory
  , tag::helpkw,     tk::ctr::HelpKw
  , tag::group,      std::string
  , tag::error,      std::vector< std::string >
>;

//! CmdLine is a TaggedTuple specialized to UnitTest
//! \details The stack is a tagged tuple
//! \see Base/TaggedTuple.h
//! \see Control/UnitTest/Types.h
class CmdLine : public tk::TaggedTuple< CmdLineMembers > {

  public:
    //! \brief UnitTest command-line keywords
    //! \see tk::grm::use and its documentation
    using keywords = tk::cmd_keywords< kw::verbose
                                     , kw::charestate
                                     , kw::help
                                     , kw::helpkw
                                     , kw::group
                                     , kw::screen
                                     , kw::quiescence
                                     , kw::trace
                                     , kw::version
                                     , kw::license
                                     >;

    //! Set of tags to ignore when printing this CmdLine
    using ignore =
      brigand::set< tag::cmdinfo
                  , tag::ctrinfo
                  , tag::helpkw >;

    //! \brief Constructor: set defaults.
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type. While there is a
    //!   ctrinfo parameter, it is unused here, since unittest does not have a
    //!   control file parser.
    //! \see walker::ctr::CmdLine
    CmdLine() {
      get< tag::io, tag::screen >() =
        tk::baselogname( tk::unittest_executable() );
      get< tag::verbose >() = false; // Use quiet output by default
      get< tag::chare >() = false; // No chare state output by default
      get< tag::trace >() = true; // Output call and stack trace by default
      get< tag::version >() = false; // Do not display version info by default
      get< tag::license >() = false; // Do not display license info by default
      // Initialize help: fill from own keywords
      brigand::for_each< keywords::set >( tk::ctr::Info(get<tag::cmdinfo>()) );
    }

    /** @name Pack/Unpack: Serialize CmdLine object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) { tk::TaggedTuple< CmdLineMembers >::pup(p); } 
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c CmdLine object reference
    friend void operator|( PUP::er& p, CmdLine& c ) { c.pup(p); }
    //@}

    //! Compute and return log file name
    //! \param[in] def Default log file name (so we don't mess with user's)
    //! \param[in] nrestart Number of times restarted
    //! \return Log file name
    std::string logname( const std::string& def, int nrestart ) const {
      if (get< tag::io, tag::screen >() != def)
        return get< tag::io, tag::screen >();
      else
        return tk::logname( tk::unittest_executable(), nrestart );
    }
};

} // ctr::
} // unittest::

#endif // UnitTestCmdLine_h
