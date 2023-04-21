// *****************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/CmdLine.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     MeshConv's command line definition
  \details   This file defines the heterogeneous stack that is used for storing
     the data from user input during the command-line parsing of the mesh
     converter, MeshConv.
*/
// *****************************************************************************
#ifndef MeshConvCmdLine_h
#define MeshConvCmdLine_h

#include <string>

#include <brigand/algorithms/for_each.hpp>

#include "XystConfig.hpp"
#include "Macro.hpp"
#include "Keywords.hpp"
#include "HelpFactory.hpp"
#include "MeshConv/Types.hpp"

namespace meshconv {
//! Mesh converter control facilitating user input to internal data transfer
namespace ctr {

//! Member data for tagged tuple
using CmdLineMembers = brigand::list<
    tag::io,         ios
  , tag::verbose,    bool
  , tag::chare,      bool
  , tag::reorder,    bool
  , tag::help,       bool
  , tag::quiescence, bool
  , tag::trace,      bool
  , tag::version,    bool
  , tag::signal,     int
  , tag::cmdinfo,    tk::ctr::HelpFactory
  , tag::ctrinfo,    tk::ctr::HelpFactory
  , tag::helpkw,     tk::ctr::HelpKw
  , tag::error,      std::vector< std::string >
>;

//! \brief CmdLine is a TaggedTuple specialized to MeshConv
//! \details The stack is a tagged tuple, a hierarchical heterogeneous data
//!    structure where all parsed information is stored.
//! \see Base/TaggedTuple.h
//! \see Control/MeshConv/Types.h
class CmdLine : public tk::TaggedTuple< CmdLineMembers > {

  public:
    //! \brief MeshConv command-line keywords
    //! \see tk::grm::use and its documentation
    using keywords = tk::cmd_keywords< kw::verbose
                                     , kw::charestate
                                     , kw::help
                                     , kw::helpkw
                                     , kw::input
                                     , kw::output
                                     , kw::screen
                                     , kw::reorder_cmd
                                     , kw::quiescence
                                     , kw::trace
                                     , kw::version
                                     , kw::raise_signal
                                     >;

    //! Set of tags to ignore when printing this CmdLine
    using ignore =
      brigand::set< tag::cmdinfo
                  , tag::ctrinfo
                  , tag::helpkw >;

    //! \brief Constructor: set defaults.
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type. While there is a
    //!   ctrinfo parameter, it is unused here, since meshconv does not have a
    //!   control file parser.
    //! \see walker::ctr::CmdLine
    CmdLine() {
      get< tag::io, tag::screen >() =
        tk::baselogname( tk::meshconv_executable() );
      get< tag::verbose >() = false; // Use quiet output by default
      get< tag::chare >() = false; // No chare state output by default
      get< tag::reorder >() = false; // Do not reorder by default
      get< tag::trace >() = true; // Output call and stack trace by default
      get< tag::version >() = false; // Do not display version info by default
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
        return tk::logname( tk::meshconv_executable(), nrestart );
    }
};

} // ctr::
} // meshconv::

#endif // MeshConvCmdLine_h
