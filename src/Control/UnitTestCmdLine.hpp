// *****************************************************************************
/*!
  \file      src/Control/UnitTestCmdLine.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     UnitTest's command line
*/
// *****************************************************************************
#pragma once

#include <getopt.h>

#include "NoWarning/charm++.hpp"

#include "XystConfig.hpp"
#include "Exception.hpp"
#include "Print.hpp"
#include "TaggedTuple.hpp"
#include "TaggedTupleDeepPrint.hpp"
#include "Writer.hpp"

namespace tag {
DEFTAG( group );
DEFTAG( quiescence );
} // tag::

namespace unittest {
//! UnitTest control facilitating user input to internal data transfer
namespace ctr {

//! Member data for tagged tuple
using CmdLineMembers = brigand::list<
    tag::group, std::string     // groups to test
  , tag::quiescence, bool       // enable quiescence detection
>;

//! CmdLine is a TaggedTuple specialized to MeshConv
class CmdLine : public tk::TaggedTuple< CmdLineMembers > {

  public:
    //! Contructor: initialize members with default constructors
    explicit CmdLine() = default;

    //! Contructor: parse unittest command line
    explicit CmdLine( int argc, char** argv ) {
      tk::Print print;

      // Process command line arguments
      int c; 
      while ((c = getopt( argc, argv, "h?g:qv" )) != -1) {
        switch (c) {
          case '?':
          case 'h':
          default:
            help( argv );
            CkExit();
            break;
          case 'g':
            get< tag::group >() = optarg;
            break;
          case 'q':
            get< tag::quiescence >() = true;
            break;
          case 'v':
            print << '\n';
            print.version( tk::unittest_executable(), tk::git_commit() );
            CkExit();
            break;
        }
      }

      if (optind != argc) {
        print << "\nA non-option was supplied";
        help( argv );
        CkExit( EXIT_FAILURE );
      }
    
      // Output state to file
      auto logfilename = tk::unittest_executable() + "_cmdline.log";
      tk::Writer log( logfilename );
      tk::print( log.stream(), *this );
    }

    //! Print help
    void help( char** argv ) {
      tk::Print() <<
        "\nUsage: " << argv[0] << " [OPTION]...\n"
        "\n"
        "  -h, -?       Print out this help\n"
        "  -g <group>   Select test-group to run\n"
        "  -q           Enable quiescence detection\n"
        "  -v           Print revision information\n"
        "\n";
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
};


} // ctr::
} // unittest::