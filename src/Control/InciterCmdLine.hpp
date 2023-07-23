// *****************************************************************************
/*!
  \file      src/Control/InciterCmdLine.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's command line
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
DEFTAG( input );
DEFTAG( control );
DEFTAG( output );
DEFTAG( diag );
DEFTAG( checkpoint );
DEFTAG( quiescence );
DEFTAG( virt );
DEFTAG( nonblocking );
DEFTAG( benchmark );
DEFTAG( feedback );
DEFTAG( lbfreq );
DEFTAG( rsfreq );
} // tag::

namespace inciter {
//! Inciter control facilitating user input to internal data transfer
namespace ctr {

//! Member data for tagged tuple
using CmdLineMembers = brigand::list<
    tag::input, std::string       // input mesh
  , tag::control, std::string     // control file
  , tag::output, std::string      // field-output base filename
  , tag::diag, std::string        // diagnostics output filename
  , tag::checkpoint, std::string  // checkpoint directory
  , tag::quiescence, bool         // enable quiescence detection
  , tag::virt, double             // virtualization [0...1]
  , tag::nonblocking, bool        // non-blocking migration
  , tag::benchmark, bool          // no large file I/O
  , tag::feedback, bool           // extra feedback
  , tag::lbfreq, uint64_t         // load-balancing frequency
  , tag::rsfreq, uint64_t         // checkpoint/restart frequency
>;

//! CmdLine is a TaggedTuple specialized to Inciter
class CmdLine : public tk::TaggedTuple< CmdLineMembers > {

  public:
    //! Contructor: initialize members with default constructors
    explicit CmdLine() = default;

    //! Contructor: parse inciter command line
    explicit CmdLine( int argc, char** argv ) {
      // Defaults
      get< tag::output >() = "out";
      get< tag::diag >() = "diag";
      get< tag::checkpoint >() = "restart";
      get< tag::lbfreq >() = 1;
      get< tag::rsfreq >() = 1000;

      if (argc == 1) {
        help( argv );
        CkExit( EXIT_FAILURE );
      }
      tk::Print print;

      // Process command line arguments
      int c; 
      while ((c = getopt( argc, argv, "bc:d:fh?i:l:no:qr:s:u:v" )) != -1) {
        switch (c) {
          case '?':
          case 'h':
          default:
            help( argv );
            CkExit();
            break;
          case 'b':
            get< tag::benchmark >() = true;
            break;
          case 'c':
            get< tag::control >() = optarg;
            break;
          case 'd':
            get< tag::diag >() = optarg;
            break;
          case 'f':
            get< tag::feedback >() = true;
            break;
          case 'i':
            get< tag::input >() = optarg;
            break;
          case 'l':
            get< tag::lbfreq >() = std::stoul( optarg );
            break;
          case 'n':
            get< tag::nonblocking >() = true;
            break;
          case 'o':
            get< tag::output >() = optarg;
            break;
          case 'r':
            get< tag::rsfreq >() = std::stoul( optarg );
            break;
          case 'q':
            get< tag::quiescence >() = true;
            break;
          case 'u':
            get< tag::virt >() = std::stod( optarg );
            break;
          case 'v':
            print << '\n';
            print.version( tk::inciter_executable(), tk::git_commit() );
            CkExit();
            break;
        }
      }

      if (optind != argc) {
        print << "\nA non-option was supplied";
        help( argv );
        CkExit( EXIT_FAILURE );
      }
      ErrChk( not get< tag::input >().empty(),
              "Mandatory input mesh file not specified. Use -i <filename>." );
      ErrChk( not get< tag::control >().empty(),
              "Mandatory control file not specified. Use -c <filename>." );
    
      // Output state to file
      auto logfilename = tk::inciter_executable() + "_cmdline.log";
      tk::Writer log( logfilename );
      tk::print( log.stream(), *this );
    }

    void help( char** argv ) {
      tk::Print() <<
        "\nUsage: " << argv[0] << " -i <in.msh> -c <config.q> [OPTION]...\n"
        "\n"
        "  -h, -?        Print out this help\n"
        "  -b            Benchmark mode, "
                         "default: " << get< tag::benchmark >() << "\n" <<
        "  -c <config.q> Specify control file\n"
        "  -d <diag>     Specify diagnostics file, "
                         "default: " << get< tag::diag >() << "\n" <<
        "  -f            Extra feedback, "
                         "default: " << get< tag::feedback >() << "\n" <<
        "  -i <in.msh>   Specify input mesh file\n"
        "  -l <int>      Load balancing frequency, "
                         "default: " << get< tag::lbfreq >() << "\n" <<
        "  -n            Non-blocking migration, "
                         "default: " << get< tag::nonblocking >() << "\n" <<
        "  -o <outfile>  Base-filename for field output, "
                         "default: " << get< tag::output >() << "\n" <<
        "  -r <int>      Checkpoint frequency, "
                         "default: " << get< tag::rsfreq >() << "\n" <<
        "  -q            Enable quiescence detection, "
                         "default: " << get< tag::quiescence >() << "\n" <<
        "  -u <real>     Virtualization, "
                         "default: " << get< tag::virt >() << "\n" <<
        "  -v            Print revision information\n"
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
} // inciter::
