// *****************************************************************************
/*!
  \file      src/Control/MeshConvConfig.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     MeshConv's command line
*/
// *****************************************************************************
#pragma once

#include <getopt.h>
#include <csignal>

#include "NoWarning/charm++.hpp"

#include "XystConfig.hpp"
#include "Exception.hpp"
#include "Print.hpp"
#include "TaggedTuple.hpp"
#include "PrintTaggedTupleDeep.hpp"
#include "Writer.hpp"

namespace tag {
DEFTAG( input );
DEFTAG( output );
DEFTAG( reorder );
DEFTAG( quiescence );
} // tag::

namespace meshconv {
//! Mesh converter control facilitating user input to internal data transfer
namespace ctr {

//! Member data for tagged tuple
using ConfigMembers = brigand::list<
    tag::input, std::string     // input mesh file
  , tag::output, std::string    // output mesh file
  , tag::reorder, bool          // reorder nodes
  , tag::quiescence, bool       // enable quiescence detection
>;

//! Config is a TaggedTuple specialized to MeshConv
class Config : public tk::TaggedTuple< ConfigMembers > {

  public:
    //! Parse meshconv command line
    void cmdline( int argc, char** argv ) {
      if (argc == 1) {
        help( argv );
        CkExit( EXIT_FAILURE );
      }
      tk::Print print;

      // Process command line arguments
      int c; 
      while ((c = getopt( argc, argv, "h?i:o:qrs:v" )) != -1) {
        switch (c) {
          case '?':
          case 'h':
          default:
            help( argv );
            CkExit();
            break;
          case 'i':
            get< tag::input >() = optarg;
            break;
          case 'o':
            get< tag::output >() = optarg;
            break;
          case 'q':
            get< tag::quiescence >() = true;
            break;
          case 'r':
            get< tag::reorder >() = true;
            break;
          case 's':
            std::raise( std::stoi( optarg ) );
            break;
          case 'v':
            print << '\n';
            print.version( tk::meshconv_executable(), tk::git_commit() );
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
              "Mandatory input file not specified. Use -i <filename>." );
      ErrChk( not get< tag::output >().empty(),
              "Mandatory output file not specified. Use -o <filename>." );
    
      // Output state to file
      auto logfilename = tk::meshconv_executable() + "_cmdline.log";
      tk::Writer log( logfilename );
      tk::print( log.stream(), *this );
    }

    void help( char** argv ) {
      tk::Print() <<
        "\nUsage: " << argv[0] << " -i <in.msh> -o <out.exo> [OPTION]...\n"
        "\n"
        "  -h, -?       Print out this help\n"
        "  -i <in.msh>  Specify input file\n"
        "  -o <out.exo> Specify output file\n"
        "  -r           Reorder mesh nodes\n"
        "  -q           Enable quiescence detection\n"
        "  -v           Print revision information\n"
        "\n";
    }

    /** @name Pack/Unpack: Serialize Config object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) { tk::TaggedTuple< ConfigMembers >::pup(p); }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c Config object reference
    friend void operator|( PUP::er& p, Config& c ) { c.pup(p); }
    //@}
};

} // ctr::
} // meshconv::
