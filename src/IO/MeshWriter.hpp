// *****************************************************************************
/*!
  \file      src/IO/MeshWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             2022-2023 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ group for outputing mesh data to file
  \details   Charm++ group declaration used to output data associated to
     unstructured meshes to file(s). Charm++ chares (work units) send mesh and
     field data associated to mesh entities to the MeshWriter class defined here
     to write the data to file(s).
*/
// *****************************************************************************
#ifndef MeshWriter_h
#define MeshWriter_h

#include <vector>
#include <string>
#include <tuple>
#include <map>

#include "Types.hpp"
#include "Centering.hpp"
#include "UnsMesh.hpp"

#include "NoWarning/meshwriter.decl.h"

namespace tk {

//! Charm++ group used to output particle data to file in parallel
class MeshWriter : public CBase_MeshWriter {

  public:

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif

    //! Constructor: set some defaults that stay constant at all times
    //! \param[in] benchmark True of benchmark mode. No field output happens in
    //!   benchmark mode. This (and associated if tests) are here so client code
    //!   does not have to deal with this.
    //! \param[in] nmesh Total number of meshes
    MeshWriter( bool benchmark, std::size_t nmesh ) :
      m_benchmark( benchmark ),
      m_nchare( 0 ),
      m_nmesh( nmesh ) {}

    //! Migrate constructor
    explicit MeshWriter( CkMigrateMessage* m ) : CBase_MeshWriter( m ) {}

    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Set the total number of chares
    //! \param[in] n Total number of chares across the whole problem
    void nchare( int n ) { m_nchare = n; }

    //! Output unstructured mesh into file
    void write( std::size_t meshid,
                bool meshoutput,
                bool fieldoutput,
                uint64_t itr,
                uint64_t itf,
                tk::real time,
                int chareid,
                const std::string& basefilename,
                const std::vector< std::size_t >& inpoel,
                const UnsMesh::Coords& coord,
                const std::map< int, std::vector< std::size_t > >& bface,
                const std::map< int, std::vector< std::size_t > >& bnode,
                const std::vector< std::size_t >& triinpoel,
                const std::vector< std::string >& elemfieldnames,
                const std::vector< std::string >& nodefieldnames,
                const std::vector< std::string >& elemsurfnames,
                const std::vector< std::string >& nodesurfnames,
                const std::vector< std::vector< tk::real > >& elemfields,
                const std::vector< std::vector< tk::real > >& nodefields,
                const std::vector< std::vector< tk::real > >& elemsurfs,
                const std::vector< std::vector< tk::real > >& nodesurfs,
                const std::set< int >& outsets,
                CkCallback c );

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \note This is a Charm++ group, pup() is thus only for
    //!    checkpoint/restart.
    void pup( PUP::er &p ) override {
      p | m_benchmark;
      p | m_nchare;
      p | m_nmesh;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] m MeshWriter object reference
    friend void operator|( PUP::er& p, MeshWriter& m ) { m.pup(p); }
    //@}

  private:
    //! True if benchmark mode
    bool m_benchmark;
    //! Total number chares across the whole problem
    int m_nchare;
    //! Total number of meshes
    std::size_t m_nmesh;

    //! Compute filename
    std::string filename( const std::string& basefilename,
                          std::size_t meshid,
                          uint64_t itr,
                          int chareid,
                          int surfid = 0 ) const;
};

} // tk::

#endif // MeshWriter_h
