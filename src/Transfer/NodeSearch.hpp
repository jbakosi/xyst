// *****************************************************************************
/*!
  \file      src/Transfer/NodeSearch.hpp
  \copyright 2020-2021 Charmworks, Inc.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Declarations pertaining to node search between 3d meshes
*/
// *****************************************************************************
#pragma once

#include "Fields.hpp"
#include "NoWarning/nodesearch.decl.h"

namespace transfer {

//! ...
class PotentialCollision {
  public:
    std::size_t src, dst;
    CkVector3d point;
    void pup( PUP::er& p ) {
      p | src;
      p | dst;
      p | point;
    }
};

//! Solution data interpolated from src to dst mesh
class SolutionData {
  public:
    std::size_t dst;
    std::vector< double > sol;
    void pup( PUP::er& p ) {
      p | dst;
      p | sol;
    }
};

//! NodeSearch chare array holding part of a mesh
class NodeSearch : public CBase_NodeSearch {

  public:
    //! Constructor
    explicit NodeSearch( CkArrayID p, MeshData mesh, CkCallback cb );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    explicit NodeSearch( CkMigrateMessage* m ) : CBase_NodeSearch( m ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Set the source mesh data
    void setSourceTets( const std::vector< std::size_t >& inpoel,
                        const std::array< std::vector< double >, 3 >& coord,
                        const tk::Fields& u,
                        const std::vector< double >& flag,
                        bool dir,
                        CkCallback cb );

    //! Set the destination mesh data
    void setDestPoints( const std::array< std::vector< double >, 3 >& coord,
                        tk::Fields& u,
                        std::vector< double >& flag,
                        bool trflag,
                        bool dir,
                        CkCallback cb );

    //! Process potential collisions in the destination mesh
    void processCollisions( const MeshData& src,
                            int nColls,
                            Collision* colls );

    //! Identify actual collisions in the source mesh
    void determineActualCollisions( CProxy_NodeSearch proxy,
                                    int index,
                                    int nColls,
                                    PotentialCollision* colls );

    //! Transfer the interpolated solution data to destination mesh
    void transferSolution( const std::vector< SolutionData >& sol );

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_firstchunk;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i NodeSearch object reference
    friend void operator|( PUP::er& p, NodeSearch& i ) { i.pup(p); }
    //@}

  private:
    //! The ID of my first chunk (used for collision detection library)
    int m_firstchunk;
    //! Pointer to element connectivity
    std::vector< std::size_t >* m_inpoel;
    //! Pointer to point coordinates
    std::array< std::vector< double >, 3 >* m_coord;
    //! Pointer to solution in mesh nodes
    tk::Fields* m_u;
    //! Pointer to transfer flags
    std::vector< double >* m_flag;
    //! Transfer flags if true
    bool m_trflag;
    //! Transfer direction: 0: background to overset, 1: overset to background
    bool m_dir;
    //! The number of messages sent by the dest mesh
    int m_numsent;
    //! The number of messages received by the dest mesh
    int m_numreceived;
    //! Set to nonzero once source mesh is notified that source side is done
    int m_srcnotified = 0;
    //! Callback to call when transfer is done
    CkCallback m_done;

    //! Initialize dest mesh solution with background data
    void background();

    //! Contribute vertex information to the collsion detection library
    void collideVertices();

    //! Contribute tet information to the collision detection library
    void collideTets() const;
};

} // transfer::
