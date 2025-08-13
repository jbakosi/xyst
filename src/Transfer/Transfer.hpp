// *****************************************************************************
/*!
  \file      src/Transfer/Transfer.hpp
  \copyright 2020-2021 Charmworks, Inc.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Declarations pertaining to transfer between 3D meshes
*/
// *****************************************************************************

#include "NoWarning/transfer.decl.h"
#include "NoWarning/collidecharm.h"
#include "Fields.hpp"

namespace transfer {

//! Single Charm++ chare used to initialize mesh-to-mesh transfer
class LibTransfer : public CBase_LibTransfer {
  public:
    //! Constructor: initialize mesh-to-mesh transfer
    explicit LibTransfer( CkArgMsg* msg );
};

//! Mesh configuration for a mesh involved in solution transfer
//! \details Lightweight data structure identifying a mesh
class MeshData {
  public:
    //! Host proxy of mesh
    CProxy_NodeSearch proxy;
    //! Starting chare ID of mesh partition
    int firstchunk;
    //! Number of mesh partitions
    int nchare;
    //! Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      p | proxy;
      p | firstchunk;
      p | nchare;
    }
};

//! API for registering a mesh to be part of mesh-to-mesh transfer
void addMesh( CkArrayID p, int nchare, CkCallback cb );

//! API for configuring source mesh
void setSourceTets( CkArrayID p,
                    int chare,
                    const std::vector< std::size_t >& inpoel,
                    const std::array< std::vector< double >, 3 >& coord,
                    const tk::Fields& u,
                    const std::vector< double >& flag,
                    bool dir,
                    CkCallback cb );

//! API for configuring destination mesh
void setDestPoints( CkArrayID p,
                    int chare,
                    const std::array< std::vector< double >, 3 >& coord,
                    tk::Fields& u,
                    std::vector< double >& flag,
                    bool trflag,
                    bool dir,
                    CkCallback cb );

//! ...
class Transfer : public CBase_Transfer {

  public:
    //! Constructor
    explicit Transfer() = default;

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! ...
    explicit Transfer( CkMigrateMessage* m ) : CBase_Transfer( m ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Register a mesh to be part of mesh-to-mesh transfer
    void addMesh( CkArrayID p, int nchare, CkCallback cb );

    //! ...
    void setMesh( CkArrayID p, const MeshData& d );

    //! Configure source mesh
    void setSourceTets( CkArrayID p,
                        int chare,
                        const std::vector< std::size_t >& inpoel,
                        const std::array< std::vector< double >, 3 >& coord,
                        const tk::Fields& u,
                        const std::vector< double >& flag,
                        bool dir,
                        CkCallback cb );

    //! Configure destination mesh
    void setDestPoints( CkArrayID p,
                        int chare,
                        const std::array< std::vector< double >, 3 >& coord,
                        tk::Fields& u,
                        std::vector< double >& flag,
                        bool trflag,
                        bool dir,
                        CkCallback cb );

    //! ...
    void distributeCollisions( int nColl, Collision* colls );

  private:
    //! Mesh configuration for each mesh involved in solution transfer
    std::unordered_map< CmiUInt8, MeshData > m_proxyMap;
    //! ...
    int m_current_chunk = 0;
    //! Source mesh id
    CmiUInt8 m_src;
    //! Destination mesh id
    CmiUInt8 m_dst;
};

} // transfer::
