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

//! ...
class LibMain : public CBase_LibMain {
  public:
    LibMain(CkArgMsg* msg);
};

//! ...
class MeshData {
  public:
    CProxy_NodeSearch m_proxy;
    int m_firstchunk;
    int m_nchare;
    void pup( PUP::er& p ) {
      p | m_proxy;
      p | m_firstchunk;
      p | m_nchare;
    }
};

//! Register a mesh to be part of mesh-to-mesh transfer
void addMesh( CkArrayID p, int nchare, CkCallback cb );

//! ...
void setSourceTets( CkArrayID p,
                    int index,
                    std::vector< std::size_t >* inpoel,
                    std::array< std::vector< double >, 3 >* coords,
                    const tk::Fields& u );

//! ...
void setDestPoints( CkArrayID p,
                    int index,
                    std::array< std::vector< double >, 3 >* coords,
                    tk::Fields& u,
                    CkCallback cb );

//! ...
class Transfer : public CBase_Transfer {

  public:
    //! ...
    Transfer();

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

    //! ...
    void setSourceTets( CkArrayID p,
                        int index,
                        std::vector< std::size_t >* inpoel,
                        std::array< std::vector< double >, 3 >* coords,
                        const tk::Fields& u );

    //! ...
    void setDestPoints( CkArrayID p,
                        int index,
                        std::array< std::vector< double >, 3 >* coords,
                        tk::Fields& u,
                        CkCallback cb );

    //! ...
    void distributeCollisions( int nColl, Collision* colls );

  private:
    //! ...
    std::unordered_map< CmiUInt8, MeshData > proxyMap;
    //! ...
    int current_chunk;
    //! ...
    CmiUInt8 m_sourcemesh, m_destmesh;
};

} // transfer::
