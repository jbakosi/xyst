// *****************************************************************************
/*!
  \file      src/Transfer/Transfer.cpp
  \copyright 2020-2021 Charmworks, Inc.,
             2022-2025 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Definitions pertaining to transfer between 3D meshes
*/
// *****************************************************************************

#include <cassert>

#include "Transfer.hpp"
#include "NodeSearch.hpp"

namespace transfer {

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

CProxy_Transfer g_transferProxy;
CollideHandle g_collideHandle;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

static void
collisionHandler( [[maybe_unused]] void *param, int nColl, Collision *colls )
// *****************************************************************************
//! ...
// *****************************************************************************
{
  g_transferProxy.ckLocalBranch()->distributeCollisions( nColl, colls );
}

void
addMesh( CkArrayID p, int nchare, CkCallback cb )
// *****************************************************************************
//  Register a mesh to be part of mesh-to-mesh transfer
//! \param[in] p Charm++ host proxy that holds the mesh
//! \param[in] nchare Number of mesh partitions
//! \param[in] cb Callback to continue with once finished
// *****************************************************************************
{
  g_transferProxy[ 0 ].addMesh( p, nchare, cb );
}

void
setSourceTets( CkArrayID p,
               int index,
               std::vector< std::size_t >* inpoel,
               std::array< std::vector< double >, 3 >* coords,
               const tk::Fields& u )
// *****************************************************************************
//! ...
// *****************************************************************************
{
  g_transferProxy.ckLocalBranch()->setSourceTets(p, index, inpoel, coords, u);
}

void setDestPoints( CkArrayID p,
                    int index,
                    std::array< std::vector< double >, 3 >* coords,
                    tk::Fields& u,
                    CkCallback cb )
// *****************************************************************************
//! ...
// *****************************************************************************
{
  g_transferProxy.ckLocalBranch()->setDestPoints( p, index, coords, u, cb );
}

//! ...
LibMain::LibMain( CkArgMsg* msg ) {
  delete msg;
  g_transferProxy = CProxy_Transfer::ckNew();

  // TODO: Need to make sure this is actually correct
  CollideGrid3d gridMap(CkVector3d(0, 0, 0),CkVector3d(2, 100, 2));
  g_collideHandle = CollideCreate( gridMap,
                      CollideSerialClient( collisionHandler, nullptr ) );
}

//! ...
Transfer::Transfer() : current_chunk(0) {}

void
Transfer::addMesh( CkArrayID p, int nchare, CkCallback cb )
// *****************************************************************************
//  Register a mesh to be part of mesh-to-mesh transfer
//! \param[in] p Charm++ host proxy that holds the mesh
//! \param[in] nchare Number of mesh partitions
//! \param[in] cb Callback to continue with once finished
// *****************************************************************************
{
  auto id = static_cast<std::size_t>(CkGroupID(p).idx);
  if (proxyMap.count(id) == 0) {
    CkArrayOptions opts;
    opts.bindTo(p);
    opts.setNumInitial( nchare );
    MeshData mesh;
    mesh.m_nchare = nchare;
    mesh.m_firstchunk = current_chunk;
    mesh.m_proxy = CProxy_NodeSearch::ckNew(p, mesh, cb, opts);
    proxyMap[id] = mesh;
    current_chunk += nchare;
  } else {
    CkAbort("Uhoh...\n");
  }
}

void
Transfer::setMesh( CkArrayID p, const MeshData& d )
// *****************************************************************************
//! ...
// *****************************************************************************
{
  proxyMap[static_cast<std::size_t>(CkGroupID(p).idx)] = d;
}

void
Transfer::setDestPoints( CkArrayID p,
                            int index,
                            std::array< std::vector< double >, 3 >* coords,
                            tk::Fields& u,
                            CkCallback cb)
// *****************************************************************************
//! ...
// *****************************************************************************
{
  m_destmesh = static_cast<std::size_t>(CkGroupID(p).idx);
  NodeSearch* w = proxyMap[m_destmesh].m_proxy[index].ckLocal();
  assert( w );
  w->setDestPoints(coords, u, cb);
}

void
Transfer::setSourceTets( CkArrayID p,
                            int index,
                            std::vector< std::size_t >* inpoel,
                            std::array< std::vector< double >, 3 >* coords,
                            const tk::Fields& u )
// *****************************************************************************
//! ...
// *****************************************************************************
{
  m_sourcemesh = static_cast<std::size_t>(CkGroupID(p).idx);
  NodeSearch* w = proxyMap[m_sourcemesh].m_proxy[index].ckLocal();
  assert(w);
  w->setSourceTets(inpoel, coords, u);
}

void
Transfer::distributeCollisions( int nColl, Collision* colls )
// *****************************************************************************
//  Called when all potential collisions have been found, and now need to be
//  destributed to the chares in the destination mesh to determine actual
//  collisions.
//! \param[in] nColl Number of potential collisions found
//! \param[in] colls The list of potential collisions
// *****************************************************************************
{
  //CkPrintf("Collisions found: %i\n", nColl);
  auto first = static_cast< std::size_t >( proxyMap[m_destmesh].m_firstchunk );
  auto nchare = static_cast< std::size_t >( proxyMap[m_destmesh].m_nchare );
  std::vector< std::vector< Collision > > separated( nchare );

  // Separate collisions based on the destination mesh chare they belong to
  for (std::size_t i=0; i<static_cast<std::size_t>(nColl); ++i) {
    auto achunk = static_cast< std::size_t >( colls[i].A.chunk );
    auto bchunk = static_cast< std::size_t >( colls[i].B.chunk );
    if (achunk >= first && achunk < first + nchare) {
      separated[ achunk - first ].push_back(colls[i]);
    } else {
      separated[ bchunk - first ].push_back(colls[i]);
    }
  }

  // Send out each list to the destination chares for further processing
  for (std::size_t i=0; i<nchare; ++i) {
    //CkPrintf("Dest mesh chunk %i has %lu\n", i, separated[i].size());
    proxyMap[ m_destmesh ].m_proxy[ static_cast<int>(i) ].processCollisions(
        proxyMap[ m_sourcemesh ].m_proxy,
        proxyMap[ m_sourcemesh ].m_nchare,
        proxyMap[ m_sourcemesh ].m_firstchunk,
        static_cast< int >( separated[ i ].size() ),
        separated[ i ].data() );
  }
}

} // transfer::

#include "NoWarning/transfer.def.h"
