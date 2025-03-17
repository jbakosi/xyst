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
#include <iostream>     // NOT NEEDED

#include "Transfer.hpp"
#include "NodeSearch.hpp"
#include "ContainerUtil.hpp"

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

LibTransfer::LibTransfer( CkArgMsg* msg )
// *****************************************************************************
//  Constructor: initialize mesh-to-mesh transfer
// *****************************************************************************
{
  delete msg;

  g_transferProxy = CProxy_Transfer::ckNew();

  // Need to make sure this is actually correct
  CollideGrid3d gridMap( CkVector3d(0, 0, 0), CkVector3d(2, 100, 2) );
  g_collideHandle = CollideCreate( gridMap,
                      CollideSerialClient( collisionHandler, nullptr ) );
}

void
addMesh( CkArrayID p, int nchare, CkCallback cb )
// *****************************************************************************
//  API for registering a mesh to be part of mesh-to-mesh transfer
//! \param[in] p Charm++ host proxy participating in mesh-to-mesh transfer
//! \param[in] nchare Number of mesh partitions
//! \param[in] cb Callback to continue with once finished
// *****************************************************************************
{
  g_transferProxy[0].addMesh( p, nchare, cb );
}

void
setSourceTets( CkArrayID p,
               int chare,
               const std::vector< std::size_t >& inpoel,
               const std::array< std::vector< double >, 3 >& coord,
               const tk::Fields& u,
               CkCallback cb )
// *****************************************************************************
//  API for configuring source mesh
//! \param[in] inpoel Source mesh connectivity
//! \param[in] coord Source mesh node coordinates
//! \param[in] u Source solution data
// *****************************************************************************
{
  g_transferProxy.ckLocalBranch()->setSourceTets( p, chare, inpoel, coord, u, cb );
}


void setDestPoints( CkArrayID p,
                    int chare,
                    const std::array< std::vector< double >, 3 >& coord,
                    tk::Fields& u,
                    CkCallback cb )
// *****************************************************************************
//  API for configuring destination mesh
//! \param[in] coord Destination mesh node coordinates
//! \param[in,out] u Destination mesh solution data
//! \param[in] cb Callback to call once this chare received all solution data
// *****************************************************************************
{
  g_transferProxy.ckLocalBranch()->setDestPoints( p, chare, coord, u, cb );
}

void
Transfer::addMesh( CkArrayID p, int nchare, CkCallback cb )
// *****************************************************************************
//  Register a mesh to be part of mesh-to-mesh transfer
//! \param[in] p Charm++ host proxy participating in mesh-to-mesh transfer
//! \param[in] nchare Number of mesh partitions
//! \param[in] cb Callback to continue with once finished
// *****************************************************************************
{
  auto id = static_cast<std::size_t>(CkGroupID(p).idx);
//std::cout << "addMesh: " << nchare << ", id:" << id << '\n';
  assert( m_proxyMap.count(id) == 0 );
  CkArrayOptions opts;
  opts.bindTo( p );
  opts.setNumInitial( nchare );
  MeshData mesh;
  mesh.nchare = nchare;
  mesh.firstchunk = m_current_chunk;
  mesh.proxy = CProxy_NodeSearch::ckNew( p, mesh, cb, opts );
  m_proxyMap[ id ] = mesh;
  m_current_chunk += nchare;
}

void
Transfer::setMesh( CkArrayID p, const MeshData& mesh )
// *****************************************************************************
//! \param[in] p Charm++ host proxy participating in mesh-to-mesh transfer
// *****************************************************************************
{
  m_proxyMap[static_cast<std::size_t>(CkGroupID(p).idx)] = mesh;
//std::cout << "setMesh: pe:" << CkMyPe() << ": "; for (const auto& [i,m] : m_proxyMap) std::cout << " groupid:" << i << " firstchunk:" << m.firstchunk << " nchare:" << m.nchare << ' '; std::cout << '\n';
}

void
Transfer::setSourceTets( CkArrayID p,
                         int chare,
                         const std::vector< std::size_t >& inpoel,
                         const std::array< std::vector< double >, 3 >& coord,
                         const tk::Fields& u,
                         CkCallback cb )
// *****************************************************************************
//  Configure source mesh
//! \param[in] inpoel Source mesh connectivity
//! \param[in] coord Source mesh node coordinates
//! \param[in] u Source solution data
// *****************************************************************************
{
  m_src = static_cast< std::size_t >( CkGroupID(p).idx );
  assert( m_proxyMap.count(m_src) );
  NodeSearch* w = tk::cref_find( m_proxyMap, m_src ).proxy[ chare ].ckLocal();
  assert( w );
  w->setSourceTets( inpoel, coord, u, cb );
}

void
Transfer::setDestPoints( CkArrayID p,
                         int chare,
                         const std::array< std::vector< double >, 3 >& coord,
                         tk::Fields& u,
                         CkCallback cb )
// *****************************************************************************
//  Configure destination mesh
//! \param[in] coord Pointer to the coordinate data for the destination mesh
//! \param[in,out] u Pointer to the solution data for the destination mesh
//! \param[in] cb Callback to call once this chare received all solution data
// *****************************************************************************
{
  m_dst = static_cast< std::size_t >( CkGroupID(p).idx );
  assert( m_proxyMap.count(m_dst) );
  NodeSearch* w = tk::cref_find( m_proxyMap, m_dst ).proxy[ chare ].ckLocal();
  assert( w );
  w->setDestPoints( coord, u, cb );
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
  const auto& dst = tk::cref_find( m_proxyMap, m_dst );
  auto first = static_cast< std::size_t >( dst.firstchunk );
  auto nchare = static_cast< std::size_t >( dst.nchare );
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
    const auto& src = tk::cref_find( m_proxyMap, m_src );
    dst.proxy[ static_cast<int>(i) ].processCollisions( src,
      static_cast< int >( separated[ i ].size() ),
      separated[ i ].data() );
  }
}

} // transfer::

#include "NoWarning/transfer.def.h"
