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

  // TODO: Need to make sure this is actually correct
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
               const tk::Fields& u )
// *****************************************************************************
//  API for configuring source mesh
//! \param[in] inpoel Source mesh connectivity
//! \param[in] coord Source mesh node coordinates
//! \param[in] u Source solution data
// *****************************************************************************
{
  g_transferProxy.ckLocalBranch()->setSourceTets( p, chare, inpoel, coord, u );
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

//! ...
Transfer::Transfer() : current_chunk(0) {}

void
Transfer::addMesh( CkArrayID p,
                   int nchare,
                   CkCallback cb )
// *****************************************************************************
//  Register a mesh to be part of mesh-to-mesh transfer
//! \param[in] p Charm++ host proxy participating in mesh-to-mesh transfer
//! \param[in] nchare Number of mesh partitions
//! \param[in] cb Callback to continue with once finished
// *****************************************************************************
{
  auto id = static_cast<std::size_t>(CkGroupID(p).idx);
std::cout << "addMesh: " << nchare << ", id:" << id << '\n';
  assert( proxyMap.count(id) == 0 );
  CkArrayOptions opts;
  opts.bindTo( p );
  opts.setNumInitial( nchare );
  MeshData mesh;
  mesh.nchare = nchare;
  mesh.firstchunk = current_chunk;
  mesh.proxy = CProxy_NodeSearch::ckNew( p, mesh, cb, opts );
  proxyMap[ id ] = mesh;
  current_chunk += nchare;
std::cout << "ps addMesh: " << proxyMap.size() << '\n';
}

void
Transfer::setMesh( CkArrayID p, const MeshData& mesh )
// *****************************************************************************
//! \param[in] p Charm++ host proxy participating in mesh-to-mesh transfer
// *****************************************************************************
{
  proxyMap[static_cast<std::size_t>(CkGroupID(p).idx)] = mesh;
std::cout << "ps setMesh: " << proxyMap.size() << '\n';
}

void
Transfer::setSourceTets( CkArrayID p,
                         int chare,
                         const std::vector< std::size_t >& inpoel,
                         const std::array< std::vector< double >, 3 >& coord,
                         const tk::Fields& u )
// *****************************************************************************
//  Configure source mesh
//! \param[in] inpoel Source mesh connectivity
//! \param[in] coord Source mesh node coordinates
//! \param[in] u Source solution data
// *****************************************************************************
{
  m_sourcemesh = static_cast< std::size_t >( CkGroupID(p).idx );
  NodeSearch* w = proxyMap[ m_sourcemesh ].proxy[ chare ].ckLocal();
std::cout << "ps setSourceTets: " << proxyMap.size() << '\n';
  assert( w );
  w->setSourceTets( inpoel, coord, u );
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
  m_destmesh = static_cast< std::size_t >( CkGroupID(p).idx );
  NodeSearch* w = proxyMap[ m_destmesh ].proxy[ chare ].ckLocal();
std::cout << "ps setDestPoints: " << proxyMap.size() << '\n';
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
std::cout << "ps distColl: " << proxyMap.size() << '\n';
  //CkPrintf("Collisions found: %i\n", nColl);
  auto first = static_cast< std::size_t >( proxyMap[m_destmesh].firstchunk );
  auto nchare = static_cast< std::size_t >( proxyMap[m_destmesh].nchare );
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
    proxyMap[ m_destmesh ].proxy[ static_cast<int>(i) ].processCollisions(
        proxyMap[ m_sourcemesh ].proxy,
        proxyMap[ m_sourcemesh ].nchare,
        proxyMap[ m_sourcemesh ].firstchunk,
        static_cast< int >( separated[ i ].size() ),
        separated[ i ].data() );
  }
}

} // transfer::

#include "NoWarning/transfer.def.h"
